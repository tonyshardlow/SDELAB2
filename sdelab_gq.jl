__precompile__()
#
module sdelab_gq
#
export  sde_gq_solution, sde_gq_solution_extrap, set_fcn_mf
#
import NLsolve
import ForwardDiff
import Roots
#
using gauss_quadrature
#import gauss_quadrature
#
##########################
function init(fcn, opt)
  #=
  =#
end
############################
#

type Fcnmf
  #
  drift # u -> f(u) for u (n+1)-dimensions, f(u) scalar
  drift_ext # (z,q)-> q ( f(z,r(.))) integrate out 2nd arg
  Jdrift # Jacobian of drift (auto computed)
  Jdrift_ext # (z,q) -> q (Jdrift(z,r))
  Hdrift # Hessian of drift (auto computed)
  #
  diff # u-> g(u) diffusion for u (n+1)-dimensions, g(u) scalar
  diff_ext # (z,q)-> q ( g(z,r(.))) integrate out 2nd arg
  Jdiff # Jacobian of diffuion (auto computed)
  Hdiff # Hessian of diffusion (auto computed)
  #
  mf # mean-field fn x->r(x) for x scalar, r(mf) n-dimensional
  #
  la # scipt(L)a term in notes
  #
end
#
function set_fcn_mf(drift,diff,r)
  #=
  Define mean-field sde via Fcnmf data-type
  =#
  drift_ext=(z,q)->gauss_quadrature.do_quad(y->drift([z;r(y)]),q)
  diff_ext =(z,q)->gauss_quadrature.do_quad(y->diff( [z;r(y)]),q)
  # drift (n+1)-dim to 1-dim
  Jdrift=ForwardDiff.gradient(drift) # Jacobian (n+1)-dim -> (n+1)-dim
  # only need (1,1) derivative
  drift1=u->ForwardDiff.derivative(x->drift([x;u[2:end]]),u[1])
  Jdrift_ext=(z,q)->gauss_quadrature.do_quad(y->drift1([z;r(y)]),q)
  drift2=u->ForwardDiff.derivative(x->drift1([x;u[2:end]]),u[1])
  Hdrift=drift2
  #Hdrift=ForwardDiff.hessian(drift)
  #
  Jdiff=ForwardDiff.gradient(diff) # drift (n+1)-dim -> (n+1)-dim
  Hdiff=ForwardDiff.hessian(diff)
  diff1=u->ForwardDiff.derivative(x->diff([x;u[2:end]]),u[1])
  diff2=u->ForwardDiff.derivative(x->diff1([x;u[2:end]]),u[1])
  Hdiff=diff2
  #
  dr=ForwardDiff.derivative(r)
  ddr=ForwardDiff.derivative(dr)
  #
  la= (x,mf)->(dr(x)*drift([x;mf])
              + 0.5*ddr(x)*diff([x;mf])^2)
  #
  return Fcnmf(drift,drift_ext,Jdrift,Jdrift_ext,Hdrift,
               diff, diff_ext, Jdiff, Hdiff,
               r,la)
end
#######################
function sde_gq_solution(fcn::Fcnmf, tspan, quad0::QUAD, params)
  #
  # Apply QG integration.
  #
  # Inputs:
  # fcn   = definition of MF-SDE
  # tspan = vector of times to report
  # quad0 = initial quadrature
  # params= dict of parameters
  #
  # Outputs: quadrature and
  # stats=see below (to analyse behaviour of method)
  #
  T=maximum(tspan) # get final time
  noSteps=round(Int32,T/params["MaxStepSize"])
  t=(0:noSteps)/noSteps # grid of time points
  dt=t[2]-t[1] # fixed timestep
  #
  quad=deepcopy(quad0)
  # initialise
  solver=params["Solver"]
  if solver=="MFEM0"
    params["order"]=1
    quad_dw=gauss_quadrature.Gauss_Hermite([0,dt], 2)
    get= q->get_1st(fcn, q, dt, quad_dw, params)
  elseif solver=="MFEM0e"
    params["order"]=1
    quad_dw=gauss_quadrature.Gauss_Hermite([0,dt], 2)
    get= q->get_1st_ext(fcn, q, dt, quad_dw, params)
  elseif solver=="MFEM2"
    params["order"]=2
    quad_dw=gauss_quadrature.Gauss_Hermite([0,dt], 3)
    dtmat=dt_order2(dt,quad_dw)
    get= q->get_2nd(fcn, q, dtmat, quad_dw, params)
  elseif solver=="MFEM"
    params["order"]=1
    quad_dw=gauss_quadrature.Gauss_Hermite([0,dt], 2)
    function nonlin!(un_,up,out,r)
      out[:] =(un_-up-params["Alpha"]*fcn.drift(vec([un_[1];r]))*dt) # must avoid recreating out!
    end
    function jnonlin!(un_,out,r)
      tmp=fcn.Jdrift(vec([un_[1];r]) )[1]
      out[:,:]=1-params["Alpha"]*tmp
    end
    get= q->get_1st_imp(fcn, q, dt, quad_dw, params,nonlin!,jnonlin!)
  elseif solver=="MFEMe"
    params["order"]=1
    quad_dw=gauss_quadrature.Gauss_Hermite([0,dt], 2)
    function nonlin!(un_,up,out,q::QUAD)
      out[:] =(un_-up-params["Alpha"]*fcn.drift_ext(un_[1],q)*dt) # must avoid recreating out!
    end
    function jnonlin!(un_,out,q::QUAD)
      tmp=fcn.Jdrift_ext(un_[1],q)
      out[:,:]=1-params["Alpha"]*tmp
    end
    get= q->get_1st_imp_ext(fcn, q, dt, quad_dw, params,nonlin!,jnonlin!)
  end
  # initialise
  no_pts=zeros(Int32,noSteps+1); X_min=zeros(noSteps+1); X_max=zeros(noSteps+1)
  X_min[1]=minimum(quad.X); X_max[1]=maximum(quad.X); no_pts[1]=quad.length
  m=quad.length
  factor=params["order"]^3
  for i=1:noSteps
    quad=get(quad)
    # update degree and recompute Gauss rule/threshold
    if quad.length>factor*m && i<noSteps-2
      m=no_Gauss_points(m, dt, T-t[i+1], params)
      if params["periodic"]==false # non-periodic
        centre_mass=gauss_quadrature.do_quad(x->x,quad)
        r=abs(log(dt))*params["cutoff"]
        cutoff= [centre_mass-r,centre_mass+r]
        if params["no_partition"]>1
          gauss_quadrature.partition_simplify!(quad, m, params["no_partition"], cutoff)
        else
          gauss_quadrature.simplify!(quad, m, cutoff)
        end
      else # periodic
        if params["no_partition"]>1
          gauss_quadrature.partition_periodic_simplify!(quad, m, params["no_partition"], params["length_scale"])
        else
          gauss_quadrature.periodic_simplify!(quad, m, params["length_scale"])
        end
      end # end if periodic
    end # reduce quad.length
    # stats
    X_min[i+1]=minimum(quad.X); X_max[i+1]=maximum(quad.X);
    no_pts[i+1]=quad.length;
  end # time loop
  gauss_quadrature.sortquad!(quad)
  stats=Dict()
  stats["min"]=X_min;  stats["max"]=X_max
  stats["no_pts"]=no_pts;  stats["t"]=t
  return quad, stats
end
##########################
function sde_gq_solution_extrap(fcn::Fcnmf, tspan, quad0::QUAD, params)
  # Computes the Talay-Tubaro extrapolation
  # Returns a rule with negative weights
  # Only effective for first-order integrator
  #
  # Same calling syntax as sde_pdf_solution
  quad,stats=sde_gq_solution(fcn, tspan, quad0, params)
  params1=deepcopy(params)
  params1["MaxStepSize"]=params["MaxStepSize"]/2
  quad1,stats1=sde_gq_solution(fcn, tspan, quad0, params1)
  x=[quad1.X; quad.X]; w=[2*quad1.w; -quad.w]
  new_quad=QUAD(length(x), x ,w)
  gauss_quadrature.sortquad!(new_quad)
  return new_quad, stats1
end
#############################################
function get_1st_ext(fcnmf::Fcnmf, quad_in::QUAD, dt, quad_dw::QUAD, params)
  # Inputs: fcnmf =  drift and diffusion functions
  # quad_in=current quadrature approx
  # dt=time step
  #
  # Outputs: quad_out the next timestep quad approx by Euler-Maruyama (explicit)
  #
  # same as get_1st, except mean-field calculation is completed OUTSIDE the function eval
  Y=zeros(quad_in.length,quad_dw.length)
  W=deepcopy(Y)
  drift1=z->fcnmf.drift_ext(z,quad_in)
  diff1 =z->fcnmf.diff_ext(z,quad_in)
  for j=1:quad_in.length
    x=quad_in.X[j]
    xpfx_dt = x + drift1(x)*dt
    gx      =     diff1(x)
    #
    Y[j,:]=xpfx_dt +   quad_dw.X * gx
    W[j,:]=quad_in.w[j]*quad_dw.w
  end
  if params["periodic"]==true # periodic
     Y=mod(Y,params["length_scale"])
   end
  quad_out=gauss_quadrature.QUAD(length(Y), vec(Y), vec(W) )
  return quad_out
end

#############################################
function get_1st(fcnmf::Fcnmf, quad_in::QUAD, dt, quad_dw::QUAD, params)
  # Inputs: fcnmf =  drift and diffusion functions
  # quad_in=current quadrature approx
  # dt=time step
  #
  # Outputs: quad_out the next timestep quad approx by Euler-Maruyama (explicit)
  #
  Y=zeros(quad_in.length,quad_dw.length)
  W=deepcopy(Y)
  mf=gauss_quadrature.do_quad(y->fcnmf.mf(y),quad_in)
  for j=1:quad_in.length
    x=quad_in.X[j]
    xpfx_dt = x + fcnmf.drift([x;mf])*dt
    gx=fcnmf.diff([x;mf])
    #
    Y[j,:]=xpfx_dt +   quad_dw.X * gx
    W[j,:]=quad_in.w[j]*quad_dw.w
  end
  if params["periodic"]==true # periodic
     Y=mod(Y,params["length_scale"])
   end
  quad_out=gauss_quadrature.QUAD(length(Y), vec(Y), vec(W) )
  return quad_out
end
##############################################

function get_2nd(fcnmf::Fcnmf, quad_in::QUAD, dtmat, quad_dw::QUAD, params)
  # Inputs: fcn =  drift and diffusion functions
  # quad_in=current quadrature approx by 2nd order
  #
  Y=zeros(quad_in.length, quad_dw.length)
  W=deepcopy(Y)
  mf=gauss_quadrature.do_quad(fcnmf.mf, quad_in)
  La=gauss_quadrature.do_quad(x->fcnmf.la(x,mf), quad_in)
  #
  for j=1:quad_in.length
    Y[j,:]=quad_in.X[j]+ coeffs_order2(quad_in.X[j], mf, La, fcnmf)*dtmat
    W[j,:]=quad_in.w[j]*quad_dw.w
  end
  #
  if params["periodic"]==true
     Y=mod(Y,params["length_scale"])
   end
  quad_out=gauss_quadrature.QUAD(length(Y), vec(Y), vec(W) )
  return quad_out
end
##########################################
#(fcn, q, dt, quad_dw, params,nonlin!,jnonlin!))
function get_1st_imp(fcnmf::Fcnmf, quad_in::QUAD, dt, quad_dw::QUAD, params, nonlin!, jnonlin!)
  # Inputs: fcn =  drift and diffusion functions
  # quad_in=current quadrature approx
  # dt=time step
  #
  # Outputs: quad_out the next timestep quad approx
  P=quad_in.length
  sqrt_dt=sqrt(dt)
  Y=zeros(P,1)
  Z=zeros(P,1)
  mf=gauss_quadrature.do_quad(fcnmf.mf, quad_in)
  helper=(y1,y2)->NLsolve.nlsolve(
    (u_,out)-> nonlin!(u_,y1, out,mf),                       (u_,out)->jnonlin!(u_,    out,mf),
    [y2])
  for j=1:P
    x=quad_in.X[j]
    xpfx_dt  = x + (1-params["Alpha"])*fcnmf.drift([x;mf])*dt
    gx_sqrth = fcnmf.diff([x;mf])*sqrt_dt
    #
    un=xpfx_dt+gx_sqrth
    solve=helper(un,x)
    Y[j]=solve.zero[1]
    #
    un=xpfx_dt-gx_sqrth

    solve=helper(un,x)
    Z[j]=solve.zero[1]
  end
  quad_out=QUAD(2*P, vec([Y Z]), [quad_in.w;quad_in.w]/2)
  return quad_out
end
##########################################
#(fcn, q, dt, quad_dw, params,nonlin!,jnonlin!))
function get_1st_imp_ext(fcnmf::Fcnmf, quad_in::QUAD, dt, quad_dw::QUAD, params, nonlin!, jnonlin!)
  # Inputs: fcn =  drift and diffusion functions
  # quad_in=current quadrature approx
  # dt=time step
  #
  # Outputs: quad_out the next timestep quad approx
  P=quad_in.length
  sqrt_dt=sqrt(dt)
  Y=zeros(P,1)
  Z=zeros(P,1)
  drift1=z->fcnmf.drift_ext(z,quad_in)
  diff1 =z->fcnmf.diff_ext(z, quad_in)
  helper=(y1,y2)->NLsolve.nlsolve(
    (u_,out)-> nonlin!(u_,y1, out,quad_in),                       (u_,out)->jnonlin!(u_,    out,quad_in),
    [y2])
  for j=1:P
    x=quad_in.X[j]
    xpfx_dt  = x + (1-params["Alpha"])*drift1(x)*dt
    gx_sqrth = diff1(x)*sqrt_dt
    #
    un=xpfx_dt+gx_sqrth
    solve=helper(un,x)
    Y[j]=solve.zero[1]
    #
    un=xpfx_dt-gx_sqrth
    solve=helper(un,x)
    Z[j]=solve.zero[1]
  end
  quad_out=QUAD(2*P, vec([Y Z]), [quad_in.w;quad_in.w]/2)
  return quad_out
end
# #####################################################
# function getYZA_implicit(mf::Fcnmf, quad_in::QUAD, dt)
#   # Inputs: fcn =  drift and diffusion functions
#   # quad_in=current quadrature approx
#   # dt=time step
#   #
#   # Outputs: quad_out the next timestep quad approx
#   P=quad_in.length
#   sqrt_dt=sqrt(dt)
#   Y=zeros(P,1)
#   Z=zeros(P,1)
#   A=zeros(P,1)
#   for j=1:P
#     x=quad_in.X[j]
#     xpfx_dt  = x + (1-mf.Alpha)*do_quad(y->mf.f([x,y]),quad_in)*dt
#     gx_sqrth = sqrt(3)*do_quad(y->mf.g([x,y],sqrt_dt),quad_in)
#     #
#     un=xpfx_dt+gx_sqrth
#     solve=NLsolve.nlsolve((u_,out)->mf.nonlin!(u_,un,out,dt,quad_in),
#                           (u_,out)->mf.jnonlin!(u_,out,dt,quad_in),
#                           [x])
#     #
#     Y[j]=solve.zero[1]
#     #
#     un=xpfx_dt-gx_sqrth
#     solve=NLsolve.nlsolve((u_,out)->mf.nonlin!(u_,un,out,dt,quad_in),
#                           (u_,out)->mf.jnonlin!(u_,out,dt,quad_in),
#                           [x])
#     #
#     Z[j]=solve.zero[1]
#     #
#     un=xpfx_dt
#     solve=NLsolve.nlsolve((u_,out)->mf.nonlin!(u_,un,out,dt,quad_in),
#                           (u_,out)->mf.jnonlin!(u_,out,dt,quad_in),
#                           [x])
#     #
#     A[j]=solve.zero[1]
#   end
#   quad_out=QUAD(3*P, vec([Y Z A]),
#               [quad_in.w/6;quad_in.w/6;quad_in.w*(2/3)])
#   return quad_out
# end
######################################################
function no_Gauss_points(m_old,dt,t_remaining,params)
  # function m=no_Gauss_points(m_old, dt, t_remaining, params)
  # return number of Gauss points (at least m_old),
  # to preserve dt^(order) accuracy
  #
  m_new=m_old
  if t_remaining>0
    ldt=abs(log(dt))
    order=params["order"]
    if params["periodic"]==true
      lldt_term=params["length_scale"]/params["no_partition"]
    else
      lldt_term=log((8*params["length_scale"]^2)*ldt)
    end
    lt_r=(0.0)*abs(log(t_remaining))
    if params["skiphalf"]==true
      order1=order+0.5
    else
      order1=order+1
    end
    constant=order1*ldt - 0.5*lldt_term -(order-1)*lt_r+params["const"]
    linear=lldt_term+lt_r
    g=y->f(y,constant,linear)
    m_large=1e4
    low  =g(m_old)
    upper=g(m_large)
    if (low*upper<0)
      m_new=Roots.fzero(g, [m_old, m_large])
    end
  end
  return round(Int16,m_new)
end
#####################################################
function f(m, constant, linear)
  # helper function for no_Gauss_points (roots gives a value for number of Gauss points)
  # Inputs: m=number of Gauss points
  #
  y= constant + m * linear - lgamma(2*m-1)
  return y
end
#####################################################
function coeffs_order2(x,mf,La,fcn::Fcnmf)
  # Coefficients for second-order method
  #
  a    = fcn.drift([x;mf])
  grada=fcn.Jdrift([x;mf])
  b     = fcn.diff([x;mf])
  gradb =fcn.Jdiff([x;mf])
  #
  Lan=[a;La]
  #
  term1=gradb[1]*b
  term2=grada[1]*b +dot(gradb,Lan) + 0.5*fcn.Hdiff([x;mf]) * b^2
  term3=dot(grada,Lan) + 0.5 * fcn.Hdrift([x;mf])*b^2
  #
  return [a  b  term1  term2  term3] # row vector
end
#
function dt_order2(dt,quad_dw::QUAD)
  # multiplier for second-order method
  # want to generate increment by coeffs_order2 * dtmat
  dtmat=zeros(5,quad_dw.length)
  dtmat[1,:]=dt
  dtmat[2,:]=quad_dw.X
  dtmat[3,:]=0.5* (quad_dw.X.^2 - dt)
  dtmat[4,:]=0.5* quad_dw.X * dt
  dtmat[5,:]=0.5* dt^2
  return dtmat
end
#
end # module definition

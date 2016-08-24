__precompile__()

####################
module SDELab2
import NLsolve
import ForwardDiff
import iteratedintegral

export set_opt, sde_strong_solution, set_fcn,  method_code, do_plotting

method_code=Dict()
method_code["SIE0"]="Strong Ito Euler "
method_code["SSEH0"]="Strong Stratonovich Euler-Heun "
method_code["SIE1"]="Strong Ito Euler (drift implicit)"
method_code["SSEH1"]="Strong Stratonovich Euler (drift implicit)"
method_code["SIE"]="Strong Ito Euler (alpha drift-implicit)"
method_code["SSEH"]="Strong Stratonovich Euler (alpha drift-implicit)"
method_code["SIM0"]="Strong Ito Milstein explicit"
method_code["SSM0"]="Strong Stratonovich Milstein explicit"
method_code["SIM1"]="Strong Ito Milstein (drift implicit)"
method_code["SSM1"]="Strong Stratonovich Milstein (drift implicit)"
method_code["SIM"]="Strong Ito Milstein (alpha drift-implicit)"
method_code["SSM"]="Strong Stratonovich Milstein (alpha drift-implicit)"
method_code["SIBDF"]="Strong Ito BDF-2"
# sdelab_gq
method_code["MFEM0"]="Mean-field SDE Gauss-quadrature Euler-Maruyama"
method_code["MFEM2"]="Mean-field SDE Gauss-quadrature second-order explicit"
method_code["MFEM"] ="Mean-field SDE Gauss-quadrature alpha-drift implicit (first order)"
method_code["MFEM0e"]="Mean-field SDE Gauss-quadrature Euler-Maruyama (exterior mean-field evaluation)"
method_code["MFEMe"] ="Mean-field SDE Gauss-quadrature alpha-drift implicit (first order, exterior mean-field)"
using sdelab_gq
using gauss_quadrature
#####
function __init__()
end
##########################################
function get_increment!(x)
  x.dW[:]=sqrt(x.dt)*randn(x.m)
end
##########################################
# explicit Euler methods (i.e., Euler-Maruyama and Heun)
# data structure to store required bits
# and one-step routines
type EM
  m
  drift
  diff_mat
  dt
  dW
end
# Ito version
function one_step_em(u,x::EM)
  get_increment!(x)
  un=u+x.drift(u)*x.dt+x.diff_mat(u,x.dW)
  return un
end
# Stratonovich (Heun method)
function one_step_heun(u,x::EM)
  get_increment!(x)
  uaux=u+x.diff_mat(u,x.dW)
  un=u+x.drift(u)*x.dt+0.5*(x.diff_mat(u,x.dW)+x.diff_mat(uaux,x.dW))
  return un
end
##################################################
# implicit Euler methods
# parameters for nonlinear solver
type NLSOLVE_param
  maxfeval
  ftol
  xtol
end
# required bits for implicit-type Euler methods
type EMImp
  m
  drift
  diff_mat
  Alpha
  dt
  dW
  nonlin
  jnonlin
  nlsolve::NLSOLVE_param
end
# Ito version of onestep
function one_step_em_alpha(u,x::EMImp)
  get_increment!(x)
  up=u + (1-x.Alpha)*x.drift(u)*x.dt + x.diff_mat(u,x.dW)
  solve=NLsolve.nlsolve((u,out)->x.nonlin(u,up,out),x.jnonlin,u, iterations=x.nlsolve.maxfeval, ftol=x.nlsolve.ftol, xtol=x.nlsolve.xtol)
  return solve.zero
end
# Stratonovich version of onestep
function one_step_heun_alpha(u,x::EMImp)
  get_increment!(x)
  uaux=u+x.diff_mat(u,x.dW)
  up=u+ (1-x.Alpha)*x.drift(u)*x.dt+0.5*(x.diff_mat(u,x.dW)+x.diff_mat(uaux,x.dW))
  solve=NLsolve.nlsolve((u,out)->x.nonlin(u,up,out),x.jnonlin,u)
  return solve.zero
end
# BDF method
# required bits
type BDF
  m
  drift
  diff_mat
  dt
  dW
  dWp
  up
  nonlin
  jnonlin
  nlsolve::NLSOLVE_param
end
# Ito BDF
function one_step_bdf(u,x::BDF)
  get_increment!(x)
  up=(4/3)*u -(1/3)*x.up +x.diff_mat(u,x.dW)-(1/3)* x.diff_mat(x.up,x.dWp)
  solve=NLsolve.nlsolve((u,out)->x.nonlin(u,up,out),x.jnonlin,u, iterations=x.nlsolve.maxfeval, ftol=x.nlsolve.ftol, xtol=x.nlsolve.xtol)
  x.up[:]=u
  x.dWp[:]=x.dW
  return solve.zero
end
##########################################################
# Explicit Milstein methods
type Mil
  d
  m
  drift
  diff_mat
  Jdiff_vecs
  S
  dt
end
# Milstein explicit, Stratonovich version
function one_step_strat_mil(u,x::Mil)
  IteratedIntegral.get_strat!(x.dt,x.S)
  J=x.S["J"];  K=x.S["K"]
  un=u+x.drift(u)*x.dt+x.diff_mat(u,J)
  for i=1:x.m
    Ki=reshape(K[i,:],(x.m,1))
    un+=x.Jdiff_vecs[i](u)*x.diff_mat(u,Ki)
  end
  return reshape(un,x.d)
end
# Milstein explicit, Ito version
function one_step_ito_mil(u,x::Mil)
  IteratedIntegral.get_ito!(x.dt,x.S)
  J=x.S["J"];  K=x.S["K"]
  un=u+x.drift(u)*x.dt+x.diff_mat(u,J)
  for i=1:x.m
    Ki=reshape(K[i,:],(x.m,1))
    un+=x.Jdiff_vecs[i](u)*x.diff_mat(u,Ki)
  end
  return reshape(un,x.d )
end
###################################################
# Implicit Milstein methods
type MilImp
  d
  m
  drift
  diff_mat
  Jdiff_vecs
  Alpha
  S
  dt
  nonlin
  jnonlin
  nlsolve::NLSOLVE_param
end
# Milstein implicit Stratonovich
function one_step_strat_mil_alpha(u,x::MilImp)
  IteratedIntegral.get_strat!(x.dt,x.S)
  J=x.S["J"];   K=x.S["K"]
  up=u + x.Alpha*x.drift(u)*x.dt + x.diff_mat(u,J)
  for i=1:x.m
    Ki=reshape(K[i,:],(x.m,1))
    up+=x.Jdiff_vecs[i](u)*x.diff_mat(u,Ki)
  end
  solve=NLsolve.nlsolve((u,out)->x.nonlin(u,up,out),x.jnonlin,u,iterations=x.nlsolve.maxfeval, ftol=x.nlsolve.ftol, xtol=x.nlsolve.xtol)
  return solve.zero
end
# Milstein implicit Ito
function one_step_ito_mil_alpha(u,x::MilImp)
  IteratedIntegral.get_ito!(x.dt,x.S)
  J=x.S["J"];   K=x.S["K"]
  up=u+ x.Alpha*x.drift(u)*x.dt +x.diff_mat(u,J)
  for i=1:x.m
    Ki=reshape(K[i,:],(x.m,1))
    up += x.Jdiff_vecs[i](u)*x.diff_mat(u,Ki)
  end
  solve=NLsolve.nlsolve((u,out)->x.nonlin(u,up,out), x.jnonlin, u, iterations=x.nlsolve.maxfeval, ftol=x.nlsolve.ftol, xtol=x.nlsolve.xtol)
  return solve.zero
end
###################################################
# Set solver options (many defaults)
function set_opt(maxStepSize,solver)
  opt=Dict()
  opt["MaxStepSize"]=maxStepSize
  opt["Solver"]=solver
  opt["Alpha"]=0.5
  opt["nlsolve_param"]=NLSOLVE_param(10,1e-6,1e-8)
  opt["max_p"]=1000
  opt["calc_final_W"]=true
  opt["is_ito"]=(solver[1:2]=="SI")
  # plotting
  opt["always_clear"]=true
  opt["do_plotting"]=false
  opt["output_sel"]=[1,2]
  return opt
end
#
#
function sde_strong_solution(fcn, tspan, y0, opt)
  #=
  many driver routine:
  fcn= SDE definition
  tspan= vector of required times
  y0=value of sample path at first time in tspan
  opt=solver options
  =#
  ts=minimum(tspan); tf=maximum(tspan)
  N=round((tf-ts)/opt["MaxStepSize"]); dt=(tf-ts)/N
  solver=opt["Solver"]
  t_out=Array(Float64,length(tspan))
  y_out=Array(Float64,(length(tspan),fcn.d))
  y=deepcopy(y0); y_out[1,:]=y # initial data
  t=ts; t_out[1]=ts # initial time
  #
  Wf=zeros(fcn.m)
  #
  if solver=="SIE0" || solver=="SSEH0"# StrongItoEuler Alpha=0 explicit
    x=EM(fcn.m,fcn.drift,fcn.diff_mat,dt,zeros(fcn.m))
    dW=x.dW
    if solver=="SIE0"
      call_onestep(y)=(one_step_em(y,x))
    else
      call_onestep(y)=(one_step_heun(y,x))
    end
  elseif solver=="SIE1" || solver=="SSEH1"# strongItoEuler Alpha=1 implicit
    nlsolve_param=opt["nlsolve_param"]
    function nonlin!(un_,up,out)
      out[:] =(un_-up-fcn.drift(un_)*dt) # must avoid recreating out!
    end
    function jnonlin!(un_,out)
      out[:,:]=eye(fcn.d)-fcn.Jdrift(un_)*dt
    end
    x=EMImp(fcn.m,fcn.drift,fcn.diff_mat,1.,dt,zeros(fcn.m),nonlin!,jnonlin!,nlsolve_param)
    dW=x.dW
    if solver=="SIE1"
      call_onestep(y)=(one_step_em_alpha(y,x))
    else
      call_onestep(y)=(one_step_heun_alpha(y,x))
    end
  elseif solver=="SIE" || solver =="SSEH"# strongItoEuler Alpha parameter
    Alpha=opt["Alpha"]
    nlsolve_param=opt["nlsolve_param"]
    function nonlin!(un_,up,out)
      out[:] =(un_-up-Alpha*fcn.drift(un_)*dt) # must avoid recreating out!
    end
    function jnonlin!(un_,out)
      out[:,:]=eye(fcn.d)-Alpha*fcn.Jdrift(un_)*dt
    end
    x=EMImp(fcn.m,fcn.drift,fcn.diff_mat,Alpha,dt,zeros(fcn.m),nonlin!,jnonlin!, nlsolve_param)
    dW=x.dW
    if solver=="SIE"
      call_onestep(y)=(one_step_em_alpha(y,x))
    else
      call_onestep(y)=(one_step_heun_alpha(y,x))
    end
  elseif solver=="SIM0" || solver=="SSM0"# StrongItoMilstein Alpha=0 explicit
    S=IteratedIntegral.init(fcn,opt)
    x=Mil(fcn.d,fcn.m,fcn.drift,fcn.diff_mat,fcn.Jdiff_vecs,S,dt)
    dW=x.S["J"]
    if solver=="SIM0"
      call_onestep(y)=(one_step_ito_mil(y,x))
    else
      call_onestep(y)=(one_step_strat_mil(y,x))
    end
  elseif solver=="SIM1" || solver=="SSM1"# StrongItoMilstein Alpha=1 implicit
      S=IteratedIntegral.init(fcn,opt)
      nlsolve_param=opt["nlsolve_param"]
      function nonlin!(un_,up,out)
        out[:] =(un_-up-fcn.drift(un_)*dt) # must avoid recreating out!
      end
      function jnonlin!(un_,out)
        out[:,:]=eye(fcn.d)-fcn.Jdrift(un_)*dt
      end
      x=MilImp(fcn.d,fcn.m,fcn.drift,fcn.diff_mat,fcn.Jdiff_vecs,0.,S,dt,nonlin!,jnonlin!,nlsolve_param)
      dW=x.S["J"]
      if solver=="SIM1"
        call_onestep(y)=(one_step_ito_mil_alpha(y,x))
      else
        call_onestep(y)=(one_step_strat_mil_alpha(y,x))
      end
    elseif solver=="SIM" || solver=="SSM"
      S=IteratedIntegral.init(fcn,opt)
      nlsolve_param=opt["nlsolve_param"]
      Alpha=opt["Alpha"]
      function nonlin!(un_,up,out)
        out[:] =(un_-up-Alpha*fcn.drift(un_)*dt) # must avoid recreating out!
      end
      function jnonlin!(un_,out)
        out[:,:]=eye(fcn.d)-Alpha*fcn.Jdrift(un_)*dt
      end
      x=MilImp(fcn.d,fcn.m,fcn.drift,fcn.diff_mat,fcn.Jdiff_vecs,Alpha,S,dt,nonlin!,jnonlin!,nlsolve_param)
      dW=x.S["J"]
      if solver=="SIM"
        call_onestep(y)=(one_step_ito_mil_alpha(y,x))
      else
        call_onestep(y)=(one_step_strat_mil_alpha(y,x))
      end
  elseif solver=="SIBDF"
    nlsolve_param=opt["nlsolve_param"]
    # starting method
    function nonlin0!(un_,up,out)
      out[:] =(un_-up-fcn.drift(un_)*dt) # must avoid recreating out!
    end
    function jnonlin0!(un_,out)
      out[:,:]=eye(fcn.d)-fcn.Jdrift(un_)*dt
    end
    x0=EMImp(fcn.m,fcn.drift,fcn.diff_mat,1.,dt,zeros(fcn.m),nonlin0!,jnonlin0!,nlsolve_param)
    y1=one_step_em_alpha(y,x0)
    t=t+dt
    Wf+=x0.dW
    # bdf method
    function nonlin!(un_,up,out)
      out[:] =(un_-up-(2/3)*fcn.drift(un_)*dt) # must avoid recreating out!
    end
    function jnonlin!(un_,out)
      out[:,:]=eye(fcn.d)-(2/3)*fcn.Jdrift(un_)*dt
    end
    x=BDF(fcn.m,fcn.drift,fcn.diff_mat,dt,zeros(fcn.m),x0.dW,deepcopy(y0),nonlin!,jnonlin!,nlsolve_param)
    y=y1;
    dW=x.dW
    call_onestep(y)=(one_step_bdf(y,x))
  end

  #
  for j=2:(length(tspan)) # loop over required times
    N=round(Int,(tspan[j]-t)/dt) # steps to next tspan node
    for i=1:N # does not run if N=0
      t=t+dt
      y=call_onestep(y)
      if opt["calc_final_W"]==true
        Wf+=dW
      end
    end
    t_out[j]=t
    y_out[j,:]=y # constant approx if tspan too fine relative to dt
  end
  try
    do_plotting(t_out,y_out,opt)
  catch
    print("plotting failed")
  end
  return t_out,y_out,Wf
end
#
using PyPlot
function do_plotting(t,y,opt)
  if opt["do_plotting"]
    if opt["always_clear"]
      clf()
    end
    str=opt["output_plot_type"]
    if str=="path_plot"
      sel=opt["output_sel"]
      plot(t,y[:,sel])
    elseif str=="phase_plot"
      sel=opt["output_sel"]
      plot(y[:,sel[1]], y[:,sel[2]])
    elseif str=="time_phase"
      print("not yet implemented :-(")
    elseif str=="phase3_plot"
      print("not yet implemented :-(")
    end
  end
end
#
type Fcn
  d # state-space dimension
  m # no. Brownian motions
  drift # u -> f(u)
  Jdrift # Jacobian of drift (auto computed)
  diff_mat # (u,dw)-> g(u)*dw
  diff_vecs # [i]: u -> g_i(u) for i=1:m
  Jdiff_vecs # Jacobian of diffusion vectors (auto computed)
  noise_type # 1=unstructure; 2=diagonal; 3=commutative
end
#
function set_fcn(d,m,drift,diff_mat,diff_vecs,noise_type=1)
  #=
  Define sde via Fcn type
  =#
  Jdrift=ForwardDiff.jacobian(drift)
  if diff_vecs!=Union{}
    Jdiff_vecs=Array(Function,m)
    for i=1:m
      Jdiff_vecs[i]=ForwardDiff.jacobian(diff_vecs[i])
    end
  else
    Jdiff_vecs=Union{}
  end
  fcn=Fcn(d,m,drift,Jdrift,diff_mat,diff_vecs,Jdiff_vecs,noise_type)
end
#
end #sdelab2
#

module IteratedIntegral
export init, get_ito, get_strat
#################################################################
function init(fcn, opt)
  #=
  fcn=SDE definition (for dimension m and noise type)
  opt= options (for max_p)
  =#
  S=Dict("KL_C"=>1, "m"=>fcn.m, "NoiseType"=>fcn.noise_type)
  try
    S["max_p"]=opt["max_p"]
  catch
    S["max_p"]=1000
  end
  # get tail of sum_{p=P]^infty 1/p^2
  no_terms=100
  H2_tail=Array(Float16,no_terms)
  H2_tail[1]=pi^2/6 # zeta(2)
  for i=2:no_terms
    H2_tail[i]=H2_tail[i-1]-1/(i^2)
  end
  #S["H2_tail"]=H2_tail
  # map for upper-triangular to vector
  T=Int16
  m=fcn.m
  m2=round(T,m*(m-1)/2)
  ri=Array(T,m2)
  ci=Array(T,m2)
  k=1
  for i = 1:(m-1)
      for j = (i+1):m
        ri[k]=i
        ci[k]=j
        k=k+1
      end
  end
  S=merge(S,Dict("ri"=>ri, "ci"=>ci, "m2"=>m2))
  # initialise storage
  S["tmp_m2"]=Array(Float64,(m,m))
  S["tail"]=Array(Float64,(m,m))
  S["J"]=Array(Float64,m)
  S["K"]=Array(Float64,(m,m))
  #
  return S
  end
######################
function get_p(dt,S)
  #=
  get number of terms for a KL expansion in terms of dt
  and KL_C
  (not actually used as tail correction is implemented here)
  =#
  thresh=2*pi^2*S["KL_C"]^2 * dt
  p=findfirst(x->x<thresh,S["H2_tail"])
  if p<S["max_p"]
      p= 1 / (2 * (pi* S["KL_C"])^2 * dt)
  end
  return round(Int,p)
end
####################################
function Gamma_1!(out,b,S)
  #=
  helper function: see notes \Gamma_1
  in: b m2-vector
      S dictionary (see init function for defs of m and m2 and ri,ci)
  out: out (m x m)-matrix
  =#
  r=S["ri"]
  c=S["ci"]
  for i in 1:length(b)
    out[r[i],c[i]]=b[i]
    out[c[i],r[i]]=-b[i]
  end
end
######################################
function Gamma_2!(out,A,dW,S)
  #=
  helper function: see notes
  in: A (m x m)-matrix
  dW m-vector
  S dictionary (see init function)
  out: out (m x m)-matrix
  =#
  for i=1:length(dW)
    for j=1:length(dW)
      out[i,j]=vecdot(A[i,:],dW)*dW[j]
    end
  end
end
########################################
function Gamma_3!(out,A,S)
  #=
  helper function: see notes
  in: A (m x m)-matrix
  S dictionary (see init function)
  out: m2-vector
  =#
  r=S["ri"]
  c=S["ci"]
  for i in length(r)
    out[i]=A[r[i],c[i]] - A[c[i],r[i]]
  end
end
####################################
function compute_tail!(tail,dt,dW,ap,S)
  #=
  Tail approximation for KL expansion of stochastic integrals
  in: dt timestep
  dW m-vector
  ap coefficient (see notes)
  S dictionary (see init function)
  out: tail (m x m)-matrix
  =#
  scale=dt^2*ap/(2*pi^2)
  alpha=sqrt(1+sum(dW.*dW))
  noiseVec=randn(S["m2"])*sqrt(scale)
  tmp_M=noiseVec/((1+alpha)*dt)
  tmp_m2=S["tmp_m2"]
  Gamma_1!(tmp_m2,tmp_M,S)
  Gamma_2!(tail,tmp_m2,dW,S)
  Gamma_3!(tmp_M,tail,S)
  Gamma_1!(tail,tmp_M+noiseVec,S)
end
#
function get_strat!(dt,S)
  #=
  Main routine for generating integrals (Stratonvich)
  in: dt timestep
  S (see init function)
  out: J,K updated in S
  =#
  m=S["m"]
  J=S["J"]
  K=S["K"]
  if S["NoiseType"]==2 # diagonal
    J[:]=sqrt(dt)*randn(m)
    fill!(K,0)
    for i=1:m
      K[i,i]=0.5*J[i]^2
    end
  elseif S["NoiseType"]==3 # commutative
    J[:]=sqrt(dt)*randn(m)
    for i=1:m
      for j=1:m
       K[i,j]=0.5*J[i]*J[j]
    end
  end
  else # unstructured
    BBT!(dt,S)
  end
end
######################################
function get_ito!(dt,S)
  #=
  Main routine for generating integrals (Ito)
  in: dt timestep
  S (see init function)
  out: J,K updated in S
  =#
  get_strat!(dt,S)
  K=S["K"]
  for i=1:S["m"]
    K[i,i]-=0.5*dt
  end
end
########################################
function BBT!(dt,S)
  #=
  Routine for generating integrals (Stratonovich) unstructured case
  in: dt timestpe
  S dictionary (see init function)
  out: J,K updated in S
  =#
  m=S["m"]
  J=S["J"]
  K=S["K"]
  #
  p_const=sqrt(m*(m-1)/(24*dt)/(pi*S["KL_C"]))
  a_p=pi^2/6
  J[:]=sqrt(dt)*randn(m)
  for i=1:m
    K[i,i]=0.5*J[i]^2
  end
  p=round(p_const*sqrt(m*4*sum(J.*J)/dt))
  for r=1:(p+1)
    a_p=a_p-(1/r^2)
    sd_sq=dt/(2*pi*r)
    sd=sqrt(sd_sq)
    eta=randn(m)*sd
    zeta=randn(m)*sd
    K[:,1]-= zeta*(2/sqrt(pi*r))
    for i=2:m
      for j=1:i
        K[i,j]+=zeta[i]*eta[j]-eta[i]*zeta[j]
      end
    end
  end
  for i=2:m
    for j=1:i
      tmp1=K[i,j]-0.5*(K[j,1]*J[i]-K[i,1]*J[j])
      tmp2=0.5*J[i]*J[j]
      K[i,j]=tmp2+tmp1
      K[j,i]=tmp2-tmp1
    end
  end
  tail=S["tail"]
  compute_tail!(tail,dt,J,a_p,S)
  K+=tail
end
end
####################
module SDELab2
import NLsolve
import ForwardDiff
import IteratedIntegral
export set_opt, sde_strong_solution, set_fcn, method_code, do_plotting
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
        out[:] =(un_-up-(1-Alpha)*fcn.drift(un_)*dt) # must avoid recreating out!
      end
      function jnonlin!(un_,out)
        out[:,:]=eye(fcn.d)-(1-Alpha)*fcn.Jdrift(un_)*dt
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
end #sdelab2

# generalised OU
include("common.jl")
import SDELab2
################### problem parameters
T=1 # time step, final time
x0=1. # initial data
alp=-1/2 # drift parameter
bet= 4/5 # drift parameter
sig=sqrt(1/2) # diffusion param
quad0=SDELab2.delta_measure(x0)
###########################################
function drift(u,alp,bet)
  x=u[1]; mf=u[2]
  return alp*x+bet*mf
end
##########################################
function diff(u,sig)
  return sig
end
#######################################
function mf_r(x)
  return [x]
end
##########################################
fcn=SDELab2.set_fcn_mf(u->drift(u,alp,bet),
          u->diff(u,sig),
          x->mf_r(x))
# test function
phi=x->[x;x.^2]
###################################### exact
exact_mean=x0*exp((alp+bet));
exact_var=(sig^2)/(2*alp)*exp(2*alp)-sig^2/(2*alp);
exact_mom2=exact_var+exact_mean^2;
exact=[exact_mean,exact_mom2];
println("ode soln ",exact)
##################################### method parameters
dt=1e-3
solver="MFEM0"
params=SDELab2.set_opt(dt,solver)
params["cutoff"]=10 # radius R for Alg 3.1
params["const"]=1e3 # for no_Gauss_points
params["length_scale"]=1 # for no_Gauss_points
params["skiphalf"]=true
params["periodic"]=false
params["no_partition"]=1
dtpts=[0.2,0.1,0.05,0.025,0.01,0.005]#,0.0025,0.001]
errors=zeros(length(dtpts),2)
errors_no=deepcopy(errors)
errors_nd=deepcopy(errors)
cpu_times=zeros(length(dtpts))
cpu_times_no=deepcopy(cpu_times)
cpu_times_nd=deepcopy(cpu_times)


params["Solver"]="MFEM0"
dtpts_no=dtpts/4
for i=1:(length(dtpts))
    params["MaxStepSize"]=dtpts_no[i]
    qoi,cpu_time,stats,q=get_soln(fcn,quad0,params,phi,false)
    errors_no[i,:],cpu_times_no[i]=helper(qoi, exact, cpu_time, "normal")
end
########################
params["Solver"]="MFEM0"
for i=1:(length(dtpts))
    params["MaxStepSize"]=dtpts[i]
    qoi,cpu_time,stats,q=get_soln(fcn,quad0,params,phi,true)
    errors[i,:],cpu_times[i]=helper(qoi, exact, cpu_time, "normal")
end
################
params["Solver"]="MFEM2"
for i=1:(length(dtpts))
    params["MaxStepSize"]=dtpts[i]
    qoi,cpu_time,stats,q=get_soln(fcn,quad0,params,phi,false)
    errors_nd[i,:],cpu_times_nd[i]=helper(qoi, exact, cpu_time, "normal")
end
#########################################

figure(1);
# first functional
i=1
time_step_err(dtpts,errors[:,i],
  dtpts_no,errors_no[:,i],
  dtpts,errors_nd[:,i],[1e0,1e0],utf8("mom1.pdf"))
err_cpu(errors[:,i],cpu_times[:],
  errors_no[:,i],cpu_times_no[:],
  errors_nd[:,i],cpu_times_nd[:],[1e0,1e0], utf8("mom1_.pdf"))
# second functional
i=2
time_step_err(dtpts,errors[:,i],
      dtpts_no,errors_no[:,i],
      dtpts,errors_nd[:,i],[1e0,1],utf8("mom2.pdf"))
err_cpu(errors[:,i], cpu_times[:],
      errors_no[:,i],cpu_times_no[:],
      errors_nd[:,i],cpu_times_nd[:],[1e0,1], utf8("mom2_.pdf"))
#########################################
println("finished preparing plots")

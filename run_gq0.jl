# polynomial drift
include("common.jl")
import SDELab2
#
T=1 #  final time
x0=1. # initial data
alp=2. # parameter
sig=1 #
quad0=SDELab2.delta_measure(x0)
###########################################
function drift_ff(u,alp)
  x=u[1]; mf=u[2]; mf_sq=u[3]
  return alp*x+mf-x*mf_sq
end
##########################################
function diff_gg(u,sig)
  x=u[1]
  return sig*x
end
#########################################
function mf_r(x)
  return [x; x^2]
end
#########################################
fcn=SDELab2.set_fcn_mf(u->drift_ff(u,alp),
          u->diff_gg(u,sig),
          x->mf_r(x))
######################################
function odef(t, r)
    # Extract the coordinates from the r vector
    (x, y) = r
    dx_dt = (alp+1)*x- x*y
    dy_dt = (2*alp+1)*y + 2*x*x -2*y*y
    # Return the derivatives as a vector
    return [dx_dt; dy_dt]
end
####################################
using ODE
const odedt = 0.1
tvec = 0:odedt:T
tic()
tvec, pos = ode45(odef, [x0;x0], tvec,abstol=1e-9,reltol=1e-9)
exact=pos[end];
println("ode soln ",exact,"; cpu time ",toq()," secs.")
#############################################

dt=1e-2
solver="MFEM0"
params=SDELab2.set_opt(dt,solver)
params["cutoff"]=40 # range R for Alg 3.1 cutoff
params["const"]=1 # consts for no_Gauss_points
params["length_scale"]=10 # for no_Gauss_points
params["skiphalf"]=true
params["periodic"]=false
params["no_partition"]=1
dtpts=[0.05,0.025,0.01,0.005,0.0025]#,0.001]
errors=zeros(length(dtpts),2)
errors_no=deepcopy(errors);errors_nd=deepcopy(errors)
cpu_times=zeros(length(dtpts));
cpu_times_no=deepcopy(cpu_times);cpu_times_nd=deepcopy(cpu_times)
phi=x->[x;x.^2]
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
str1="--" # extrap
str2="k-." # 2nd-order
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

import SDELab2
#
T=1 # time step, final time
x0=1. # initial data
alp=-1 # parameter
bet=0
sig=0.5
quad0=SDELab2.delta_measure(x0)
###########################################
function drift_fff(u,alp,bet)
  x=u[1]; mf=u[2]
  return alp*x+bet*mf
end
##########################################
function diff_ggg(u,sig)
  x=u[1]
  return sig*x
end
#########################################
function mf_r(x)
  return [x]
end
#########################################
fcnq=SDELab2.set_fcn_mf(u->drift_fff(u,alp,bet),
          u->diff_ggg(u,sig),
          x->mf_r(x))

# test function
phi=x->[x;x.^2]
######################################
exact_mean=x0*exp((alp+bet));
exact_var=(sig^2)/(2*alp)*exp(2*alp)-sig^2/(2*alp);
exact_mom2=exact_var+exact_mean^2;
#exact_mom2=(exp(2*alp)*x0^2
#            +x0^2*exp(2*alp)*(exp(2*bet)-1)
#            +sig^2/(2*alp)*(exp(2*alp)-1))
exact=[exact_mean,exact_mom2];
println("ode soln ",exact)
##############################################
# MLMC implementation in MATLAB
mlmc_cpus= 1e3*[0.000142083385202,  0.000453974977685,   0.012572159661283, 0.051641483155560,   1.913720991516199,   8.940239873090615]
mlmc_errs=[ 0.008701679831122,   0.003964115036436,   0.000721699175333, 0.000323972021408,   0.000089514867964,   0.000056234777698]
#
dt=1e-3
solver="MFEM0"
params=SDELab2.set_opt(dt,solver)
params["eta"]=1 # determines timesteps
params["cutoff"]=10 # for diameter reduction
params["const"]=1e3 # for no_Gauss_points
params["length_scale"]=1 # for no_Gauss_points
params["skiphalf"]=true
params["periodic"]=false
params["no_partition"]=1
dtpts=[0.2,0.1,0.05,0.025,0.01,0.005,0.0025,0.001,0.0005]
errors=zeros(length(dtpts),2)
errors_no=deepcopy(errors)
errors_nd=deepcopy(errors)
cpu_times=zeros(length(dtpts))
cpu_times_no=deepcopy(cpu_times)
cpu_times_nd=deepcopy(cpu_times)




params["Solver"]="MFEM0"
for i=1:(length(dtpts))
    params["MaxStepSize"]=dtpts[i]/4
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
  dtpts,errors_no[:,i],
  dtpts,errors_nd[:,i],[1e0,1e0],utf8("mom1.pdf"))
err_cpu(errors[:,i],cpu_times[:],
  errors_no[:,i],cpu_times_no[:],
  errors_nd[:,i],cpu_times_nd[:],[1e0,1e0], utf8("mom1_.pdf"))
# second functional
i=2
time_step_err(dtpts,errors[:,i],
      dtpts,errors_no[:,i],
      dtpts,errors_nd[:,i],[1e0,1],utf8("mom2.pdf"))
err_cpu(errors[:,i], cpu_times[:],
      errors_no[:,i],cpu_times_no[:],
      errors_nd[:,i],cpu_times_nd[:],[1e0,1], utf8("mom2_.pdf"))
#########################################
println("finished preparing plots")
#

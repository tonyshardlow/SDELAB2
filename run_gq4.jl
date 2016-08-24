# plane rotator
import SDELab2
include("common.jl")
########################## problem parameters
T=1
K=1;
tau=1/8; sig=sqrt(2*tau)
init_mu_sig_sq=[pi/2 3*pi/4]
quad0=gauss_quadrature.Gauss_Hermite_per(init_mu_sig_sq,40,2*pi)
#########################################
function drift_ff(u,K)
  x=u[1]; mfs=u[2]; mfc=u[3]
  return K*mfs*cos(x)-(K*mfc+1)*sin(x)
end
#########################################
function diff_gg(u,sig)
  return sig
end
#########################################
function mf_r(x)
  return [sin(x); cos(x)]
end
#########################################
fcn=SDELab2.set_fcn_mf(u->drift_ff(u,K),
          u->diff_gg(u,sig),
          x->mf_r(x))
#########################################
######################################### method parameters
solver="MFEM2"
dt=1e-2
params=SDELab2.set_opt(dt,solver)
params["const"]=1           # for no_Gauss_points
params["length_scale"]=2*pi # lam=1/length^2
params["skiphalf"]=true
params["periodic"]=true
params["no_partition"]=5
#########################################
###################### computing reference value
print("Computing a reference value...")
phi=x->[sin(x)^2;sin(x)]
tic()
params["extrap"]=false
dtpts=0.1./(2.0).^(1:6)
params["MaxStepSize"]=dtpts[end]
exact,cpu_time,stats,q1,q0=get_soln(fcn,quad0,params,phi,false)
println(stats["no_pts"])
println("exact=", exact)
#########################################
initial_final_cdf_pdf(q0,q1,[100,40],[0.08,0.01])

#########################################
dtpts=dtpts[1:end-1] # omit last value, used for reference
errors=zeros(length(dtpts),2); errors_no=deepcopy(errors); errors_nd=deepcopy(errors)
cpu_times=zeros(length(dtpts)); cpu_times_no=deepcopy(cpu_times);cpu_times_nd=deepcopy(cpu_times)
params["MaxStepSize"]=dtpts[1]
######################
params["Solver"]="MFEMe"
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
  dtpts,errors_nd[:,i],[1,1],utf8("mom1.pdf"))
err_cpu(errors[:,i],cpu_times[:],
  errors_no[:,i],cpu_times_no[:],
  errors_nd[:,i],cpu_times_nd[:],[1,1], utf8("mom1_.pdf"))
# second functional
i=2
time_step_err(dtpts,errors[:,i],
      dtpts_no,errors_no[:,i],
      dtpts,errors_nd[:,i],[1,1],utf8("mom2.pdf"))
err_cpu(errors[:,i], cpu_times[:],
      errors_no[:,i],cpu_times_no[:],
      errors_nd[:,i],cpu_times_nd[:],[1,1], utf8("mom2_.pdf"))
#########################################
println("finished preparing plots")

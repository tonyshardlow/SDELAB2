# Burgers
import SDELab2
include("common.jl")
########################## problem parameters
T=1
diffusion=0.1
sig=sqrt(2*diffusion)

quad0=gauss_quadrature.delta_measure(0.)
#quad0=gauss_quadrature.Gauss_Hermite([0.,1e-8],4)
#########################################
function drift_ff(u)
  x=u[1]; mf=u[2]
  if x>mf
    return 0.
  elseif x==mf
    return 0.5
  else
    return 1.
  end
end
#
############################
function drift_fr(u,ell)
  x=u[1]; mf=u[2]
  return 0.5*erfc((x-mf)/ell)
end
#########################################
function diff_gg(u,sig)
  return sig
end
#########################################
function mf_r(x)
  return [x]
end
#########################################
fcn=SDELab2.set_fcn_mf(u->drift_fr(u,ell),
          u->diff_gg(u,sig),
          x->mf_r(x))
#########################################
function exact_cdf2(x,t,sig)
  sig_sq=sig^2
  scale=sig*sqrt(2*t)
 return 1- ( erfc(-x/scale)./
 (
 erfc(-x/scale)
 + exp((t-2*x)/(2*sig_sq)) .* (2-erfc((t-x)/scale))
 )    )
end
#
function exact_cdf3(x,t,sig)
  sig_sq=sig^2
  scale=sig*sqrt(2*t)
  T1=erfc(-x/scale)
  T2=exp((t-2*x)/(2*sig_sq))
  T3=erfc((t-x)/scale)
 return T2*(2-T3)./(T1+T2.*(2-T3))
end
######################################### method parameters
solver="MFEM0e"
dt=1e-2
params=SDELab2.set_opt(dt,solver)
params["const"]=1e0       # for no_Gauss_points
params["length_scale"]=1 # lam=1/length^2
params["cutoff"]=1  #
params["skiphalf"]=true
params["periodic"]=false
params[ "no_partition"]=20
#########################################
###################### computing reference value
print("Computing a reference value...")
phi=x->[x^2;x]
tic()
params["extrap"]=false
dtpts=0.1./(2.0).^(1:8 )
params["MaxStepSize"]=dtpts[end]
ell=sqrt(dtpts[end])

#ell=0.5#(dtpts[end])
ell=0.1
#ell=0.01
#ell=0.001
quad0=gauss_quadrature.delta_measure(0.)#-dtpts[end]/2)
fcn=SDELab2.set_fcn_mf(u->drift_fr(u,ell),
          u->diff_gg(u,sig),
          x->mf_r(x))
#quad0=gauss_quadrature.delta_measure(sqrt(dtpts[end]))
exact,cpu_time,stats,q1=get_soln(fcn,quad0,params,phi,true)
println(stats["no_pts"])
println("exact=", exact)
#########################################
#initial_final_cdf_pdf(q0,q1,[100,40],[0.08,0.01])
#

cum1=cumsum(q1.w)
x1=(q1.X[2:end]+q1.X[1:end-1])/2
#x1=q1.X
#x1=q1.X[1:end-1]
cum_exact=1-exact_cdf2(q1.X,1.,sig)
q2=gauss_quadrature.QUAD(length(x1),x1,q1.w)
println("l1 ",gauss_quadrature.l1cdf(q2,z->(1-exact_cdf2(z,1.,sig))))
using PyPlot
figure(2)
clf()
plot(q1.X,cum_exact,"b")
plot(x1,cum1[2:end],"g")#
axis([-1.5,1.5,0,1])
xlabel("X")
ylabel("cdf")
gcf()[:set_size_inches](width_inches,height_inches)
tight_layout()
  savefig(utf8("cdf_cmp.pdf"))
  clf()
  plot(q1.X,abs(cum_exact-cum1),"b")
  #plot(q1.X,cum1[2:end],"g")#
    axis([-1.5,1.5,0,0.04])
  xlabel("X")
  ylabel("error")
  gcf()[:set_size_inches](width_inches,height_inches)
  tight_layout()
  savefig(utf8("cdf_cmp_del.pdf"))
#########################################
dtpts=dtpts[1:end-1] # omit last value, used for reference
errors=zeros(length(dtpts),2); errors_no=deepcopy(errors); errors_nd=deepcopy(errors)
cpu_times=zeros(length(dtpts)); cpu_times_no=deepcopy(cpu_times);cpu_times_nd=deepcopy(cpu_times)
params["MaxStepSize"]=dtpts[1]
######################
params["Solver"]="MFEM0e"
dtpts_no=dtpts
for i=1:(length(dtpts))
    params["MaxStepSize"]=dtpts_no[i]
    qoi,cpu_time,stats,q=get_soln(fcn,quad0,params,phi,false)
    errors_no[i,:],cpu_times_no[i]=helper(qoi, exact, cpu_time, "normal")
end
########################
params["Solver"]="MFEM0e"
for i=1:(length(dtpts))
    params["MaxStepSize"]=dtpts[i]
    qoi,cpu_time,stats,q=get_soln(fcn,quad0,params,phi,true)
    errors[i,:],cpu_times[i]=helper(qoi, exact, cpu_time, "normal")
end
################
#########################################
figure(1);
# first functional
i=1

time_step_err(dtpts,errors[:,i],
  dtpts_no,errors_no[:,i],
  dtpts,errors_nd[:,i],[2,2],utf8("mom1.pdf"))
err_cpu(errors[:,i],cpu_times[:],
  errors_no[:,i],cpu_times_no[:],
  errors_nd[:,i],cpu_times_nd[:],[1,1], utf8("mom1_.pdf"))
# second functional
i=2
time_step_err(dtpts,errors[:,i],
      dtpts,errors_no[:,i],
      dtpts,errors_nd[:,i],[1,1],utf8("mom2.pdf"))
err_cpu(errors[:,i], cpu_times[:],
      errors_no[:,i],cpu_times_no[:],
      errors_nd[:,i],cpu_times_nd[:],[1,1], utf8("mom2_.pdf"))
#########################################
println("finished preparing plots")

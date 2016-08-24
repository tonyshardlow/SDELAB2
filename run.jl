using SDELab2 # SDELab2.jl should be findable via LOAD_PATH
println("run.jl")

function drift_f(u,r)
  return (r*u)
end
#
function diff_g_mat(u,dw,alpha)
  return alpha*(u.*dw)
end
#
function diff_g_vecs(alpha)
  m=2
  A=Array(Function,m)
  g1(u)=[u[1]*alpha,0]
  g2(u)=[0,u[2]*alpha]
  A[1]=g1
  A[2]=g2
  return A
end
#
function gbm_exact(y0,r,sig,t,w,ito_flag)
  if ito_flag==true
    out=exp((r-0.5*sig^2)*t+sig*w).*y0
  else
    out=exp(r*t+sig*w).*y0
  end
end
#
function get_error(dt,solver,test_case)
  r=test_case["r"]
  sig=test_case["sig"]
  noise_type=test_case["noise_type"]
  d=2
  m=2
  tf=1
  u=test_case["initial_data"]
  fcn=set_fcn(d,m,x->drift_f(x,r),(x,dw)->diff_g_mat(x,dw,sig),diff_g_vecs(sig), noise_type)
  opt=set_opt(dt,solver)
  t,y,W=sde_strong_solution(fcn, [0,tf], u, opt)
  exact=gbm_exact(u,r,sig,t[end],W,opt["is_ito"])
  yf=reshape(y[end,:],fcn.d)
  return norm(exact-yf)/norm(exact)# relative error
end
test_case1=Dict("r"=>-0.1, "sig"=>0.5, "noise_type"=>2, "initial_data"=>[1.,1.])
#

#fcn=set_fcn(2,2,x->drift_f(x,r),(x,dw)->diff_g_mat(x,dw,sig),diff_g_vecs(1), 2)
#opt=set_opt(0.01,"SIBDF")
#opt["nlsolve_param"].xtol=1e-5

#t,y,W=sde_strong_solution(fcn, linspace(0.0,1.0,10000), u, opt)
#exact=gbm_exact(u,r,sig,t[end],W,opt["is_ito"])
#
#plot(t,y)
#print(W,"u=",u)
dt_vals=[0.01,0.005,0.002,0.001,0.0005,1e-4]
error=similar(dt_vals)
all_error=Dict()
M=10
for key in keys(SDELab2.method_code)
  print("Integrator:", method_code[key])
  for i in 1:length(dt_vals)
    sum=0
    for m in 1:M
      err=get_error(dt_vals[i],key,test_case1)
      sum+=err^2
   end
  error[i]=sqrt(sum/M)
  end
  all_error[key]=deepcopy(error)
end
using PyPlot
for key in keys(method_code)
  loglog(dt_vals,all_error[key],label=method_code[key])
  println(method_code[key],all_error[key])
end
legend()

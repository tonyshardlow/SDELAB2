using SDELab2 # SDELab2.jl should be findable via LOAD_PATH
println("run_vdp2.jl")
type vdp_param
  alpha
  beta
  A
  B
end

function drift_f(u,p)
  return [u[2], p.alpha*u[1]+p.beta*u[2]-p.A*u[1]^3-p.B*u[1]^2*u[2]]
end
#
function diff_g_mat(u,dw,alpha)
  return [0;alpha*u[1]*dw]
end
#
function diff_g_vecs(alpha)
  m=1
  A=Array(Function,m)
  g1(u)=[0,u[1]*alpha]
  A[1]=g1
  return A
end
# dimensions
d=2
m=1
noise_type=3 # commutative
# params
p=vdp_param(-1,0.1,1.,1.)
sig=0.1
# initial data
u=[0.,0.001]
t0=0
tf=500
#
fcn=set_fcn(d,m,x->drift_f(x,p),(x,dw)->diff_g_mat(x,dw,sig),diff_g_vecs(sig), 3)
opt=set_opt(0.01,"SIBDF")
opt["do_plotting"]=true
opt["output_plot_type"]="phase_plot"
opt["output_sel"]=[1,2]
#using PyPlot
println("Start sde_strong_soln")
tic(); t,y,W=sde_strong_solution(fcn, linspace(t0,tf,10000), u, opt); run_time=toc()

#

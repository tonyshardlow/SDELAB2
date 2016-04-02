using SDELab2 # SDELab2.jl should be findable via LOAD_PATH
println("run_HH.jl")
#
# Hodgkin-Huxley example from
# H. Alzubaidi and T. Shardlow. Numerical simulations of SDEs and SPDEs.
# In: Stochastic Methods in Neuroscience. Ed. by C. Laing and G. Lord. OUP,
# 2009. Chap. 12, pp. 344–366.
# doi: 10.1093/acprof:oso/9780199235070.003.0012 .
#
function drift_f(y,Mu)
  g_k=36; g_Na=120; g_L=0.3
  C_m=1 # membrane capabilities
  C1=g_k/C_m; C2=g_Na/C_m; C3=g_L/C_m
  # constants of resting potential of such iterations
  v_k=-12; v_Na=115; v_L=10.613
  a_n=(10-y[1]) / (100*(exp((10-y[1])/10)-1))
  b_n=exp(-y[1]/80)/8
  a_m=(25-y[1])/(10*(exp((25-y[1])/10)-1.))
  b_m=4*exp(-y[1]/18)
  a_h=(7*exp(-y[1]/20))/100
  b_h=1/(exp((30-y[1])/10)+1)
  #
  # compute drift function for HH model
  z1=(C1*(v_k-y[1])*(y[2]^4)) + (C2*(v_Na-y[1])*(y[3]^3)*y[4]) + (C3*(v_L-y[1])) + Mu
  z2=a_n*(1-y[2])-b_n*y[2]
  z3=a_m*(1-y[3])-b_m*y[3]
  z4=a_h*(1-y[4])-b_h*y[4]
  return [z1,z2,z3,z4]
end

#
function diff_g_mat(u,dw,sigma)
  return   [sigma*dw;0;0;0]
end
#
function diff_g_vecs(sig)
  m=1
  A=Array(Function,m)
  g1(u)=[sig,0,0,0]
  A[1]=g1
  return A
end
# dimensions
d=4
m=1
noise_type=3 # commutative
# params
sig=0.1
Mu=2.5
# initial data
a0_n=0.1/(exp(1.)-1.); b0_n=0.125
a0_m=2.5/(exp(2.5)-1); b0_m=4;
a0_h=0.07; b0_h=1/(exp(3.)+1)
y0=[0; a0_n/(a0_n+b0_n); a0_m/(a0_m+b0_m); a0_h/(a0_h+b0_h)]
#
t0=0
tf=40

#
fcn=set_fcn(d,m,x->drift_f(x,Mu),(x,dw)->diff_g_mat(x,dw,sig),diff_g_vecs(sig), noise_type)
opt=set_opt(0.01,"SIBDF")
opt["do_plotting"]=true
opt["output_plot_type"]="path_plot"
opt["output_sel"]=[1,2]
#using PyPlot
println("Start sde_strong_soln")
tic();t,y,W=sde_strong_solution(fcn, linspace(t0,tf,100), y0, opt);run_time=toc()


#

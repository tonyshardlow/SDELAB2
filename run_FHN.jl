#push!(LOAD_PATH,utf8("T:\\Research\\sdelap") )
using SDELab2 # SDELab2.jl should be findable via LOAD_PATH
println("run_FHN.jl")

#
# FitzHugh-Nagumo example from
#
# H. Alzubaidi and T. Shardlow. Numerical simulations of SDEs and SPDEs.
# In: Stochastic Methods in Neuroscience. Ed. by C. Laing and G. Lord. OUP,
# 2009. Chap. 12, pp. 344–366.
# doi: 10.1093/acprof:oso/9780199235070.003.0012 .
#
type fhn_params
  D
  Mu
  a
  b
  gamma
end


function drift_f(y,p::fhn_params)
  d=round(Int,length(y)/2)
  A=SymTridiagonal(-2*ones(d),ones(d-1))
  u=y[1:d] # fasts
  v=y[(d+1):end]
  F=(u.*(1-u).*(u-p.a))-v
  B=p.D*((d+1)^2)*A*u
  z1=B+F+p.Mu
  z2=p.b*(u-p.gamma*v)
  return[z1;z2]
end
#
#
function diff_g_mat(u,dw,sigma)
  d=length(dw)
  Q=sqrt(d+1) # orthognoal matrix does not change distribution
  return   [sigma*Q*dw;zeros(d)]
end
#
function diff_g_vecs(sig,m)
  A=Array(Function,m)
  for i=1:m
    A[i]=(u->1)
  end
  # do not use Milstein methods (no point for additive!)
  return A
end
# dimensions
d=9
d_all=2*d
noise_type=2 # diagonal
# params
p=fhn_params(0.01,0.5,0.05,0.008,0.5)
sig=0.05
# initial data
a0_n = -.1/(exp(1.)-1.); b0_n=0.125
v0=a0_n/(a0_n+b0_n)
y0=[zeros(d);v0*ones(d)]
#
t0=0
tf=200
#
fcn=set_fcn(d_all,d, x->drift_f(x,p),(x,dw)->diff_g_mat(x,dw,sig),diff_g_vecs(sig,d), noise_type)
opt=set_opt(0.01,"SIBDF")
opt["do_plotting"]=true
opt["output_plot_type"]="path_plot"
opt["output_sel"]=1:d
using PyPlot
println("Start sde_strong_soln")
tic();t,y,W=sde_strong_solution(fcn, linspace(t0,tf,100), y0, opt);run_time=toc()

#

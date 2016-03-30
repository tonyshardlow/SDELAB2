# SDELAB2

This is a reimplementation of SDELAB, a MATLAB package for generating sample paths of initial-value problems for SDEs numerically,
as described in 

H. Gilsing and T. Shardlow. SDELab: a package for solving stochastic
differential equations in MATLAB. In: J. Comput. Appl. Math. 205.2 (2007),
pp. 1002–1018. issn: 0377-0427.
doi: 10.1016/j.cam.2006.05.037 .

SDELAB2 is implemented in Julia and provides the functionality described in this paper, though there are changes in syntax.
The key differences are
- drift and diffusion are autonomous (non-automous functions needs to disguised as an autonomous SDE)
- you specify the drift function (u->f(u)), 
diffusion-matrix function ((u,dw)->g(u) dw),
and diffusion-vector function (dif_vec[i]=(u->gi(u)) (i.e., a vector of functions u->gi(u0).

- automatic differentiation creates Jacobian fns for J (for nonlinear solve) 
and for the diffusion-vector function (for Milstein methods)

- the main driver routine is
t,y,W=sde_strong_solution(fcn, tspan, u0, opt)

where fcn is a special data type for describing the fcn, created by 

fcn=set_fcn(d,m,x->drift_f(x,p),(x,dw)->diff_g_mat(x,dw,sig),diff_g_vecs(sig), noise_type)

for d=phase-space dimension, m=number of Brownian motions, 
drift_f, diff_g_mat, diff_g_vecs are the coefficient functions for the SDE, 
noise_type=1 (unstructured), =2(diagonal), =3(commutative) (see the above paper for the meaning of this).

tspan=well-ordered vector of times where the solution is required

u0=d-vector of initial data

opt=set_opt(dt,method), where dt=maximum time step, 
method="SIE0" or "SIE1" or "SIE" (Euler-Maruyama for explicit, implicit, or alpha-implicit),
"SSEH0" or "SSEH1" or "SSEH" (explicit Stratonovich Heun for explicit, implicit, or alpha-implicit), 
"SIM0" or "SIM1" or "SIM" (for Ito Milstein),
"SSM0" or "SSM1" or "SSM" (Stratonovich Milstein), or 
"SIBDF" (second-order BDF). See above-referenced paper. All use auto-computed derivatives.

opt can be further adjusted to set
opt["Alpha"]=parametr for alpha-implicit methods (default alpha=0.5; trapezium rule)
opt["nlsolve_param"].maxfeval and .ftol and .xtol    to control nonlinear sovler (default 10,1e-6,1e-8)
opt["calc_final_w"] gives final position W of driving Wiener process (useful for calculating exact soln of geometric BM)
opt["do_plotting"]=true or false, for automatic plotting via PyPlot
opt["always_clear"]=true or false, to clear current plot each time
opt["output_sel"]= vector of indicies to plot (default to [1,2])
opt["output_plot_type"]="path_plot" or "phase_plot" ("time_phase" and "phase3_plot" not yet implemented).

THis depends on the Julia package NLsolve and PyPlot.

The package is currently being refined. The two examples presented (geometric brownian motion and van der pol, illustrate its use).


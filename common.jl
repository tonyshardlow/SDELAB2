using PyPlot
using SDELab2
# preferred pdf size
width_inches=4; height_inches=3
#######################################
function get_soln(fcn,quad0,params,phi,extrap_flag)
  println("Computing with dt=",params["MaxStepSize"])
  println(method_code[params["Solver"]])
  tic()
  tspan=[0,T]
  if extrap_flag==true
    quad,stats=SDELab2.sde_gq_solution_extrap(fcn,tspan,
                                      quad0,params)
  else
    quad,stats=SDELab2.sde_gq_solution(fcn,tspan,
                                      quad0,params)

  end
  cpu_time=toq()
  println("Mean number of points ",mean(stats["no_pts"]))
  println(stats["no_pts"])
  println("cpu time ", cpu_time)
  qoi=gauss_quadrature.do_quad(phi,quad)
  return qoi,cpu_time,stats,quad
end
#########################################
function helper(qoi,exact,cpu_time,str)
  println(str)
  println("qoi",qoi)
  error=abs(exact-qoi)
  println("error ", error)
  return error, cpu_time
end
##########################################
function initial_final_cdf_pdf(q0,q1, N, s)
  str2="k-"
  str2b="b:"
  str1="b-"
  clf()
  x,pdf,cdf,x1,cdf1=gauss_quadrature.quad_to_pdf(q1,N[1], s[1])
  q00=gauss_quadrature.QUAD((q0.length+1), [0;q0.X],[0;q0.w])
  x0,pdf0,cdf0,x10,cdf10=gauss_quadrature.quad_to_pdf(q00,N[2],s[2])
  plot(x,pdf,str1)
  plot(x0,pdf0,str2)
  xlabel(utf8("X"))
  ylabel(utf8("pdf"))
  axis([minimum(x),maximum(x),0,maximum([pdf;pdf0]) ])
  gcf()[:set_size_inches](width_inches,height_inches)
  tight_layout()
  savefig(utf8("pdf.pdf"))
  #
  clf()
  plot(x,cdf,str1)
  plot(x1,cdf1,str2b)
  plot(x10,cdf10,str2)
  axis([minimum(x),maximum(x),0,maximum([cdf;cdf0]) ])
  xlabel(utf8("X"))
  ylabel(utf8("cdf"))
  gcf()[:set_size_inches](width_inches,height_inches)
  tight_layout()
  savefig(utf8("cdf.pdf"))
end
#

function initial_final_cdf_pdf(q0,q1)
  str2="k-"
  str2b="b:"
  str1="b-"
  clf()
  x,pdf,cdf,x1,cdf1=gauss_quadrature.quad_to_pdf(q1)#100,0.08)
  q00=gauss_quadrature.QUAD((q0.length+1), [0;q0.X],[0;q0.w])
  x0,pdf0,cdf0,x10,cdf10=gauss_quadrature.quad_to_pdf(q00)#40,0.01)
  plot(x,pdf,str1)
  plot(x0,pdf0,str2)
  xlabel(utf8("X"))
  ylabel(utf8("pdf"))
  axis([minimum(x),maximum(x),0,maximum([pdf;pdf0]) ])
  gcf()[:set_size_inches](width_inches,height_inches)
  tight_layout()
  savefig(utf8("pdf.pdf"))
  #
  clf()
  plot(x,cdf,str1)
  plot(x1,cdf1,str2b)
  plot(x10,cdf10,str2)
  axis([minimum(x),maximum(x),0,maximum([cdf;cdf0]) ])
  xlabel(utf8("X"))
  ylabel(utf8("cdf"))
  gcf()[:set_size_inches](width_inches,height_inches)
  tight_layout()
  savefig(utf8("cdf.pdf"))
end
#
function loglog_safe(a,b,s)
  if length(a)==length(b)
    loglog(a,b,s)
  end
end
#
function time_step_err(dt1,e1,dt2,e2,dt3,e3,s,filename)
  clf()
  str1="b--" # extrap
  str2="g-"
  str3="k-." # 2nd-order
  axis("auto")
  plt[:rcParams]["xtick.labelsize"]="small"
  plt[:rcParams]["ytick.labelsize"]="small"

  loglog_safe(dt1,e1,str1)#dtpts,errors[:,i],str1)
  loglog_safe(dt2,e2,str2)#dtpts,errors_no[:,i])
  loglog_safe(dt3,e3,str3)#dtpts,errors_nd[:,i],str2)
  loglog_safe(dt1,s[1]*e1[1]*(dt1/dt1[1]).^(2),"y-")
  loglog_safe(dt2,s[2]*e2[1]*(dt2/dt2[1]).^(1),"y-")
  xlabel(utf8("time step"))
  ylabel(utf8("error"))
  axis([minimum([dt1;dt2;dt3]), maximum([dt1;dt2;dt3]), minimum([e1;e2;e3]), maximum([e1;e2;e3]) ] )
  axis([5e-4,1e-1,1e-8,1e-2])
  gcf()[:set_size_inches](width_inches,height_inches)
  tight_layout()
  savefig(filename)
end
#
function err_cpu(e1,c1,e2,c2,e3,c3,s,filename)
  clf()
  axis("auto")
  str1="b--" # extrap
  str2="g-"
  str3="k-." # 2nd-order
  plt[:rcParams]["xtick.labelsize"]="small"
  plt[:rcParams]["ytick.labelsize"]="small"
  loglog_safe(e1,c1,str1)#dtpts,errors[:,i],str1)
  loglog_safe(e2,c2,str2)#dtpts,errors_no[:,i])
  loglog_safe(e3,c3,str3)#dtpts,errors_nd[:,i],str2)
  loglog_safe(e1,s[1]*c1[1]*(e1/e1[1]).^(-1/2),"y-")
  loglog_safe(e2,s[2]*c2[1]*(e2/e2[1]).^(-1),"y-")
  xlabel(utf8("error"))
  ylabel(utf8("cpu time (secs)"))

  axis([minimum([e1;e2;e3]), maximum([e1;e2;e3]), minimum([c1;c2;c3]), maximum([c1;c2;c3]) ] )
  gcf()[:set_size_inches](width_inches,height_inches)
  tight_layout()
  savefig(filename)
end

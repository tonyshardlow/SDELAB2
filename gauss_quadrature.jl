__precompile__()
module gauss_quadrature
using Dierckx
export  do_quad, delta_measure,Gauss_Hermite, Gauss_Hermite_per,QUAD, Golub_Welsch!,periodic_simplify!,simplify!,sortquad!,quad_to_pdf
##########################
function init(fcn, opt)
  #=
  =#
end
#######################
type QUAD
  length
  X
  w
end
#######################
function do_quad(test,quad::QUAD)
  # evaluate test fn wrt to quad
  #
  sum=0.0
  for i=1:quad.length
    sum += test(quad.X[i]) * quad.w[i]
  end
  return sum
end
##############################################
function lanczos(N,xw)
  # Adapted from Gautschi's OPQ Matlab routines
  Ncap=size(xw,1)
  if N<=0 || N>Ncap
      print("N:",N,"; Ncap:",Ncap)
     error("N out of range")
  end
  p0=xw[:,1]; p1=zeros(Ncap,1); p1[1]=xw[1,2]
  for n=1:Ncap-1
    pn=xw[n+1,2]; gam=1; sig=0; t=0; xlam=xw[n+1,1]
    for k=1:n+1
      rho=p1[k]+pn; tmp=gam*rho; tsig=sig

      if (rho)<=0
        gam=1; sig=0
      else
        gam=p1[k]/rho
        sig=pn/rho
      end
      tk=sig*(p0[k]-xlam)-gam*t
      p0[k]=p0[k]-(tk-t); t=tk
      if (sig)<=0
        pn=tsig*p1[k]
      else
        pn=(t^2)/sig
      end
      tsig=sig; p1[k]=tmp
    end
  end
  ab=[p0[1:N] p1[1:N]]
  return ab
end
######################################
function gauss(N,ab)
  # Adapted from Gautschi's OPQ Matlab routines
  N0=size(ab,1)
  if N0<N
    error("input array ab too short")
  end
  J=zeros(N,N)
  for n=1:N
     J[n,n]=ab[n,1]
  end
  for n=2:N
    J[n,n-1]=sqrt(ab[n,2])
    J[n-1,n]=J[n,n-1]
  end
  D,V=eig(J)
  #V=V(:,I);
  xw=[D ab[1,2]*V[1,:]'.^2]
  return xw
end
########################################
function Golub_Welsch!(quad::QUAD, m)
  # Compute m-point GQ
  # No change if input quad has too few points
  if quad.length> m
    m=round(Int16,m)
    ab=lanczos(m, [quad.X quad.w])
    xw=gauss(m,ab)
    quad.length=m
    quad.X=xw[:,1]
    quad.w=xw[:,2]
  end
end
###########################
function sortquad!(quad::QUAD)
  # sort quad points
  perm=sortperm(quad.X)
  quad.X=quad.X[perm]
  quad.w=quad.w[perm]
end
###########################
function simplify!(quad::QUAD, m, cutoff)
# Return reduced quadrature
# Group points outside interval [cutoff[1],cutoff[2]]
# Apply m-point GQ on interior
# No action if input quad less than m
  if quad.length > m
     ind=find(x->(x<=cutoff[1]), quad.X) # small nodes
     m1=length(ind); w1=sum(quad.w[ind])
     ind=find(x->(x>=cutoff[2]), quad.X) # big nodes
     m2=length(ind); w2=sum(quad.w[ind])
     if (quad.length-m1-m2>m)
       ind=find(y->(y>cutoff[1] && y<cutoff[2]), quad.X) # interior nodes
       ab=lanczos(m, [quad.X[ind] quad.w[ind]])
       xw=gauss(m,ab)
       quad.length=m+2
       quad.X=[cutoff[1]; xw[:,1]; cutoff[2]]
       quad.w=[w1;        xw[:,2]; w2] # update quad
     end
  end
end
#################################################
function partition_simplify!(quad::QUAD, m, d, cutoff)
# Apply periodic conditions and find m-point GQ
# on d sub-intervals.
  if quad.length > m
   ind=find(x->(x<cutoff[1]), quad.X) # small nodes
   m1=length(ind); w1=sum(quad.w[ind])
   ind=find(x->(x>cutoff[2]), quad.X) # big nodes
   m2=length(ind); w2=sum(quad.w[ind])
   if (quad.length-m1-m2>m)
    partition=collect(linspace(cutoff[1], cutoff[2], d+1))
    x=[]; w=[];
    for i=1:d
      inds=find(y->(y>partition[i] && y<=partition[i+1]), quad.X)
      quad1=QUAD(length(inds),quad.X[inds],quad.w[inds])
      Golub_Welsch!(quad1, max(5,m*min(1,10*sum(quad1.w) ) ) )
      x=[x; quad1.X]
      w=[w; quad1.w]
    end
    x=[cutoff[1]; x; cutoff[2]]
    w=[w1; w; w2]
    quad.length=length(x)
    quad.X=x
    quad.w=w
  end
 end
end
#####
function l1(q1::QUAD, q2::QUAD)
  mysum=0.
  i1=2; i2=2; c1=cumsum(q1.w); c2=cumsum(q2.w)
  if (q1.X[i1] < q2.X[i2]) # find smallest point
    b=q1.X[i1];  i1=i1+1
  else
    b=q2.X[i2];  i2=i2+1
  end
  while (i1<=q1.length && i2<=q2.length)
    if (q1.X[i1]<q2.X[i2]) # q1 is next smallest
      delta=q1.X[i1]-b; # work out delta X
      # interpolate on q2
      if (q2.X[i2]>q2.X[i2-1])
        lambda=(q1.X[i1]-q2.X[i2-1])/(q2.X[i2]-q2.X[i2-1])
        cinterp=lambda*c2[i2]+(1-lambda)*c2[i2-1]#
        c_delta=abs(cinterp-c1[i1]);
      else
        c_delta=0
      end
      # record current x
      b=q1.X[i1]; i1+=1
    else
      delta=q2.X[i2]-b; # work out delta X
      # interpolate on q2
      if q1.X[i1-1]<q1.X[i1]
        lambda=(q2.X[i2]-q1.X[i1-1])/(q1.X[i1]-q1.X[i1-1])
        cinterp=lambda*c1[i1]+(1-lambda)*c1[i1-1]
        c_delta=abs(cinterp-c2[i2])
      else
        c_delta=0
      end
      # record current x
      b=q2.X[i2]; i2+=1
    end
    mysum+=c_delta*delta
  end
  if (i1<q1.length)
    mysum+=sum(abs(c1[i1:end-1]).*(q1.w[i1+1:end]-q1.w[i1:end-1]))
  end
  if (i2<q2.length)
    mysum+=sum(abs(c2[i2:end-1]).*(q2.w[i2+1:end]-q2.w[i2:end-1]))
  end
  return mysum
end
#

#####
function l1cdf(q1::QUAD, cdf)
  mysum=0.
  c1=cumsum(q1.w)
  for i=2:q1.length,
    delta1=abs(c1[i-1]-cdf(q1.X[i-1] ))
    delta2=abs(c1[i]  -cdf(q1.X[i]) )
    mysum+= (q1.X[i]-q1.X[i-1])*(delta1+delta2)/2
  end
  return mysum
end
#################################################
function periodic_simplify!(quad::QUAD, m, length_scale)
# Apply periodic conditions and find m-point GQ
  quad.X=mod(quad.X, length_scale)
  Golub_Welsch!(quad, m)
end
#################################################
function partition_periodic_simplify!(quad::QUAD, m, d, length_scale)
# Apply periodic conditions and find m-point GQ
# on d sub-intervals.
  theta=mod(quad.X, length_scale)
  partition=collect(linspace(-1e-6, length_scale, d+1))
  x=[]; w=[];
  for i=1:d
    inds=find(y->(y>partition[i] && y<=partition[i+1]), theta)
    quad1=QUAD(length(inds),theta[inds],quad.w[inds])
    Golub_Welsch!(quad1, m)
    x=[x; quad1.X]
    w=[w; quad1.w]
  end
  quad.length=length(x)
  quad.X=x
  quad.w=w
end
###########################
function Gauss_Hermite(mu_sig_sq, n)
  # Returns m-point Gauss-Hermite quadrature rule
  # for N(m,s2), for mu_sig_sq=[m,s2]
  #
  # using Hermite recursion rule for w(x)=e^(-x^2)
  # and change of variable
  #
  println("GH: [mean,variance]=",mu_sig_sq)
  v =sqrt((1:n-1)/2.0);
  J=diagm(v,1)+diagm(v,-1)
  D,V=eig(J)
  xw=[D sqrt(pi)*V[1,:]'.^2]
  # change of variable
  f=x->mu_sig_sq[1]+x*sqrt(2*mu_sig_sq[2])
  weight_scale=sqrt(pi)
  quad=QUAD(n,f(xw[:,1]),xw[:,2]/weight_scale)
  #
  #
  println(quad.length,"-pt Gauss-Hermite rule.")
  sortquad!(quad)
  return quad
end
###########################
function Gauss_Hermite_per(mu_sig_sq, m, length_scale)
  # Returns m-point Gauss-Hermite quadrature rule
  # for N(m,s2), for mu_sig_sq=[m,s2]
  #
  # using Hermite recursion rule for w(x)=e^(-x^2)
  # and change of variable
  #
  # Apply periodic conditions on domain with length_scale length
  println("GH: [mean,variance]=",mu_sig_sq)
  nn=10*m
  v =sqrt((1:nn-1)/2.0);
  J=diagm(v,1)+diagm(v,-1)
  D,V=eig(J)
  xw=[D sqrt(pi)*V[1,:]'.^2]
  # change of variable
  f=x->mu_sig_sq[1]+x*sqrt(2*mu_sig_sq[2])
  weight_scale=sqrt(pi)
  quad=QUAD(nn,f(xw[:,1]),xw[:,2]/weight_scale)
  #
  periodic_simplify!(quad, m, length_scale)
  #
  println(quad.length,"-pt Gauss-Hermite rule.")
  sortquad!(quad)
  return quad
end
################################
function Gauss_Hermite_sech(mu_nu, n)
  # Returns m-point Gauss-Hermite quadrature rule
  # for sech^2( (x-mu)/(4 nu))/(8 nu)
  # for mu_nu=[mu,nu]
  # via Gauss-Hermite
  #
  nu=mu_nu[2]
  mu=mu_nu[1]
  s=nu*sqrt(8)
  mu_sig_sq=[mu, s^2]
  quad=Gauss_Hermite(mu_sig_sq, n)
  c=4*sqrt(pi*nu^2);
  helper=z-> c/(4*nu* (exp(-0.25*z/nu)+exp(0.25*z/nu))*exp(-(0.25*z/nu)^2) )
  for i=1:quad.length
    quad.w[i]=quad.w[i]*helper(quad.X[i]-mu)
  end
  quad.w=quad.w/sum(quad.w)
  return quad
end
#############################
function delta_measure(x)
  return  QUAD(1,x,1)
end
######################
function quad_to_pdf(quad::QUAD,n,s_in)
  sortquad!(quad)
  cdf_all=cumsum(quad.w)
  kk=3
  cdf_all=[zeros(kk);cdf_all;ones(kk)]
  eps=1e-2
  X=[quad.X[1]+linspace(-kk*eps,-eps,kk);quad.X;quad.X[end]+linspace(eps,kk*eps,kk)]
  n=min(quad.length,n)
  k=3
  range=collect(linspace(quad.X[1],quad.X[end],n))
  #print(length(X)," ", length(cdf_all), " ", length(range))
  #spl=Spline1D(X,cdf_all,range;w=ones(length(X)), k=k, bc="nearest")
  spl=Spline1D(X,cdf_all;w=ones(length(X)), k=k, bc="nearest",s=s_in)
  pdf=derivative(spl,range)
  cdf=spl(range)
  return range,pdf,cdf,X,cdf_all
end

######################
function quad_to_pdf(quad::QUAD)
  sortquad!(quad)
  cdf_all=cumsum(quad.w)
  kk=3
  cdf_all=[zeros(kk);cdf_all;ones(kk)]
  eps=1e-2
  X=[quad.X[1]+linspace(-kk*eps,-eps,kk);quad.X;quad.X[end]+linspace(eps,kk*eps,kk)]
  n=min(quad.length,n)
  k=3
  range=collect(linspace(quad.X[1],quad.X[end],n))
  #print(length(X)," ", length(cdf_all), " ", length(range))
  #spl=Spline1D(X,cdf_all,range;w=ones(length(X)), k=k, bc="nearest")
  spl=Spline1D(X,cdf_all)
  pdf=derivative(spl,range)
  cdf=spl(range)
  return range,pdf,cdf,X,cdf_all
end
#
end # module defn

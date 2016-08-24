__precompile__()
module iteratedintegral
export init, get_ito, get_strat
#################################################################
function init(fcn, opt)
  #=
  fcn=SDE definition (for dimension m and noise type)
  opt= options (for max_p)
  =#
  S=Dict("KL_C"=>1, "m"=>fcn.m, "NoiseType"=>fcn.noise_type)
  try
    S["max_p"]=opt["max_p"]
  catch
    S["max_p"]=1000
  end
  # get tail of sum_{p=P]^infty 1/p^2
  no_terms=100
  H2_tail=Array(Float16,no_terms)
  H2_tail[1]=pi^2/6 # zeta(2)
  for i=2:no_terms
    H2_tail[i]=H2_tail[i-1]-1/(i^2)
  end
  #S["H2_tail"]=H2_tail
  # map for upper-triangular to vector
  T=Int16
  m=fcn.m
  m2=round(T,m*(m-1)/2)
  ri=Array(T,m2)
  ci=Array(T,m2)
  k=1
  for i = 1:(m-1)
      for j = (i+1):m
        ri[k]=i
        ci[k]=j
        k=k+1
      end
  end
  S=merge(S,Dict("ri"=>ri, "ci"=>ci, "m2"=>m2))
  # initialise storage
  S["tmp_m2"]=Array(Float64,(m,m))
  S["tail"]=Array(Float64,(m,m))
  S["J"]=Array(Float64,m)
  S["K"]=Array(Float64,(m,m))
  #
  return S
end
#####
function __init__()
end
######################
function get_p(dt,S)
  #=
  get number of terms for a KL expansion in terms of dt
  and KL_C
  (not actually used as tail correction is implemented here)
  =#
  thresh=2*pi^2*S["KL_C"]^2 * dt
  p=findfirst(x->x<thresh,S["H2_tail"])
  if p<S["max_p"]
      p= 1 / (2 * (pi* S["KL_C"])^2 * dt)
  end
  return round(Int,p)
end
####################################
function Gamma_1!(out,b,S)
  #=
  helper function: see notes \Gamma_1
  in: b m2-vector
      S dictionary (see init function for defs of m and m2 and ri,ci)
  out: out (m x m)-matrix
  =#
  r=S["ri"]
  c=S["ci"]
  for i in 1:length(b)
    out[r[i],c[i]]=b[i]
    out[c[i],r[i]]=-b[i]
  end
end
######################################
function Gamma_2!(out,A,dW,S)
  #=
  helper function: see notes
  in: A (m x m)-matrix
  dW m-vector
  S dictionary (see init function)
  out: out (m x m)-matrix
  =#
  for i=1:length(dW)
    for j=1:length(dW)
      out[i,j]=vecdot(A[i,:],dW)*dW[j]
    end
  end
end
########################################
function Gamma_3!(out,A,S)
  #=
  helper function: see notes
  in: A (m x m)-matrix
  S dictionary (see init function)
  out: m2-vector
  =#
  r=S["ri"]
  c=S["ci"]
  for i in length(r)
    out[i]=A[r[i],c[i]] - A[c[i],r[i]]
  end
end
####################################
function compute_tail!(tail,dt,dW,ap,S)
  #=
  Tail approximation for KL expansion of stochastic integrals
  in: dt timestep
  dW m-vector
  ap coefficient (see notes)
  S dictionary (see init function)
  out: tail (m x m)-matrix
  =#
  scale=dt^2*ap/(2*pi^2)
  alpha=sqrt(1+sum(dW.*dW))
  noiseVec=randn(S["m2"])*sqrt(scale)
  tmp_M=noiseVec/((1+alpha)*dt)
  tmp_m2=S["tmp_m2"]
  Gamma_1!(tmp_m2,tmp_M,S)
  Gamma_2!(tail,tmp_m2,dW,S)
  Gamma_3!(tmp_M,tail,S)
  Gamma_1!(tail,tmp_M+noiseVec,S)
end
#
function get_strat!(dt,S)
  #=
  Main routine for generating integrals (Stratonvich)
  in: dt timestep
  S (see init function)
  out: J,K updated in S
  =#
  m=S["m"]
  J=S["J"]
  K=S["K"]
  if S["NoiseType"]==2 # diagonal
    J[:]=sqrt(dt)*randn(m)
    fill!(K,0)
    for i=1:m
      K[i,i]=0.5*J[i]^2
    end
  elseif S["NoiseType"]==3 # commutative
    J[:]=sqrt(dt)*randn(m)
    for i=1:m
      for j=1:m
       K[i,j]=0.5*J[i]*J[j]
    end
  end
  else # unstructured
    BBT!(dt,S)
  end
end
######################################
function get_ito!(dt,S)
  #=
  Main routine for generating integrals (Ito)
  in: dt timestep
  S (see init function)
  out: J,K updated in S
  =#
  get_strat!(dt,S)
  K=S["K"]
  for i=1:S["m"]
    K[i,i]-=0.5*dt
  end
end
########################################
function BBT!(dt,S)
  #=
  Routine for generating integrals (Stratonovich) unstructured case
  in: dt timestpe
  S dictionary (see init function)
  out: J,K updated in S
  =#
  m=S["m"]
  J=S["J"]
  K=S["K"]
  #
  p_const=sqrt(m*(m-1)/(24*dt)/(pi*S["KL_C"]))
  a_p=pi^2/6
  J[:]=sqrt(dt)*randn(m)
  for i=1:m
    K[i,i]=0.5*J[i]^2
  end
  p=round(p_const*sqrt(m*4*sum(J.*J)/dt))
  for r=1:(p+1)
    a_p=a_p-(1/r^2)
    sd_sq=dt/(2*pi*r)
    sd=sqrt(sd_sq)
    eta=randn(m)*sd
    zeta=randn(m)*sd
    K[:,1]-= zeta*(2/sqrt(pi*r))
    for i=2:m
      for j=1:i
        K[i,j]+=zeta[i]*eta[j]-eta[i]*zeta[j]
      end
    end
  end
  for i=2:m
    for j=1:i
      tmp1=K[i,j]-0.5*(K[j,1]*J[i]-K[i,1]*J[j])
      tmp2=0.5*J[i]*J[j]
      K[i,j]=tmp2+tmp1
      K[j,i]=tmp2-tmp1
    end
  end
  tail=S["tail"]
  compute_tail!(tail,dt,J,a_p,S)
  K+=tail
end
end

# input d>=2, eta>0 (eta=1 for uniform)
# output correlation matrix rr[][], density proportional to
#   det(R)^{eta-1}
#' @importFrom stats rbeta
rcorcvine<-function(d,eta=1)
{
  d<-as.integer(d)
  if(d<=0 || !is.integer(d))
  { stop("The dimension 'd' should be a positive integer!\n") }
  if(eta<=0)
  { stop("'eta' should be positive!\n") }

  #handling of d=1 and d=2
  if(d==1)
  { rr<-matrix(1,1,1); return(rr) }
  if(d==2)
  { rho<-2*rbeta(1,eta,eta)-1
    rr<-matrix(c(1,rho,rho,1),2,2); return(rr)
  }
  rr<-matrix(0,d,d)
  # matrix of partial correlations as generated
  prr<-matrix(0,d,d)
  diag(rr)<-1
  for(i in 2:d)
  { alp<-eta+(d-2)/2
    rr[1,i]<-2*rbeta(1,alp,alp)-1
    rr[i,1]<-rr[1,i]
    prr[1,i]<-rr[1,i]
  }
  for(m in 2:(d-1))
  { alp<-eta+(d-1-m)/2
    for(i in (m+1):d)
    { prr[m,i]<-2*rbeta(1,alp,alp)-1
      # back calculate thru lower order partials
      tem<-prr[m,i]
      for(k in (m-1):1)
      { tem<-prr[k,m]*prr[k,i]+tem*sqrt((1-prr[k,m]^2)*(1-prr[k,i]^2)) }
      rr[m,i]<-tem
      rr[i,m]<-rr[m,i]
    }
  }
  return(rr)
}

# input d>=2, eta>0 (eta=1 for uniform)
# output correlation matrix rr[][], density proportional to
#   det(R)^{eta-1}
#' @importFrom stats rbeta rnorm
rcoronion<-function(d,eta=1)
{
  d<-as.integer(d)
  if(d<=0 || !is.integer(d))
  { stop("The dimension 'd' should be a positive integer!\n") }
  if(eta<=0)
  { stop("'eta' should be positive!\n") }

  #handling of d=1 and d=2
  if(d==1)
  { rr<-matrix(1,1,1); return(rr) }
  if(d==2)
  { rho<-2*rbeta(1,eta,eta)-1
    rr<-matrix(c(1,rho,rho,1),2,2); return(rr)
  }
  rr<-matrix(0,d,d)
  beta<-eta+(d-2)/2
  # step 1
  r12<-2*rbeta(1,beta,beta)-1
  rr<-matrix(c(1,r12,r12,1),2,2)
  # iterative steps
  for(m in 2:(d-1))
  { beta<-beta-0.5
    y<-rbeta(1,m/2,beta)
    z<-rnorm(m,0,1)
    znorm<-sqrt(sum(z^2))
    # random on surface of unit sphere
    z<-z/znorm
    w=sqrt(y)*z
    # can spped up by programming incremental Cholesky?
    rhalf<-chol(rr)
    qq<-w%*%rhalf
    rr<-cbind(rr,t(qq))
    rr<-rbind(rr,c(qq,1))
  }
  # return rr
  rr
}

#' @importFrom matrixcalc is.positive.definite
rcorcvine.wrap<-function(nmat=1, d, eta=1, onlypos=FALSE){
  r.mat<-matrix(NA,d*d,0)
  for(r.indx in nmat){
    r<-rcorcvine(d,eta)
    if (onlypos) {
      r = (r+1)/2
    }
    r.mat<-cbind(r.mat,c(r))
    while(!is.positive.definite(matrix(r,d,d))){
      r<-c(rcorcvine(d,eta))
      if (onlypos) {
        r = (r+1)/2
      }
    }
  }
  return(r.mat)
}

#' @importFrom matrixcalc is.positive.definite
rcoronion.wrap<-function(nmat=1, d, eta=1, onlypos=FALSE){
  r.mat<-matrix(NA,d*d,0)
  for(r.indx in nmat){
    r<-rcoronion(d,eta)
    if (onlypos) {
      r = (r+1)/2
    }
    while(!is.positive.definite(matrix(r,d,d))){
      r<-c(rcoronion(d,eta))
      if (onlypos) {
        r = (r+1)/2
      }
    }
    r.mat<-cbind(r.mat,c(r))
  }
  return(r.mat)
}

mcmc.jump.defaults<-function(o){
  if (o==3) {
    jmpsize=.56
  } else if (o==4) {
    jmpsize=.48
  } else if (o==5) {
    jmpsize=.4
  } else if (o==6) {
    jmpsize=.34
  } else if (o==7) {
    jmpsize=.3
  } else if (o==8) {
    jmpsize=.26
  } else if (o==9) {
    jmpsize=.23
  } else if (o==10) {
    jmpsize=.21
  } else if (o==11) {
    jmpsize=.19
  } else if (o==12 & o <15) {
    jmpsize=.175
  } else if (o==15) {
    jmpsize=.12
  } else if (o>=16) {
    jmpsize=.1
  }

  return(jmpsize)
}

#' @importFrom stats runif rnorm
mcmc <- function(nmat = 1, d, eta=1, reject=NULL, onlypos=NULL) {

  eta_init<-eta
  thin<-floor(eta/length(nmat))

  r.mat<-matrix(NA,d*d,0)
  #eta here is the number of iterations to go through before outputting the matrix (i.e., thinning)
  o = d
  r = matrix(0,o,o)

  jmpsize<-mcmc.jump.defaults(o)

  if(is.null(reject)){
    reject <- 0

    nr<-as.integer(o*(o-1)/2)
    if(onlypos){
      xcand <- rep(.5,times=nr)
    } else {
      xcand <- rep(0,times=nr)
    }
    xtmp <- rep(0,times=nr)
  }

  while(eta>=0){
    tmp <- rnorm(nr)
    xtemp <- runif(1)
    xtmp<-(xtemp^(1/nr))*tmp/as.vector(sqrt(crossprod(tmp,tmp)))

    ycand <- xcand + jmpsize*xtmp

    # this is the only place I see where positive correlations are enforced, and this seems unnecessary
    if(onlypos){
      yflag<-ifelse(all(ycand>0 & ycand<1), 0, 1)
    } else {
      yflag<-ifelse(all(ycand>-1 & ycand<1), 0, 1)
    }
    if (yflag==1) {
      reject<-reject+1 # this used for anything?
    } else {
    #if (yflag!=1) {
      i=1 # put contents of y in l.t. of r (columnwise), check eigenvalues.
      for (k in 1:(o-1)){
        for (j in (k+1):o) {
          r[j,k]=ycand[i]
          r[k,j]=ycand[i]
          i=i+1
        }
      }
      diag(r)<-1

      ev<-eigen(r)$values
      #if(!(any(ev<.08))){
      if(!(any(ev<.01))){
      #if (is.positive.definite(r)) {

        xcand<-ycand
        eta<-eta-1
        # multiple of thinning
        if((eta %% thin & eta < eta_init)==0){
          r.mat<-cbind(r.mat,c(r))
        }
      }
    }
  }
  return(r.mat)
}

# Helper function that generates random matrices
#' @importFrom matrixcalc is.positive.definite
genmat<-function(d,eta_val,onlypos,mat_func){

  temp=c(mat_func(d,eta_val))
  if (onlypos) {
    temp = (temp+1)/2
  }

  while(!is.positive.definite(matrix(temp,d,d))){
    temp=c(mat_func(d,eta_val))
    if (onlypos) {
      temp = (temp+1)/2
    }
  }
  temp
}

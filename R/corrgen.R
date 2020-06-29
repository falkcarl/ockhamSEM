# Functions for generation of random correlation matrices

mcmc.args.default<-function(nchains, d, args=NULL){

  control<-list(
    iter = 5000000, # total number of iterations to run
    miniter = 10000, # min number of iterations per chain
    dim = d, # dimension of covariance matrix
    jmpsize = mcmc.jump.defaults(d)
  )

  # override if any custom arguments provided
  control[names(args)] <- args

  # divide iterations per chain
  control$iter<-max(ceiling(control$iter/nchains),control$miniter)

  return(control)
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
  } else if (o>=12 & o <15) {
    jmpsize=.175
  } else if (o==15) {
    jmpsize=.12
  } else if (o>=16) {
    jmpsize=.1
  }

  return(jmpsize)
}

onion.args.default<-function(d,args=NULL){

  control<-list(
    dim = d, # dimension of covariance matrix
    covMethod = "onion",
    eta = 1, # uniform
    rangeVar = c(1,1) # unit variance per variable
  )

  # override if any custom arguments provided
  control[names(args)] <- args

  return(control)

}

clustergen.args.default<-function(d,args=NULL){

  control<-list(
    dim = d # dimension of covariance matrix
  )

  # override if any custom arguments provided
  control[names(args)] <- args

  return(control)
}

#' @importFrom stats runif rnorm
mcmc <- function(nmat, dim, iter, jmpsize, reject=NULL, onlypos=NULL) {

  iter_init<-iter
  thin<-floor(iter/length(nmat))

  r.mat<-matrix(NA,dim*dim,0)

  o = dim
  r = matrix(0,o,o)

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

  while(iter>=0){
    tmp <- rnorm(nr)
    xtemp <- runif(1)
    xtmp<-(xtemp^(1/nr))*tmp/as.vector(sqrt(crossprod(tmp,tmp)))

    ycand <- xcand + jmpsize*xtmp

    # this is the only place I see where positive correlations are enforced
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
      # other ways to check positive definiteness that have been tried:
      #if(!(any(ev<.08))){
      #if (is.positive.definite(r)) {
      if(!(any(ev<.01))){

        xcand<-ycand
        iter<-iter-1
        # multiple of thinning
        if((iter %% thin & iter < iter_init)==0){
          r.mat<-cbind(r.mat,c(r))
        }
      }
    }
  }
  return(r.mat)
}

# Helper function that generates random matrices
#' @importFrom matrixcalc is.positive.definite is.symmetric.matrix
#' @importFrom clusterGeneration genPositiveDefMat
genmat<-function(nmat=1, rmethod=c("mcmc","onion","clustergen"), control, onlypos=FALSE){

  rmethod <- match.arg(rmethod)

  if(rmethod=="mcmc"){
    r.mat <- mcmc(nmat, control$dim, control$iter, control$jmpsize, onlypos=onlypos)
  } else if (rmethod=="onion" | rmethod=="clustergen"){
    r.mat<-matrix(NA,control$dim*control$dim,0)
    for(r.indx in nmat){
      r <- do.call("genPositiveDefMat",control)$Sigma
      if (onlypos) {
        r = (r+1)/2 # ad-hoc correction to ensure positive manifold
      }
      if(!is.symmetric.matrix(r)){
        r<-round(r,5) # ad-hoc fix
      }
      while(!is.positive.definite(r)){
        r<-c(do.call("genPositiveDefMat",control)$Sigma)
        if (onlypos) {
          r = (r+1)/2 # ad-hoc correction to ensure positive manifold
        }
        if(!is.symmetric.matrix(r)){
          r<-round(r,5) # ad-hoc fix
        }
      }
      r.mat<-cbind(r.mat,as.vector(r))
    }
  }
  return(r.mat)
}

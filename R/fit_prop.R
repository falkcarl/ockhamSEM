# input d>=2, eta>0 (eta=1 for uniform)
# output correlation matrix rr[][], density proportional to
#   det(R)^{eta-1}
#' @importFrom stats rbeta
rcorcvine<-function(d,eta=1,nmat=1)
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

#' @importFrom stats runif rnorm
mcmc <- function(nmat = 1, d, eta=1, reject=NULL, onlypos=NULL) {

  eta_init<-eta
  thin<-floor(eta/length(nmat))

  r.mat<-matrix(NA,d*d,0)
  #eta here is the number of iterations to go through before outputting the matrix (i.e., thinning)
  o = d
  r = matrix(0,o,o)

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

#' Run fit propensity analyses
#'
#' @param ... models to be fit to the same data of class lavaan
#' @param fit.measure character string or vector that indicates which fit measure to extract from fitted models
#' @param saveModel logical value indicating whether the to save fitted models for later examination
#' @param saveR logical value indicating whether to save random correlation matrices
#' @param control list of parameters for generation of random correlation matrices
#' @param cluster a cluster created by makeCluster() from the parallel package
#' @return An object of class fitprop
#' @export
#' @importFrom lavaan lavInspect lavNames parTable lavaan
#' @importFrom parallel clusterSetRNGStream clusterSplit parLapply
fit_prop <- function(..., fit.measure=c("srmr"),
                      saveModel=FALSE, saveR=FALSE,
                      control=list(eta=1, reps=1000,rmethod="onion", onlypos=FALSE, seed=1234),
                      cluster=NULL) {

  fun.call<-match.call()
  models <- list(...)
  n.models <- length(models)

  for(lavmodel in models){
    if(!class(lavmodel)=="lavaan"){
      stop("Input lavmodel must be a fitted lavaan model")
    }

    # Check for multiple group model
    if(lavInspect(lavmodel,"ngroups")>1){
      stop("Only single group models are currently supported")
    }
  }

  # Extract variable names
  vnames<-lavNames(models[[1]])

  # number of variables
  d<-length(vnames)

  # number of fit measures
  n.fit<-length(fit.measure)

  # check fit.measure
  if(!is.character(fit.measure)){
    stop("fit.measure should be a character string or vector indicating the measure to be extracted from fitMeasures()")
  }

  #This bit is only for mcmc generation method
  o=d
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

  # check control
  if(is.null(control[["eta"]])){control[["eta"]]<-1}
  if(is.null(control[["reps"]])){control[["reps"]]<-5000}
  if(is.null(control[["rmethod"]])){control[["rmethod"]]<-"onion"}
  if(is.null(control[["onlypos"]])){control[["onlypos"]]<-FALSE}
  if(is.null(control[["seed"]])){control[["seed"]]<-1234}
  if(is.null(control[["iter"]])){control[["iter"]]<-5000000}
  if(is.null(control[["miniter"]])){control[["miniter"]]<-10000}
  onlypos = as.logical(control["onlypos"])

  #Set matrix for storing outputs
  nsim = as.numeric(control["reps"])
  fit_list = vector("list", n.models)
  for (i in 1:length(fit_list)) {
    fit_list[[i]] = matrix(0,nsim,length(fit.measure))
    colnames(fit_list[[i]]) <- fit.measure
  }

  eta_val = as.numeric(control["eta"])

  mat_func = rcoronion.wrap
  if(control["rmethod"]=="onion") {
    mat_func = rcoronion.wrap
  } else if (control["rmethod"]=="vine") {
    mat_func = rcorcvine.wrap
  } else if (control["rmethod"]=="mcmc") {
    mat_func = mcmc
    if(!is.null(cluster)){
      eta_val <- max(ceiling(control$iter/length(cluster)),control$miniter)
    } else {
      eta_val <- control$iter
    }

  } else {
    stop(paste('Bad matrix function:', control["rmethod"]))
  }

  #Generate the matrices
  print("Generate matrices")
  if(!is.null(cluster)){
    clusterSetRNGStream(cluster, control$seed)
    nmat<-clusterSplit(cluster, 1:nsim)
    out_mat<-parLapply(cluster,nmat,mat_func,d=d,eta=eta_val,onlypos=onlypos)
  } else {
    set.seed(control$seed)
    out_mat<-lapply(list(1:nsim),mat_func,d=d,eta=eta_val,onlypos=onlypos)
  }
  out_mat<-do.call(cbind,out_mat)
  row.names(out_mat)<-NULL

  # Do elements of @Data and @SampleStats depend on whether data is raw or cov matrix is used?
  # If so, pre-fit each model to cov matrix once
  # If not, this code isn't necessary
  j<-1
  for (lavmodel in models) {
    temp_matrix<-lavInspect(lavmodel,"sampstat")$cov
    ptab<-parTable(lavmodel)
    opt<-lavmodel@Options
    opt$se<-"none"
    nobs<-lavInspect(lavmodel,"nobs")
    my_fitted_model<-lavaan(model=ptab,
                            sample.cov=temp_matrix,
                            #se="none",
                            sample.nobs= nobs,
                            slotOptions = opt)
                            #slotModel = lavmodel@Model)
    #my_fitted_model<-update(lavmodel, sample.cov=temp_matrix, se="none")
    models[[j]]<-my_fitted_model
    j<-j+1
  }

  #Fit model to each matrix
  j = 0
  print("Fitting models")
  mod_list<-list()
  for (lavmodel in models) {
    j = j + 1
    if(!is.null(cluster)){
      #clusterEvalQ(cluster,require(lavaan))
      fit_tmp<-parLapply(cluster,1:nsim,fitmod,lavmodel=lavmodel,out_mat=out_mat,vnames=vnames,j=j,saveModel=saveModel,d=d,fit.measure=fit.measure)
    } else {
      fit_tmp<-lapply(1:nsim,fitmod,lavmodel=lavmodel,out_mat=out_mat,vnames=vnames,j=j,saveModel=saveModel,d=d,fit.measure=fit.measure)
    }
    fit_list[[j]]<-matrix(unlist(fit_tmp),ncol=length(fit.measure),byrow=TRUE)
    colnames(fit_list[[j]])<-fit.measure

    if(saveModel){
      mod_list[[j]]<-lapply(fit_tmp,attr,"mod")
    }
  }

  # Process number of missing values (if any) per model and fit index
  na_list<-list()
  complete_mat<-matrix(1,nrow=control[["reps"]],ncol=n.fit)
  for(j in 1:length(fit_list)){
    na_list[[j]]<-is.na(fit_list[[j]])
    complete_mat[which(is.na(fit_list[[j]]))]<-NA
  }

  out<-list()
  out$fun.call<-fun.call
  out$eta=control$eta
  out$reps=control$reps
  out$rmethod=control$rmethod
  out$onlypos=control$onlypos
  out$nfit=n.fit

  out$na_list<-na_list
  out$complete_mat<-complete_mat

  out$fit_list<-fit_list

  if(saveR){
    out$R<-list()
    for(j in 1:nsim){
      out$R[[j]]<-matrix(out_mat[,j],d,d)
    }
  }

  if(saveModel){
    out$mod_list<-mod_list
  }

  class(out)<-"fitprop"
  return(out)
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

#' @importFrom lavaan lavaan lavInspect fitMeasures
fitmod <- function(indx,lavmodel,out_mat,vnames,j,saveModel,d,fit.measure){
  temp_matrix <- matrix(out_mat[,indx],d,d)
  colnames(temp_matrix) <- rownames(temp_matrix) <- vnames
  # added from ShortForm Tabu code; nuke all starting values? Does this work? Not with older lavaan version
  lavmodel@ParTable$est<-NULL
  lavmodel@ParTable$start<-NULL

  # replaced with update(); might be slower, but more stable? (doesn't work with older lavaan version)
  #res<-try({my_fitted_model<-update(lavmodel, sample.cov=temp_matrix)
  res<-try ({
    lavmodel@Options$se="none"
    lavmodel@Options$start="default"
    lavmodel@Options$do.fit=TRUE
    #browser()
  # lavmodel@Data
  # lavmodel@SampleStats
    #lavmodel@SampleStats@cov[[1]]<-temp_matrix
    #lavmodel@SampleStats@icov[[1]]<-solve(temp_matrix)
    #lavmodel@SampleStats@WLS.obs[[1]]<-lav_matrix_vech(temp_matrix)

    my_fitted_model<-lavaan(sample.cov=temp_matrix,
                                    sample.nobs=lavInspect(lavmodel,"nobs"),
                            #slotData = lavmodel@Data,
                            #slotSampleStats = lavmodel@SampleStats,
                            #slotModel = lavmodel@Model,
                            slotOptions = lavmodel@Options,
                            slotParTable = lavmodel@ParTable,
                            slotCache = lavmodel@Cache)

    if(my_fitted_model@optim$converged){
      #Store fit values
      fit_out <- fitMeasures(my_fitted_model, fit.measure)
    } else {
      fit_out <-rep(NA,length(fit.measure))
    }

    # regardless, save fitted lavaan model for later inspection
    if(saveModel) {
      attr(fit_out,"model")<-my_fitted_model
    }

  }, silent=T)

  if(class(res)!="try-error"){
    return(fit_out)
  } else {
    return(rep(NA,length(fit.measure)))
  }
}

#' Plot function for fitprop objects
#' @param x Object of class fitprop, as created by fit_prop function.
#' @param ... Does nothing but to hopefully make this generic function pass R CMD check
#' @param type What type of plot to produce?
#' @param whichmod Index number corresponding to which model(s) to include on the plot
#' @param whichfit Character vector indicating which indices of model fit to include
#' @param savePlot Logical value indicating whether to save plot to a list (TRUE) or just produce plot to output.
#' @param xlim Vector of length 2 indicating the limits for the x-axis of the plot
#' @param samereps Logical value indicating whether to use only results from replications in which all selected models yielded results
#' @param cutoff Numeric vector indicating what cut value of the fit indice(s) to use for euler plots
#' @param lower.tail Logical vector indicating whether lower values of each fit index corresponds to good fit
#' @param mod.lab Optional character vector of labels for each model
#' @param mod.brewer.pal Optional character corresponding to the palette from RColorBrewer to use for the different models
#' @export
#' @importFrom ggplot2 ggplot aes stat_ecdf scale_color_brewer theme element_blank xlab ylab xlim ylim
#' @importFrom stats na.omit
#' @importFrom tidyr gather
#' @importFrom eulerr euler
#' @importFrom nVennR plotVenn
#' @importFrom rlang .data
plot.fitprop<-function(x,...,type="ecdf",whichmod=NULL,whichfit=colnames(x$fit_list[[1]]),savePlot=FALSE,
                       xlim=c(0,1),samereps=TRUE,cutoff=rep(.1,length(whichfit)),lower.tail=rep(TRUE,length(whichfit)),
                       mod.lab=NULL,mod.brewer.pal="Set1"){

  # TO DO: probably match.args for "type", some other input checking

  data<-x$fit_list
  nmod<-length(data) # number of models
  nfit<-ncol(data[[1]]) # number of fit measures
  nrep<-nrow(data[[1]]) # number of available replications

  if(is.null(whichfit)){
    whichfit<-colnames(data[[1]])
  }
  if(is.null(whichmod)){
    whichmod<-1:nmod
  }
  if(is.null(mod.lab)){
    mod.lab<-paste0("Model ",whichmod)
  }

  plots<-list()
  # loop over fit indices
  m<-1
  for(fm in whichfit){

    # extract and format data
    dat<-matrix(,nrep,nmod)
    j<-1
    for(mod in 1:nmod){
      dat[,j]<-data[[mod]][,fm]
      j<-j+1
    }
    dat<-as.data.frame(dat)
    colnames(dat)<-mod.lab
    dat$id<-1:nrow(dat)

    # generate plots
    if(type=="ecdf"){
      #dat<-melt(dat,na.rm=TRUE,id.vars="id")
      if(samereps){
        dat<-na.omit(dat)
      }
      dat<-gather(dat,"variable","value",whichmod)
      graph<-ggplot(dat,aes(x=.data$value))+
        stat_ecdf(aes(linetype=.data$variable, color=.data$variable),na.rm=TRUE,size=.7) +
        scale_color_brewer(palette=mod.brewer.pal)

      if(lower.tail[m]){
        graph<-graph+ xlim(xlim[1],xlim[2])
      } else {
        graph<-graph+ xlim(xlim[2],xlim[1])
      }
      graph<-graph+theme(legend.title=element_blank())+
        xlab(fm) + ylab("Cumulative Probability")
    } else if (type=="euler"){
      if(lower.tail[m]){
        dat[,whichmod]<-dat[,whichmod]<cutoff # how many meet cutoff criterion
      } else {
        dat[,whichmod]<-dat[,whichmod]>cutoff # how many meet cutoff criterion
      }
      tmp<-cbind(data.frame(Total=rep(TRUE,nrow(dat))),dat[,whichmod])
      if(samereps){
        tmp<-na.omit(tmp)
      }
      colnames(tmp)<-c("Total",whichmod)
      eulerfit<-euler(tmp)
      graph<-plot(eulerfit)
    } else if (type=="nVennR"){
      #stop("Not yet functional")
      #browser()
      dat<-na.omit(dat)
      tmp<-list()
      indx<-1
      for(j in whichmod){
        if(lower.tail[m]){
          tmp[[indx]]<-dat$id[dat[,j]<cutoff]
        } else {
          tmp[[indx]]<-dat$id[dat[,j]>cutoff]
        }
        indx<-indx+1
      }
      tmp[[indx]]<-dat$id # last group is the total
      #tmp<-cbind(data.frame(Total=rep(TRUE,nrow(dat))),dat[,whichmod])
      #if(samereps){
      #  tmp<-na.omit(tmp)
      #}
      #colnames(tmp)<-c("Total",whichmod)
      #eulerfit<-euler(tmp)
      #graph<-plot(eulerfit)
      graph<-plotVenn(tmp,nCycles=5000,showPlot=F,sNames=c(mod.lab,"Total"))
    }

    # do something with plot
    if(savePlot){
      plots[[m]]<-graph
    } else {
      if(m>1){
        invisible(readline(prompt="Press [enter] to continue"))
      }
      print(graph)

    }
    m<-m+1
  }
  if(savePlot){
    invisible(plots)
  }
}

#' @export
print.fitprop<-function(x,...){
  data<-x$fit_list
  nmod<-length(data) # number of models
  nfit<-ncol(data[[1]]) # number of fit measures
  nrep<-x$reps # number of available replications

  cat(
"\n",
#"Function call =           ", deparse(x$fun.call), "\n",
"Number of fitted models = ", nmod, "\n",
"Number of fit measures =  ", nfit, "\n",
"Fit measures =            ", colnames(data[[1]]), "\n",
"Number of replications =   ", nrep, "\n",
"Options for R generation\n",
"  method =                ", x$rmethod, "\n",
"  eta =                   ", x$eta, "\n",
"  Only positive R?        ", x$onlypos, "\n"
)

}

#' Summary function for fitprop objects
#' @param object Object of class fitprop, as created by fit_prop function.
#' @param ... Does nothing but to hopefully make this generic function pass R CMD check
#' @param probs passed to quantile to determine what probabilities to report
#' @param samereps Logical value indicating whether to use only results from replications in which all selected models yielded results
#' @param lower.tail Logical vector indicating whether lower values of each fit index corresponds to good fit
#' @export
#' @importFrom stats quantile median ks.test na.omit
#' @importFrom utils str
#' @importFrom effsize cliff.delta cohen.d
summary.fitprop<-function(object,...,probs=seq(0,1,.1),samereps=TRUE,lower.tail=rep(TRUE,ncol(object$fit_list[[1]]))){

  data<-object$fit_list
  nmod<-length(data) # number of models
  nfit<-ncol(data[[1]]) # number of fit measures
  nrep<-object$reps # number of available replications

  stats<-list()
  quantiles<-list()
  efs<-list()

  cat(
    "\n",
    "Quantiles for each model and fit measure:\n"
  )
  for(mod in 1:nmod){
    cat("\n Model ",mod,"\n")
    if(samereps){
      qtmp<-NULL
      for(j in 1:nfit){
        qtmp<-cbind(qtmp, quantile(data[[mod]][!is.na(object$complete_mat[,j]),j],...,probs=probs,na.rm=TRUE))
      }
      colnames(qtmp)<-colnames(data[[mod]])
    } else {
      qtmp<-apply(data[[mod]],2,quantile,...,probs=probs,na.rm=TRUE)
    }
    for(m in 1:nfit){
      qtmp[,m]<-sort(qtmp[,m],decreasing=!lower.tail[m])
    }
    quantiles[[mod]]<-qtmp
    print.default(round(qtmp,3))
    #printCoefmat(qtmp,digits=3,dig.tst=3,cs.ind=1)
  }

  cat(
    "\n",
    "Information about replications for each model and fit measure:\n"
  )

  for(mod in 1:nmod){
    cat("\n Model ",mod,"\n")
    finite<-apply(data[[mod]],2,function(x){sum(is.finite(x))})
    na<-apply(data[[mod]],2,function(x){sum(is.na(x))})
    nrep<-nrow(data[[mod]])
    if(samereps){
      means<-NULL
      medians<-NULL
      for(j in 1:nfit){
        means<-c(means,mean(data[[mod]][!is.na(object$complete_mat[,j]),j]))
        medians<-c(medians,median(data[[mod]][!is.na(object$complete_mat[,j]),j]))
      }
      names(means)<-names(medians)<-colnames(data[[mod]])
    } else {
      means<-colMeans(data[[mod]],na.rm=TRUE)
      medians<-apply(data[[mod]],2,median,na.rm=TRUE)
    }

    stats[[mod]]<-list()
    stats[[mod]]$finite<-finite
    stats[[mod]]$finite<-na
    stats[[mod]]$means<-means
    stats[[mod]]$medians<-medians

    cat("\nMean across replications\n")
    print.default(round(means,3))
    cat("\nMedian across replications\n")
    print.default(round(medians,3))
    cat("\nNumber of finite values\n")
    print.default(finite)
    cat("\nNumber of NA values\n")
    print.default(na)
  }

  cat(
    "\n",
    "Effect Sizes for Differences in Model Fit:\n"
  )

  # K-S tests and other effect sizes
  # Compare all possible pairs of models, and all fit measures
  for(j in 1:nfit){

    efs[[j]]<-list()
    cat("\n ",colnames(data[[1]])[j],"\n")

    indx<-1
    for(mod in 1:(nmod-1)){
      for(mod2 in (mod+1):nmod){
        efs[[j]][[indx]]<-list()
        efs[[j]][[indx]]$title<-paste0("Model ",mod, " vs. Model ", mod2)
        if(samereps){
          tmp<-c(data[[mod]][!is.na(object$complete_mat[,j]),j],data[[mod2]][!is.na(object$complete_mat[,j]),j])
          grp.tmp<-c(rep(1,sum(!is.na(object$complete_mat[,j]))),rep(2,sum(!is.na(object$complete_mat[,j]))))

          efs[[j]][[indx]]$d<-cohen.d(tmp,grp.tmp)$estimate
          efs[[j]][[indx]]$delta<-cliff.delta(data[[mod]][!is.na(object$complete_mat[,j]),j],data[[mod2]][!is.na(object$complete_mat[,j]),j])$estimate
          efs[[j]][[indx]]$ks<-ks.test(data[[mod]][!is.na(object$complete_mat[,j]),j],data[[mod2]][!is.na(object$complete_mat[,j]),j])$statistic
        } else {
          tmp<-c(data[[mod]][,j],data[[mod2]][,j])
          tmp<-na.omit(tmp)
          grp.tmp<-c(rep(1,nrep),rep(2,nrep))
          grp.tmp<-grp.tmp[-attr(tmp,"na.action")]
          efs[[j]][[indx]]$d<-cohen.d(tmp,grp.tmp)$estimate

          efs[[j]][[indx]]$delta<-cliff.delta(data[[mod]][,j],data[[mod2]][,j])$estimate
          efs[[j]][[indx]]$ks<-str(ks.test(data[[mod]][,j],data[[mod2]][,j]))$statistic
        }

        cat("\n", efs[[j]][[indx]]$title, "\n")
        cat("   Cohen's d:          ", round(efs[[j]][[indx]]$d,3), "\n")
        cat("   Cliff's delta:      ", round(efs[[j]][[indx]]$delta,3), "\n")
        cat("   Komolgorov Smirnov: ", round(efs[[j]][[indx]]$ks,3), "\n")

        indx<-indx+1
      }
    }
  }

  out<-list()
  out[["quantiles"]]<-quantiles
  out[["stats"]]<-stats
  out[["efs"]]<-efs
  invisible(out)
}

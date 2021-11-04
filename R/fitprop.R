# Main UI functions

#' Run fit propensity analyses
#'
#' @param ... Models of class lavaan for which the user would like to compare fit propensity.
#' @param fit.measure Character vector that indicates which fit measure to extract from fitted
#'   models. Possible options include anything returned by \code{\link[lavaan]{fitMeasures}} from the lavaan
#'   package when applied to fitted models.
#' @param rmethod String indicating the type of random correlation generation approach.
#'   Choices are \code{"mcmc"} (default), \code{"onion"}, and \code{"clustergen"}. See details.
#' @param reps Number of random correlation matrices to generate for fit propensity analysis.
#' @param onlypos Logical value indicating whether to generate correlation matrices. Note that if there
#'   are many variables, generation of correlation matrices and fitting models to them will be very, very
#'   computationally intensive.
#'   with positive manifold (\code{TRUE}); i.e., only positive relationships among variables.
#' @param seed Random number seed used by set.seed or parallel package.
#' @param mcmc.args Named list of arguments that controls options for
#'   \code{"mcmc"} correlation matrix generation. See details.
#' @param clustergen.args Named list of arguments that controls generation of
#'   correlation matrices if \code{"onion"} or \code{"clustergen"} is used. See details.
#' @param saveModel Logical value indicating whether the to save fitted models
#'   for later examination.
#' @param saveR Logical value indicating whether to save randomly generated correlation matrices.
#' @param cluster (Optional) A cluster created by \code{\link[parallel]{makeCluster}} from the parallel package. If
#'   provided, computations will be parallelized as much as possible.
#' @details Inspired by work by Preacher (2003, 2006) and Bonifay & Cai (2017),
#'   this function performs three steps for analyses to assess the fit propensity of competing
#'   structural equation models: 1. Randomly generate correlation (or covariance matrices);
#'   2. Fit models to each correlation matrix; and 3. Save a indices that could be used for
#'   evaluating model fit in subsequent summaries. Conceptually, models that exhibit better fit
#'   to such randomly generated data may have better fit propensity, and are therefore potentially
#'   less parsimonious.
#'
#'   Analyses are performed with the \code{lavaan} package, and fitted lavaan models of
#'   \code{\link[lavaan]{lavaan-class}} (e.g., created from
#'   \code{\link[lavaan]{cfa}}, \code{\link[lavaan]{sem}}, or \code{\link[lavaan]{lavaan}} functions)
#'   for the competing models must be passed as initial arguments
#'   to the function. Currently, only single-group models and those relying on ML estimation are
#'   supported. Otherwise, the underlying options for the fitted lavaan models will be re-used by
#'   the run.fitprop function for the fit propensity analyses. It is optional to save the randomly
#'   generated matrices from Step 1 and the models fit in Step 2. Follow-up summaries of results
#'   saved in Step 3 are provided by plot.fitprop and summary.fitprop functions.
#'
#'   Generation of random correlation matrices is provided using several approaches. The \code{"mcmc"}
#'   algorithm implements a Markov Chain Monte Carlo approach and was ported from Fortran code
#'   in Preacher (2003). For details on the algorithm's actual implementation, see Preacher (2003),
#'   Falk and Muthukrishna (in prep), or the source code for the mcmc function. If this algorithm
#'   is chosen, \code{mcmc.args} accepts a list that can modify some default settings. In particular,
#'   \code{iter} sets the total number of iterations to run (default = 5000000). If parallel processing
#'   is enabled, this number will be divided amonst the number of chains. \code{miniter} sets a
#'   minimum number of iterations per chain to avoid many processors leading to too few iterations per
#'   chain (default = 10000). \code{jmpsize} overrides the step size for each update to the candidate
#'   correlation matrix. Smaller step sizes typically lead to more acceptance and may be necessary for
#'   larger correlation matrices (default jump size depends on the number of variables). Though, in
#'   general the MCMC algorithm becomes more difficult to work well with many variables.
#'
#'   The \code{"onion"} method is one approach that relies on work of Joe (2006) and
#'   Lewandowski, Kurowick, and Joe (2009); matrices are generated recursively, one variable at a
#'   time. The onion method is computationally more efficient than the MCMC algorithm. Under the
#'   hood, the \code{\link[clusterGeneration]{genPositiveDefMat}} function in the clusterGeneration package is used, with default
#'   arguments of \code{covMethod="onion"}, \code{eta=1}, and \code{rangeVar=c(1,1)}. These arguments ensure that the
#'   Onion method is used, generation is uniform over the space of positive definite matrices
#'   (but see note on positive manifold below), and with unit variances.
#'
#'   An additional option \code{"clustergen"} is provided for direct interface with the \code{\link[clusterGeneration]{genPositiveDefMat}}
#'   function in the clusterGeneration package. A named list can be passed to \code{clustergen.args} to
#'   override any defaults used by \code{\link[clusterGeneration]{genPositiveDefMat}}, and the user is referred to documentation
#'   for that function. This allows, for example, generation using C-Vines, covariance matrices
#'   (i.e., variables that do not all have unit variances), and several other covaraince/correlation
#'   matrix generation techniques.
#'
#'   onlypos controls whether correlation matrices can have only positive correlations.
#'   The original MCMC algorith by Preacher (2003, 2006) generated correlation matrices with
#'   positive manifold only (i.e., only positive correlations). The algorithm is easily changed
#'   to allow also negative correlations. The Onion method and any functions from clusterGeneration
#'   by default generate matrices with both positive and negative correlations. To obtain
#'   matrices with positive manifold only, an ad-hoc correction is implemented for these latter
#'   approaches where the matrix is transformed: R = (R+1)/2. To our knowledge, there is no
#'   guarantee that this will result in uniform sampling from the space of all correlation matrices
#'   with positive manifold, yet fit propensity results for some examples are very similar to those
#'   of the MCMC algorithm.
#'
#' @slot fit_list A list of the same length as the number of models being compared. Each list
#'   contains a matrix with columns corresponding to each entry of \code{fit.measure} and for all replications.
#' @slot R A list of the same length as the number of replications, containing all correlation
#'   matrices that were used for the fit propensity analysis. This slot is only populated
#'   if \code{SaveR} is set to \code{TRUE}.
#' @slot mod_list A list of the same length as the number of models being compared. Each element
#'   contains a list of the same length as the number of replications and contains fitted lavaan
#'   models. Only populated if \code{saveModel} is set to \code{TRUE}.
#' @references
#' Bonifay, W. E., & Cai, L. (2017). On the complexity of item response theory models. Multivariate Behavioral Research, 52(4), 465–484. \url{http://doi.org/10.1080/00273171.2017.1309262}
#'
#' Falk, C. F., & Muthukrishna, M. (in press). Parsimony in model selection: Tools for assessing fit propensity. Psychological Methods.
#'
#' Lewandowski, D., Kurowicka, D., & Joe, H. (2009). Generating random correlation matrices based on vines and extended onion method. Journal of Multivariate Analysis,100(9), 1989–2001. \url{http://doi.org/10.1016/j.jmva.2009.04.008}
#'
#' Joe, H. (2006). Generating random correlation matrices based on partial correlations. Journal of Multivariate Analysis, 97(10), 2177–2189. \url{http://doi.org/10.1016/j.jmva.2005.05.010}
#'
#' Preacher, K. J. (2003). The role of model complexity in the evaluation of structural equation models (PhD thesis). The Ohio State University.
#'
#' Preacher, K. J. (2006). Quantifying parsimony in structural equation modeling. Multivariate Behavioral Research, 41(3), 227–259. \url{http://doi.org/10.1207/s15327906mbr4103_1}
#' @examples
#' \donttest{
#' # Set up a covariance matrix to fit models to
#' p<-3 # number of variables
#' temp_mat <- diag(p) # identity matrix
#' colnames(temp_mat) <- rownames(temp_mat) <- paste0("V", seq(1, p))
#'
#' # Define and fit two models using lavaan package
#' mod1a <- 'V3 ~ V1 + V2
#'   V1 ~~ 0*V2'
#' mod2a <- 'V3 ~ V1
#'   V2 ~ V3'
#'
#' mod1a.fit <- sem(mod1a, sample.cov=temp_mat, sample.nobs=500)
#' mod2a.fit <- sem(mod2a, sample.cov=temp_mat, sample.nobs=500)
#'
#' # Run fit propensity analysis a variety of different ways
#'
#' # Onion approach, only positive correlation matrices, save srmr
#' res <- run.fitprop(mod1a.fit, mod2a.fit, fit.measure="srmr",
#'   rmethod="onion",reps=1000,onlypos=TRUE)
#' summary(res)
#'
#' # Onion approach, save several fit indices
#' res <- run.fitprop(mod1a.fit, mod2a.fit, fit.measure=c("srmr","cfi","rmsea"),
#'   rmethod="onion",reps=1000)
#' summary(res)
#'
#' # mcmc approach, with parallel processing (4 cores)
#' # Save and then access correlation matrices and fitted models
#' # Note: this will take a very long time as the default number
#' # of iterations is set to a lot
#' cl<-makeCluster(2)
#' res <- run.fitprop(mod1a.fit, mod2a.fit, fit.measure="srmr",
#'   rmethod="mcmc",reps=1000,cluster=cl,saveModel=TRUE,saveR=TRUE)
#' stopCluster(cl)
#' summary(res)
#' res$R[[1]] # Correlation matrix for first replication
#' res$fit_list[[1]] # Saved fit indices for first model
#' res$mod_list[[1]][[2]] # Fitted lavaan model for first model, second replication
#'
#' # mcmc approach, overriding defaults
#' ctrl<-list(
#'   iter = 10000, # total number of iterations
#'   miniter = 5000, # but, min number of iterations per chain
#'   jmpsize = .3
#' )
#' cl<-makeCluster(2)
#' res <- run.fitprop(mod1a.fit, mod2a.fit, fit.measure="srmr",
#'   rmethod="mcmc",reps=1000,cluster=cl, mcmc.args=ctrl)
#' stopCluster(cl)
#' summary(res)
#'
#' }
#' @return An object of class fitprop for which plot and summary methods are available. Some slots are listed in the Slots section.
#' @seealso \code{\link[ockhamSEM]{plot.fitprop}} \code{\link[ockhamSEM]{summary.fitprop}}
#' @export
#' @importFrom lavaan lavInspect lavNames parTable lavaan
#' @import parallel
run.fitprop <- function(...,
                        fit.measure="srmr",
                        rmethod=c("onion","mcmc","clustergen"),
                        reps=1000,
                        onlypos=FALSE,
                        seed=1234,
                        mcmc.args = list(),
                        clustergen.args = list(),
                        saveModel=FALSE,
                        saveR=FALSE,
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

  #Set matrix for storing outputs
  fit_list = vector("list", n.models)
  for (i in 1:length(fit_list)) {
    fit_list[[i]] = matrix(0,reps,length(fit.measure))
    colnames(fit_list[[i]]) <- fit.measure
  }

  if(rmethod=="onion") {
    control <- onion.args.default(d,clustergen.args)
  } else if (rmethod=="mcmc") {
    if(!is.null(cluster)){
      control <- mcmc.args.default(length(cluster),d,mcmc.args)
    } else {
      control <- mcmc.args.default(1,d,mcmc.args)
    }
  } else if (rmethod == "clustergen") {
    control <- clustergen.args.default(d,clustergen.args)
  }

  #Generate the matrices
  print("Generate matrices")
  if(!is.null(cluster)){
    clusterSetRNGStream(cluster, seed)
    nmat<-clusterSplit(cluster, 1:reps)
    out_mat <- parLapply(cluster, nmat, genmat, rmethod=rmethod, control=control, onlypos=onlypos)
  } else {
    set.seed(seed)
    out_mat <- lapply(list(1:reps), genmat, rmethod=rmethod, control=control, onlypos=onlypos)
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
                            sample.nobs= nobs,
                            slotOptions = opt)

    # update may not work with older lavaan versions
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
      fit_tmp<-parLapply(cluster,1:reps,fitmod,lavmodel=lavmodel,out_mat=out_mat,vnames=vnames,j=j,saveModel=saveModel,d=d,fit.measure=fit.measure)
    } else {
      fit_tmp<-lapply(1:reps,fitmod,lavmodel=lavmodel,out_mat=out_mat,vnames=vnames,j=j,saveModel=saveModel,d=d,fit.measure=fit.measure)
    }
    fit_list[[j]]<-matrix(unlist(fit_tmp),ncol=length(fit.measure),byrow=TRUE)
    colnames(fit_list[[j]])<-fit.measure

    if(saveModel){
      mod_list[[j]]<-lapply(fit_tmp,attr,"mod")
    }
  }

  # Process number of missing values (if any) per model and fit index
  na_list<-list()
  complete_mat<-matrix(1,nrow=reps,ncol=n.fit)
  for(j in 1:length(fit_list)){
    na_list[[j]]<-is.na(fit_list[[j]])
    complete_mat[which(is.na(fit_list[[j]]))]<-NA
  }

  out<-list()
  out$fun.call<-fun.call
  out$reps=reps
  out$rmethod=rmethod
  out$onlypos=onlypos
  out$nfit=n.fit
  out$origmodels=models

  out$na_list<-na_list
  out$complete_mat<-complete_mat

  out$fit_list<-fit_list

  if(saveR){
    out$R<-list()
    for(j in 1:reps){
      out$R[[j]]<-matrix(out_mat[,j],d,d)
    }
  }

  if(saveModel){
    out$mod_list<-mod_list
  }

  class(out)<-"fitprop"
  return(out)
}

#' Plot function for fitprop objects
#' @param x Object of class fitprop, as created by run.fitprop function.
#' @param ... Does nothing but to hopefully make this generic function pass R CMD check
#' @param type What type of plot to produce? Options include \code{"ecdf"} or \code{"euler"}. \code{"nVennR"}
#'   is not yet operational.
#' @param whichmod Index number corresponding to which model(s) to include on the plot. Defaults to all models.
#' @param whichfit Character vector indicating which indices of model fit to include. Defaults to all saved indices.
#' @param savePlot Logical value indicating whether to save plot to a list (TRUE) or just produce plot to output.
#' @param xlim Numeric vector of length 2 indicating the limits for the x-axis of the plot.
#' @param samereps Logical value indicating whether to use only results from replications in which all selected models yielded results.
#' @param cutoff Numeric vector indicating what cut value of the fit indice(s) to use for euler plots.
#' @param lower.tail Logical vector indicating whether lower values of each fit index corresponds to good fit.
#' @param mod.lab Optional character vector of labels for each model.
#' @param mod.brewer.pal Optional string corresponding to the palette from RColorBrewer to use for the different models. e.g.,
#'   see \code{\link[RColorBrewer]{display.brewer.all}}.
#' @examples
#' \donttest{
#'
#' library(ggplot2)
#' # Set up a covariance matrix to fit models to
#' p<-3 # number of variables
#' temp_mat <- diag(p) # identity matrix
#' colnames(temp_mat) <- rownames(temp_mat) <- paste0("V", seq(1, p))
#'
#' # Define and fit two models using lavaan package
#' mod1a <- 'V3 ~ V1 + V2
#'   V1 ~~ 0*V2'
#' mod2a <- 'V3 ~ V1
#'   V2 ~ V3'
#'
#' mod1a.fit <- sem(mod1a, sample.cov=temp_mat, sample.nobs=500)
#' mod2a.fit <- sem(mod2a, sample.cov=temp_mat, sample.nobs=500)
#'
#' # Run fit propensity analysis
#' # Onion approach, save srmr and CFI
#' res <- run.fitprop(mod1a.fit, mod2a.fit, fit.measure=c("srmr","cfi"),
#'   rmethod="onion",reps=1000)
#' summary(res)
#'
#' # Generate a variety of plots
#' # ecdf - only srmr
#' plot(res, whichfit = "srmr")
#'
#' # ecdf - both srmr and cfi, properly indicating that higher values of cfi are better
#' plot(res, lower.tail=c(TRUE,FALSE))
#'
#' # ecdf change color palette, save plot to object, then display or change
#' myplot<-plot(res, whichfit="cfi",
#'              lower.tail=FALSE,
#'              mod.brewer.pal="Dark2",
#'              savePlot=TRUE)
#' myplot[[1]] # saved to the first slot
#' myplot[[1]] + theme_bw() # change something, like the theme
#'
#' # euler plot
#' # cfi - which models have cfi of .5 or better?
#' plot(res, type="euler", whichfit="cfi", lower.tail=FALSE, cutoff=.5)
#'
#' }
#' @export
#' @import ggplot2
#' @importFrom graphics plot
#' @importFrom stats na.omit
#' @importFrom tidyr gather
#' @importFrom eulerr euler
#' @importFrom nVennR plotVenn
#' @importFrom rlang .data
plot.fitprop<-function(x,...,type=c("ecdf","euler","nVennR"),whichmod=NULL,whichfit=colnames(x$fit_list[[1]]),savePlot=FALSE,
                       xlim=c(0,1),samereps=TRUE,cutoff=rep(.1,length(whichfit)),lower.tail=rep(TRUE,length(whichfit)),
                       mod.lab=NULL,mod.brewer.pal="Set1"){

  type<-match.arg(type)

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
      if(samereps){
        dat<-na.omit(dat)
      }
      dat<-gather(dat,"variable","value",whichmod)
      graph<-ggplot(dat,aes(x=.data$value))+
        stat_ecdf(aes(linetype=.data$variable, color=.data$variable),na.rm=TRUE,size=.7) +
        scale_color_brewer(palette=mod.brewer.pal)

      graph <- graph + guides(color=guide_legend(title="Model"))
      graph <- graph + guides(linetype=guide_legend(title="Model"))

      if(lower.tail[m]){
        graph<-graph+ xlim(xlim[1],xlim[2])
      } else {
        graph<-graph+ xlim(xlim[2],xlim[1])
      }
      graph<-graph+xlab(fm) + ylab("Cumulative Probability") #theme(legend.title=element_blank())+

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

      stop("nVennR is not yet functional")

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
"  Only positive R?        ", x$onlypos, "\n"
)

}

#' Summary function for fitprop objects
#' @param object Object of class fitprop, as created by \code{\link[ockhamSEM]{run.fitprop}} function.
#' @param ... Does nothing but to hopefully make this generic function pass R CMD check.
#' @param probs Vector passed to quantile to determine what probabilities to report.
#' @param samereps Logical value indicating whether to use only results from replications in which all selected models yielded results.
#' @param lower.tail Logical vector indicating whether lower values of each fit index corresponds to good fit.
#' @param NML (experimental) Logical value indicating whether to compute normalized maximum likelihood (NML; e.g., Rissanen, 2001). Requires
#'   that `logl` is a saved fit index.
#' @param UIF (experimental) Logical value indicating whether to compute uniform index of fit (UIF; Botha, Shapiro, Steiger, 1988).
#'   Original paper appeared to use least-squares estimation and compute UIF based on proportion of times that obtained fit function
#'   was better than fit function based on random data. Currently this option works for any fit index, but some may make sense more
#'   than others. Also requires \code{lower.tail} is appropriately set.
#' @references
#' Botha, J.D., Shapiro, A., \& Steiger, J.H. (1988). Uniform indices-of-fit for factor analysis models. Multivariate Behavioral Research, 23(4), 443-450. \url{http://doi.org/10.1207/s15327906mbr2304_2}
#'
#' Rissanen, J. (2001). Strong optimality of the normalized ML models as universal codes and information in data. IEEE Transactions on Information Theory, 47, 1712–1717.
#'
#' @examples
#' \donttest{
#'
#' # Borrow PoliticalDemocracy data
#' data(PoliticalDemocracy)
#'
#' # Define and fit two models using lavaan package
#' mod1a <- 'y5 ~ y1 + x1
#'   y1 ~~ 0*x1'
#' mod2a <- 'y5 ~ y1
#'   x1 ~ y5'
#'
#' mod1a.fit <- sem(mod1a, sample.cov=cov(PoliticalDemocracy), sample.nobs=500)
#' mod2a.fit <- sem(mod2a, sample.cov=cov(PoliticalDemocracy), sample.nobs=500)
#'
#' # Run fit propensity analysis
#' # Onion approach, save srmr and CFI
#' res <- run.fitprop(mod1a.fit, mod2a.fit, fit.measure=c("srmr","cfi"),
#'   rmethod="onion",reps=1000)
#'
#' # Generate summaries
#' summary(res)
#'
#' # sort quantiles differently for srmr and cfi
#' # Use different quantiles
#' summary(res, probs = c(0,.25,.5,.75,1), lower.tail=c(TRUE,FALSE))
#'
#' # If some models failed to converge, this would result in
#' # summaries computed on possibly different replications:
#' # Use different quantiles
#' summary(res, samereps=FALSE, lower.tail=c(TRUE,FALSE))
#'
#' # For computing NML (experimental)
#' res <- run.fitprop(mod1a.fit, mod2a.fit, fit.measure="logl",
#'   rmethod="onion",reps=2500)
#'
#' summary(res, NML=TRUE, lower.tail=FALSE)
#'
#' # For computing UIF (experimental)
#' # Orig UIF used least-squares estimation and examined fit function
#'
#' mod1a.fit <- sem(mod1a, sample.cov=cov(PoliticalDemocracy), sample.nobs=500,
#'   estimator="ULS")
#' mod2a.fit <- sem(mod2a, sample.cov=cov(PoliticalDemocracy), sample.nobs=500,
#'   estimator="ULS")
#'
#' res <- run.fitprop(mod1a.fit, mod2a.fit, fit.measure="fmin",
#'   rmethod="onion",reps=2500)
#'
#' summary(res, UIF=TRUE, lower.tail=TRUE)
#'
#'
#'
#' }
#' @export
#' @importFrom stats quantile median ks.test na.omit
#' @importFrom utils str
#' @importFrom effsize cliff.delta cohen.d
#' @importFrom matrixStats logSumExp
summary.fitprop<-function(object,...,probs=seq(0,1,.1),samereps=TRUE,lower.tail=rep(TRUE,ncol(object$fit_list[[1]])),
                          NML = FALSE, UIF=FALSE){

  data<-object$fit_list
  nmod<-length(data) # number of models
  nfit<-ncol(data[[1]]) # number of fit measures
  nrep<-object$reps # number of available replications

  stats<-list()
  quantiles<-list()
  efs<-list()
  if(UIF){uif<-list()}

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

    fitnames<-colnames(data[[mod]])

    finite<-apply(data[[mod]],2,function(x){sum(is.finite(x))})
    na<-apply(data[[mod]],2,function(x){sum(is.na(x))})
    nrep<-nrow(data[[mod]])
    if(samereps){
      means<-NULL
      medians<-NULL
      uifs<-NULL

      for(j in 1:nfit){

        datj<-data[[mod]][!is.na(object$complete_mat[,j]),j]

        means<-c(means,mean(datj))
        medians<-c(medians,median(datj))

        if(UIF){
          obt<-as.numeric(fitMeasures(object$origmodels[[mod]],fitnames[j]))
          n.uif <- length(datj)
          if(lower.tail[j]){
            uif.mod <- sum(obt<=datj)/n.uif
          } else {
            uif.mod <- sum(obt>=datj)/n.uif
          }
          uifs<-c(uifs,uif.mod)
        }
      }
    } else {
      means<-colMeans(data[[mod]],na.rm=TRUE)
      medians<-apply(data[[mod]],2,median,na.rm=TRUE)

      if(UIF){
        for(j in 1:nfit){
          obt<-as.numeric(fitMeasures(object$origmodels[[mod]],fitnames[j]))
          n.uif <- length(data[[mod]][!is.na(object$complete_mat[,j]),j])
          if(lower.tail[j]){
            uif.mod <- sum(obt<=data[[mod]][!is.na(object$complete_mat[,j]),j])/n.uif
          } else {
            uif.mod <- sum(obt>=data[[mod]][!is.na(object$complete_mat[,j]),j])/n.uif
          }
          uifs<-c(uifs,uif.mod)
        }
      }
    }

    names(means)<-names(medians)<-fitnames
    if(UIF){names(uifs)<-fitnames}

    stats[[mod]]<-list()
    stats[[mod]]$finite<-finite
    stats[[mod]]$finite<-na
    stats[[mod]]$means<-means
    stats[[mod]]$medians<-medians
    if(UIF){stats[[mod]]$uifs<-uifs}

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

          efs[[j]][[indx]]$d<-cohen.d(tmp,as.factor(grp.tmp))$estimate
          efs[[j]][[indx]]$delta<-cliff.delta(data[[mod]][!is.na(object$complete_mat[,j]),j],data[[mod2]][!is.na(object$complete_mat[,j]),j])$estimate
          efs[[j]][[indx]]$ks<-ks.test(data[[mod]][!is.na(object$complete_mat[,j]),j],data[[mod2]][!is.na(object$complete_mat[,j]),j])$statistic

        } else {
          tmp<-c(data[[mod]][,j],data[[mod2]][,j])
          tmp<-na.omit(tmp)
          grp.tmp<-c(rep(1,nrep),rep(2,nrep))
          if(!is.null(attr(tmp,"na.action"))){
            grp.tmp<-grp.tmp[-attr(tmp,"na.action")]
          }

          efs[[j]][[indx]]$d<-cohen.d(tmp,as.factor(grp.tmp))$estimate
          efs[[j]][[indx]]$delta<-cliff.delta(data[[mod]][,j],data[[mod2]][,j])$estimate
          efs[[j]][[indx]]$ks<-ks.test(data[[mod]][,j],data[[mod2]][,j])$statistic
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

  if(NML | UIF){
    cat(
      "\n",
      "Fit indices for obtained model\n"
    )

    if(NML & "logl" %in% colnames(data[[1]])){
      cat(
        "\n",
        "-2*Log-Normalized Maximum Likelihood:\n"
      )

      #NML (Rissanen)
      nmls<-vector("numeric")
      for(mod in 1:nmod){

        # remove any clearly invalid logl values (these can't be >0, actually)
        ll.tmp<-as.numeric(data[[mod]][,"logl"])
        ll.tmp<-ll.tmp[ll.tmp<0 & !is.na(ll.tmp)]
        n.ll <- length(ll.tmp)

        # ll obtained
        ll.obt<-as.numeric(fitMeasures(object$origmodels[[mod]],"logl"))

        num <- ll.obt
        denom <- logSumExp(ll.tmp) - log(n.ll)

        # what about log-nml?
        log.nml.mod <- -2*((ll.obt) - (logSumExp(ll.tmp) - log(n.ll)))
        cat("\n", paste0("Model ",mod, ": "), round(log.nml.mod, 5))
        nmls<-c(nmls,log.nml.mod)
      }
      out[["nml"]]<-nmls
    }


    if(UIF){
      cat("\nUniform Index of Fit (Botha, Shapiro, Steiger, 1988)\n")
      for(mod in 1:nmod){
        cat("\n Model ",mod,"\n")
        print.default(round(stats[[mod]]$uifs,3))
      }
    }

  }

  invisible(out)
}

# Internal function used for fitting models

#' @importFrom lavaan lavaan lavInspect fitMeasures
fitmod <- function(indx,lavmodel,out_mat,vnames,j,saveModel,d,fit.measure){

  temp_matrix <- matrix(out_mat[,indx],d,d)
  colnames(temp_matrix) <- rownames(temp_matrix) <- vnames

  # added from ShortForm Tabu code; nuke all starting values? Does this work? Not with older lavaan version
  lavmodel@ParTable$est<-NULL
  lavmodel@ParTable$start<-NULL

  # replaced with update()?; might be slower, but more stable? (doesn't work with older lavaan version)
  #res<-try({my_fitted_model<-update(lavmodel, sample.cov=temp_matrix)

  res<-try ({
    lavmodel@Options$se="none"
    lavmodel@Options$start="default"
    lavmodel@Options$do.fit=TRUE
    #lavmodel@Data
    #lavmodel@SampleStats
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

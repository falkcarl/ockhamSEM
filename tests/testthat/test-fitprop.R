test_that("example 1", {

  ## Example 1 from the paper
  mod1a <- 'V3 ~ V1 + V2
    V1 ~~ 0*V2'
  mod2a <- 'V3 ~ V1
    V2 ~ V3'
  p<-3 # number of variables
  temp_mat <- diag(p) # identity matrix

  # set row and column names
  colnames(temp_mat) <- rownames(temp_mat) <- paste0("V", seq(1, p))

  mod1a.fit <- sem(mod1a, sample.cov=temp_mat, sample.nobs=500)
  mod2a.fit <- sem(mod2a, sample.cov=temp_mat, sample.nobs=500)

  res.on <- run.fitprop(mod1a.fit, mod2a.fit, fit.measure="srmr",
                        rmethod="onion",reps=5000,onlypos=TRUE)

  cl <- makeCluster(2)
  res.mcmc <- run.fitprop(mod1a.fit,mod2a.fit,fit.measure="srmr",
                          rmethod = "mcmc", reps = 5000, onlypos=TRUE,
                          cluster=cl)
  stopCluster(cl)

  #plot1<-plot(res.on, savePlot=TRUE,
  #            mod.lab=c("Model 1A","Model 2A"),
  #            mod.brewer.pal="Set1")
  #plot2<-plot(res.mcmc, savePlot=TRUE,
  #            mod.lab=c("Model 1A","Model 2A"),
  #            mod.brewer.pal="Set1")

  expect_snapshot(summary(res.mcmc, probs=c(.25,.5,.75)))

})

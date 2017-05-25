params.wd = "."
getQualityMeasuresScript = file("getQualityMeasures.R")

SKAT = 0
CHISQ = 1

CONSISTENCY = 1
AICc = 2
BIC = 3
AICcn = 4

rgwas = file("$params.rgwas")
rcausal = file("$params.rcausal")

process runGin_SKAT_AICc {

  publishDir "$params.wd", overwrite: true

  input:
    file rgwas
    file rcausal
    file getQualityMeasuresScript

  output:
    file "gin.skat.aic.*.RData" into gin_skat_aic_rdata

  """
  #!/usr/bin/env Rscript
  library(martini)
  source("$getQualityMeasuresScript")
  load("$rgwas")
  load("$rcausal")

  test <- "gin_SKAT_AICc"
  fit <- runGin(gwas\$X, gwas\$Y, gwas\$net, list(test_statistic = $SKAT,
                                                     gridsearch_depth = 3,
                                                     selection_criterion = $AICc))
  qual <- getQualityMeasures(fit\$indicator, as.numeric(causal), test)

  save(test, fit, qual, file = paste0("gin.skat.aic.", id, ".RData"))
  """

}

process runGin_CHISQ_AICc {

  publishDir "$params.wd", overwrite: true

  input:
    file rgwas
    file rcausal
    file getQualityMeasuresScript

  output:
    file "gin.chisq.aic.*.RData" into gin_chisq_aic_rdata

  """
  #!/usr/bin/env Rscript
  library(martini)
  source("$getQualityMeasuresScript")
  load("$rgwas")
  load("$rcausal")

  test <- "gin_CHISQ_AICc"
  fit <- runGin(gwas\$X, gwas\$Y, gwas\$net, list(test_statistic = $CHISQ,
                                                     gridsearch_depth = 3,
                                                     selection_criterion = $AICc))
  qual <- getQualityMeasures(fit\$indicator, as.numeric(causal), test)

  save(test, fit, qual, file = paste0("gin.chisq.aic.", id, ".RData"))
  """

}

process runGin_SKAT_CONSISTENCY {

  publishDir "$params.wd", overwrite: true

  input:
    file rgwas
    file rcausal
    file getQualityMeasuresScript

  output:
    file "gin.skat.consistency.*.RData" into gin_skat_consistency_rdata

  """
  #!/usr/bin/env Rscript
  library(martini)
  source("$getQualityMeasuresScript")
  load("$rgwas")
  load("$rcausal")

  test <- "gin_SKAT_CONSISTENCY"
  fit <- runGin(gwas\$X, gwas\$Y, gwas\$net, list(test_statistic = $SKAT,
                                                     gridsearch_depth = 3,
                                                     selection_criterion = $CONSISTENCY))
  qual <- getQualityMeasures(fit\$indicator, as.numeric(causal), test)

  save(test, fit, qual, file = paste0("gin.skat.consistency.", id, ".RData"))
  """

}

process runGin_CHISQ_CONSISTENCY {

  publishDir "$params.wd", overwrite: true

  input:
    file rgwas
    file rcausal
    file getQualityMeasuresScript

  output:
    file "gin.chisq.consistency.*.RData" into gin_chisq_consistency_rdata

  """
  #!/usr/bin/env Rscript
  library(martini)
  source("$getQualityMeasuresScript")
  load("$rgwas")
  load("$rcausal")

  test <- "gin_CHISQ_CONSISTENCY"
  fit <- runGin(gwas\$X, gwas\$Y, gwas\$net, list(test_statistic = $CHISQ,
                                                     gridsearch_depth = 3,
                                                     selection_criterion = $CONSISTENCY))
  qual <- getQualityMeasures(fit\$indicator, as.numeric(causal), test)

  save(test, fit, qual, file = paste0("gin.chisq.consistency.", id, ".RData"))
  """

}


process runLasso {

  publishDir "$params.wd", overwrite: true

  input:
    file rgwas
    file rcausal
    file getQualityMeasuresScript

  output:
    file "lasso.*.RData" into lasso_rdata

  """
  #!/usr/bin/env Rscript
  library(glmnet)
  source("$getQualityMeasuresScript")
  load("$rgwas")
  load("$rcausal")

  test <- "LASSO"
  fit.cv <- cv.glmnet(gwas\$X, gwas\$Y, family = "binomial", type.measure = "auc")
  fit <- glmnet(gwas\$X, gwas\$Y, lambda = fit.cv\$lambda.1se)

  qual <- getQualityMeasures(as.numeric(fit\$beta != 0), as.numeric(causal), test)

  save(test, fit.cv, fit, qual, file = paste0("lasso.", id, ".RData"))
  """

}

params.wd = "."
getQualityMeasuresScript = file("getQualityMeasures.R")

SKAT = 0
CHISQ = 1

CONSISTENCY = 1
AICc = 2
BIC = 3
AICcn = 4

rdata = file("$params.rdata")

process runScones_SKAT_AICc {

  publishDir "$params.wd", overwrite: true

  input:
    file rdata
    file getQualityMeasuresScript

  output:
    file "scones.skat.aic.*.RData" into scones_skat_aic_rdata

  """
  #!/usr/bin/env Rscript
  library(rscones2)
  source("$getQualityMeasuresScript")
  load("$rdata")

  test <- "scones SKAT AICc"
  fit <- runScones(gwas\$X, gwas\$Y, gwas\$net, list(test_statistic = $SKAT,
                                                     gridsearch_depth = 3,
                                                     selection_criterion = $AICc))
  qual <- getQualityMeasures(fit\$indicator, truth, test)

  save(test, fit, qual, file = paste0("scones.skat.aic.", id, ".RData"))
  """

}

process runScones_CHISQ_AICc {

  publishDir "$params.wd", overwrite: true

  input:
    file rdata
    file getQualityMeasuresScript

  output:
    file "scones.chisq.aic.*.RData" into scones_chisq_aic_rdata

  """
  #!/usr/bin/env Rscript
  library(rscones2)
  source("$getQualityMeasuresScript")
  load("$rdata")

  test <- "scones CHISQ AICc"
  fit <- runScones(gwas\$X, gwas\$Y, gwas\$net, list(test_statistic = $CHISQ,
                                                     gridsearch_depth = 3,
                                                     selection_criterion = $AICc))
  qual <- getQualityMeasures(fit\$indicator, truth, test)

  save(test, fit, qual, file = paste0("scones.chisq.aic.", id, ".RData"))
  """

}

process runScones_SKAT_CONSISTENCY {

  publishDir "$params.wd", overwrite: true

  input:
    file rdata
    file getQualityMeasuresScript

  output:
    file "scones.skat.consistency.*.RData" into scones_skat_consistency_rdata

  """
  #!/usr/bin/env Rscript
  library(rscones2)
  source("$getQualityMeasuresScript")
  load("$rdata")

  test <- "scones SKAT CONSISTENCY"
  fit <- runScones(gwas\$X, gwas\$Y, gwas\$net, list(test_statistic = $SKAT,
                                                     gridsearch_depth = 3,
                                                     selection_criterion = $CONSISTENCY))
  qual <- getQualityMeasures(fit\$indicator, truth, test)

  save(test, fit, qual, file = paste0("scones.skat.consistency.", id, ".RData"))
  """

}

process runScones_CHISQ_CONSISTENCY {

  publishDir "$params.wd", overwrite: true

  input:
    file rdata
    file getQualityMeasuresScript

  output:
    file "scones.chisq.consistency.*.RData" into scones_chisq_consistency_rdata

  """
  #!/usr/bin/env Rscript
  library(rscones2)
  source("$getQualityMeasuresScript")
  load("$rdata")

  test <- "scones CHISQ CONSISTENCY"
  fit <- runScones(gwas\$X, gwas\$Y, gwas\$net, list(test_statistic = $CHISQ,
                                                     gridsearch_depth = 3,
                                                     selection_criterion = $CONSISTENCY))
  qual <- getQualityMeasures(fit\$indicator, truth, test)

  save(test, fit, qual, file = paste0("scones.chisq.consistency.", id, ".RData"))
  """

}


process runLasso {

  publishDir "$params.wd", overwrite: true

  input:
    file rdata
    file getQualityMeasuresScript

  output:
    file "lasso.*.RData" into lasso_rdata

  """
  #!/usr/bin/env Rscript
  library(glmnet)
  source("$getQualityMeasuresScript")
  load("$rdata")

  test <- "LASSO"
  fit.cv <- cv.glmnet(gwas\$X, gwas\$Y, family = "binomial", type.measure = "auc")
  fit <- glmnet(gwas\$X, gwas\$Y, lambda = fit.cv\$lambda.1se)

  qual <- getQualityMeasures(as.numeric(fit\$beta != 0), truth, test)

  save(test, fit.cv, fit, qual, file = paste0("lasso.", id, ".RData"))
  """

}

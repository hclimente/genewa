params.wd = "."

getQualityMeasuresScript = file("getQualityMeasures.R")

associationScores = ["chi2"]
modelScores = ["bic", "aic", "consistency"]

rgwas = file("$params.rgwas")
rcausal = file("$params.rcausal")

process run_evo {

  publishDir "$params.wd", overwrite: true

  input:
  each c from associationScores
  each m from modelScores
  file rgwas
  file rcausal
  file getQualityMeasuresScript

  output:
  file "evo.$c.$m.*.RData" into evo_rdata

  """
  #!/usr/bin/env Rscript
  library(martini)
  source("$getQualityMeasuresScript")
  load("$rgwas")
  load("$rcausal")

  test <- "evo.$c.$m"

  start.time <- Sys.time()
  cones <- search_cones(gwas, net, associationScore = "$c", modelScore = "$m")
  end.time <- Sys.time()
  time.taken <- end.time - start.time

  qual <- getQualityMeasures(as.numeric(cones\$selected), as.numeric(causal), test)

  save(test, cones, qual, time.taken, file = paste(test, id, "RData", sep = "."))
  """

}

process run_scones {

  publishDir "$params.wd", overwrite: true

  input:
    file rgwas
    file rcausal
    file getQualityMeasuresScript

  output:
    file "scones.*.RData" into evo_chisq_consistency_rdata

  """
  #!/usr/bin/env Rscript
  library(martini)
  source("$getQualityMeasuresScript")
  load("$rgwas")
  load("$rcausal")

  test <- "scones"

  start.time <- Sys.time()
  cones <- find_cones(gwas, net)
  end.time <- Sys.time()
  time.taken <- end.time - start.time

  qual <- getQualityMeasures(as.numeric(cones\$selected), as.numeric(causal), test)

  save(test, cones, qual, time.taken, file = paste(test, id, "RData", sep = "."))
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
  library(snpStats)
  source("$getQualityMeasuresScript")
  load("$rgwas")
  load("$rcausal")

  test <- "lasso"

  X <- as(gwas\$genotypes, "numeric")
  Y <- gwas\$fam\$affected

  start.time <- Sys.time()
  fit.cv <- cv.glmnet(X, Y, family = "binomial", type.measure = "auc")
  fit <- glmnet(X, Y, lambda = fit.cv\$lambda.1se)
  end.time <- Sys.time()
  time.taken <- end.time - start.time

  qual <- getQualityMeasures(as.numeric(fit\$beta != 0), as.numeric(causal), test)
  cones <- gwas\$map
  cones\$selected <- fit\$beta != 0

  save(test, fit.cv, fit, cones, qual, time.taken, file = paste(test, id, "RData", sep = "."))
  """

}

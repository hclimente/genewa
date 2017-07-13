params.wd = "."
params.scones = true
params.shake = false
params.lasso = false

getQualityMeasuresScript = file("getQualityMeasures.R")

SKAT = 0
CHISQ = 1

CONSISTENCY = 1
AICc = 2
BIC = 3
AICcn = 4

rgwas = file("$params.rgwas")
rcausal = file("$params.rcausal")

if (params.hasProperty('shake') && params.shake) {

  process run_shake_SKAT_BIC {

    publishDir "$params.wd", overwrite: true

    input:
      file rgwas
      file rcausal
      file getQualityMeasuresScript

    output:
      file "shake.skat.bic.*.RData" into shake_skat_bic_rdata

    """
    #!/usr/bin/env Rscript
    library(martini)
    source("$getQualityMeasuresScript")
    load("$rgwas")
    load("$rcausal")

    test <- "shake_SKAT_BIC"
    map <- shake(gwas, net, test_statistic = $SKAT, selection_criterion = $BIC, gridsearch_depth = 3)
    qual <- getQualityMeasures(as.numeric(map\$selected), as.numeric(causal), test)

    save(test, map, qual, file = paste(test, id, "RData", sep = "."))
    """

  }

  process run_shake_CHISQ_BIC {

    publishDir "$params.wd", overwrite: true

    input:
      file rgwas
      file rcausal
      file getQualityMeasuresScript

    output:
      file "shake.chisq.bic.*.RData" into shake_chisq_bic_rdata

    """
    #!/usr/bin/env Rscript
    library(martini)
    source("$getQualityMeasuresScript")
    load("$rgwas")
    load("$rcausal")

    test <- "shake_CHISQ_BIC"
    map <- shake(gwas, net, test_statistic = $CHISQ, selection_criterion = $BIC, gridsearch_depth = 3)
    qual <- getQualityMeasures(as.numeric(map\$selected), as.numeric(causal), test)

    save(test, map, qual, file = paste(test, id, "RData", sep = "."))
    """

  }

  process run_shake_SKAT_CONSISTENCY {

    publishDir "$params.wd", overwrite: true

    input:
      file rgwas
      file rcausal
      file getQualityMeasuresScript

    output:
      file "shake.skat.consistency.*.RData" into shake_skat_consistency_rdata

    """
    #!/usr/bin/env Rscript
    library(martini)
    source("$getQualityMeasuresScript")
    load("$rgwas")
    load("$rcausal")

    test <- "shake_SKAT_CONSISTENCY"
    map <- shake(gwas, net, test_statistic = $SKAT, selection_criterion = $CONSISTENCY, gridsearch_depth = 3)
    qual <- getQualityMeasures(as.numeric(map\$selected), as.numeric(causal), test)

    save(test, map, qual, file = paste(test, id, "RData", sep = "."))
    """

  }

  process run_shake_CHISQ_CONSISTENCY {

    publishDir "$params.wd", overwrite: true

    input:
      file rgwas
      file rcausal
      file getQualityMeasuresScript

    output:
      file "shake.chisq.consistency.*.RData" into shake_chisq_consistency_rdata

    """
    #!/usr/bin/env Rscript
    library(martini)
    source("$getQualityMeasuresScript")
    load("$rgwas")
    load("$rcausal")

    test <- "shake_CHISQ_CONSISTENCY"
    map <- shake(gwas, net, test_statistic = $CHISQ, selection_criterion = $CONSISTENCY, gridsearch_depth = 3)
    qual <- getQualityMeasures(as.numeric(map\$selected), as.numeric(causal), test)

    save(test, map, qual, file = paste(test, id, "RData", sep = "."))
    """

  }
}

process run_scones {

  publishDir "$params.wd", overwrite: true

  input:
    file rgwas
    file rcausal
    file getQualityMeasuresScript

  output:
    file "scones.*.RData" into shake_chisq_consistency_rdata

  """
  #!/usr/bin/env Rscript
  library(martini)
  source("$getQualityMeasuresScript")
  load("$rgwas")
  load("$rcausal")

  test <- "scones"
  map <- shake(gwas, net)
  qual <- getQualityMeasures(as.numeric(map\$selected), as.numeric(causal), test)

  save(test, map, qual, file = paste(test, id, "RData", sep = "."))
  """
}

if (params.hasProperty('others') && params.others) {

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

    test <- "LASSO"

    X <- as(gwas\$genotypes, "numeric")
    Y <- gwas\$fam\$affected
    fit.cv <- cv.glmnet(X, Y, family = "binomial", type.measure = "auc")
    fit <- glmnet(X, Y, lambda = fit.cv\$lambda.1se)

    qual <- getQualityMeasures(as.numeric(fit\$beta != 0), as.numeric(causal), test)

    save(test, fit.cv, fit, qual, file = paste(test, id, "RData", sep = "."))
    """

  }
}

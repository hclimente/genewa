#!/usr/bin/env nextflow

// Help message
helpMessage = """
Launch several analysis on a GWAS.

Usage:

  nextflow martiniflow/analyze/analyze_gwas.nf --rgwas gwas.RData --rnet net.RData

(rgwas, rnet) -> (*.RData)

INPUT

- rgwas                snpMatrix with a GWAS experiment and an info list.
- rnet                 An igraph network and a netType.

OUTPUT

- *.RData              RData files containing at least test, info, cones and time.taken.

PARAMETERS

- Required
  --rgwas              snpMatrix with a GWAS experiment.
  --rnet               An igraph network and a netType.
- Optional
  --out                Path where the results to be saved [Default: '.'].
"""

// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

params.out = "."

rgwas = file("$params.rgwas")
rnet = file("$params.rnet")

associationScores = ["chi2"]
modelScores = ["bic", "aic", "consistency"]

process run_evo {

  publishDir "$params.out", overwrite: true

  input:
    each c from associationScores
    each m from modelScores
    file rgwas
    file rnet

  output:
    file "cones.evo.$c.$m.*.RData" into evo_rdata

  """
  #!/usr/bin/env Rscript
  library(martini)
  load("$rgwas")
  load("$rnet")

  start.time <- Sys.time()
  cones <- search_cones(gwas, net, associationScore = "$c", modelScore = "$m")
  end.time <- Sys.time()

  info\$test <- "evo.$c.$m"
  info\$net <- netType
  info\$LD <- LD
  info\$runtime <- end.time - start.time

  save(info, cones, file = paste("cones", info\$test, info\$id, "RData", sep = "."))
  """

}

process run_scones {

  publishDir "$params.out", overwrite: true

  input:
    file rgwas
    file rnet

  output:
    file "cones.scones.*.RData" into evo_chisq_consistency_rdata

  """
  #!/usr/bin/env Rscript
  library(martini)
  load("$rgwas")
  load("$rnet")

  start.time <- Sys.time()
  cones <- find_cones(gwas, net)
  end.time <- Sys.time()

  info\$test <- "scones"
  info\$net <- netType
  info\$LD <- LD
  info\$runtime <- end.time - start.time

  save(info, cones, file = paste("cones", info\$test, info\$id, "RData", sep = "."))
  """
}

process runLasso {

  publishDir "$params.out", overwrite: true

  input:
    file rgwas
    file rnet

  output:
    file "cones.lasso.*.RData" into lasso_rdata

  """
  #!/usr/bin/env Rscript
  library(glmnet)
  library(snpStats)
  load("$rgwas")
  load("$rnet")

  X <- as(gwas\$genotypes, "numeric")
  Y <- gwas\$fam\$affected

  start.time <- Sys.time()
  fit.cv <- cv.glmnet(X, Y, family = "binomial", type.measure = "auc")
  fit <- glmnet(X, Y, lambda = fit.cv\$lambda.1se)
  end.time <- Sys.time()

  cones <- gwas\$map
  cones\$selected <- fit\$beta != 0

  info\$test <- "lasso"
  info\$net <- netType
  info\$LD <- LD
  info\$runtime <- end.time - start.time

  save(info, fit.cv, fit, cones, file = paste("cones", info\$test, info\$id, "RData", sep = "."))
  """

}

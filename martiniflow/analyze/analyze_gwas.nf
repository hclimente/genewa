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

associationScores = ["chi2","skat"]
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
  library(igraph)
  load("$rgwas")
  load("$rnet")

  start.time <- Sys.time()
  cones <- search_cones(gwas, net, associationScore = "$c", modelScore = "$m")
  end.time <- Sys.time()

  detectedGenes <- martini:::subvert(net, 'name', cones\$snp[cones\$selected])\$gene %>% unique

  info\$test <- "evo.$c.$m"
  info\$net <- netType
  info\$LD <- LD
  info\$runtime <- end.time - start.time

  save(info, cones, detectedGenes, file = paste("cones", info\$test, info\$id, "RData", sep = "."))
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
  library(igraph)
  library(glmnet)
  library(snpStats)
  library(martini)
  load("$rgwas")
  load("$rnet")

  X <- as(gwas\$genotypes, "numeric")
  Y <- gwas\$fam\$affected

  start.time <- Sys.time()
  fit.cv <- cv.glmnet(X, Y, family = "binomial", type.measure = "auc")
  fit <- glmnet(X, Y, lambda = fit.cv\$lambda.1se)
  end.time <- Sys.time()

  cones <- gwas\$map
  cones\$selected <- as.logical(fit\$beta != 0)

  detectedGenes <- martini:::subvert(net, 'name', cones\$snp.names[cones\$selected])\$gene %>% unique

  info\$test <- "lasso"
  info\$net <- netType
  info\$LD <- LD
  info\$runtime <- end.time - start.time

  save(info, fit.cv, fit, cones, detectedGenes, file = paste("cones", info\$test, info\$id, "RData", sep = "."))
  """

}

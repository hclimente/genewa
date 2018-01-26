#!/usr/bin/env nextflow

// Help message
helpMessage = """
Launch several analysis on a GWAS.

Usage:

  nextflow martiniflow/analyze/analyze_gwas.nf --rgwas gwas.RData --rnet net.RData

	(rgwas, rnet) -> (cones.evo.*.RData)

INPUT

- rgwas                 snpMatrix with a GWAS experiment and an info list.
- rnet                  An igraph network and a netType.
- associationScore      Association score (chi2 or skat).
- encoding				Genetic model (additive, recessive, dominant, codominant)
- lambda                Lambda value.
- eta                   Eta value.

OUTPUT

- cones.evo.*.RData    RData files containing at least test, info, cones and time.taken.

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

associationScore = params.associationScore
encoding = params.encoding
lambda = params.lambda
eta = params.eta

process run_evo {

  publishDir "$params.out", overwrite: true

  input:
    file rgwas
    file rnet

  output:
    file "cones.evo.${associationScore}.${encoding}.*.RData" into evo_rdata
    file "cones.evo.${associationScore}.${encoding}.*.tsv" into evo_cones

  """
  #!/usr/bin/env Rscript
  library(martini)
  library(igraph)
  load("$rgwas")
  load("$rnet")

  X <- as(gwas\$genotypes, "numeric")
  X <- martini:::encode_gwas(X, $encoding)
  Y <- gwas\$fam\$affected
  W <- igraph::as_adj(net)

  start.time <- Sys.time()
  settings <- martini:::get_evo_settings(etas = 10e$eta, lambdas = 10e$lambda)
  cones <- martini:::evo(X, Y, W, settings)
  end.time <- Sys.time()

  detectedGenes <- martini:::subvert(net, 'name', cones\$snp[cones\$selected])\$gene %>% unique

  info\$test <- "SConES"
  info\$statistic <- "$associationScore"
  info\$selection <- "L10e${lambda}-E10e${eta}"
  info\$net <- netType
  info\$LD <- LD
  info\$runtime <- end.time - start.time

  save(info, cones, detectedGenes, file = paste("cones", info\$test, info\$id, "RData", sep = "."))
  write_tsv(cones, paste("cones", info\$test, info\$id, "RData", sep = "."))
  """

}

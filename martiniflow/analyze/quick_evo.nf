#!/usr/bin/env nextflow

// Help message
helpMessage = """
Launch several analysis on a GWAS.

Usage:

  nextflow martiniflow/analyze/analyze_gwas.nf --rgwas gwas.RData --rnet net.RData

	(rgwas, rnet) -> (cones.evo.*.RData)

INPUT

- rgwas                snpMatrix with a GWAS experiment and an info list.
- rnet                 An igraph network and a netType.
- associationScore		 Association score (chi2 or skat).
- modelScore					 Model score (consistency, bic, aic or aicc)
- encoding						 Genetic model (additive, recessive, dominant, codominant)

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
modelScore = params.modelScore
encoding = params.encoding

process run_evo {

  publishDir "$params.out", overwrite: true

  input:
    file rgwas
    file rnet

  output:
    file "cones.evo.${associationScore}.${modelScore}.${encoding}.*.RData" into evo_rdata

  """
  #!/usr/bin/env Rscript
  library(martini)
  library(igraph)
  load("$rgwas")
  load("$rnet")

  start.time <- Sys.time()
  cones <- search_cones(gwas, net,
                        associationScore = "$associationScore",
                        modelScore = "$modelScore",
                        encoding = "$encoding")
  end.time <- Sys.time()

  detectedGenes <- subvert(net, 'name', cones\$snp[cones\$selected])\$gene %>% unique

  info\$test <- paste0("evo.${associationScore}.${modelScore}.${encoding}.", netType)
  info\$net <- netType
  info\$runtime <- end.time - start.time

  save(info, cones, detectedGenes, file = paste("cones", info\$test, info\$id, "RData", sep = "."))
  """

}

#!/usr/bin/env nextflow

// Help message
helpMessage = """
Read a PED file into a snpMatrix.

Usage: nextflow martiniflow/io/read_gwas.nf --ped genotype.ped --map genotype.map

(ped,map[,causal]) -> (gwas.\$ped.RData[,causal.Rdata])

OUTPUT

- gwas.\$ped.RData. A snpMatrix with the GWAS.
- causal.RData. A logical vector with causal information.

PARAMETERS

- Required
  --ped                       Path of the ped file.
  --map                       Path of the map file.
- Optional:
  --out                       Path where the results to be saved [Default: '.'].
  --causal                    Path to a table with two columns: SNPid and causal (0/1).
"""

// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

ped = file("$params.ped")
map = file("$params.map")
params.out = "."

process readGWAS {

  publishDir "$params.out", overwrite: true, mode: "copy"

  input:
    file ped
    file map
  output:
    file "gwas.${ped}.RData" into gwas_rdata

  """
  #!/usr/bin/env Rscript
  library(snpStats)

  info <- list(
    origin = "experiment",
    file = "$ped",
    id = floor(runif(1, 1, 10000000)))

  gwas <- read.pedfile("$ped", snps = "$map")

  save(gwas, info, file = paste('gwas', '$ped', 'RData', sep = '.'))
  """

}

if (params.hasProperty('causal') && params.causal) {

  causal = file("$params.causal")

  process getCausal {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
      file gwas_rdata
      file causal
    output:
      file "causal.RData" into causal_rdata

      """
      #!/usr/bin/env Rscript
      library(tidyverse)

      load("$gwas_rdata")

      causal <- read_tsv("$causal", col_types = "cd") %>%
        # reorder based of scones reading and remove filtered out snps
        .[match(gwas\$ids, .\$snp),] %>%
        .\$causalSnp %>%
        as.logical

      save(causal, info, file = "causal.RData")
      """

  }
}

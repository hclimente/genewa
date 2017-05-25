ped = file("$params.ped")
map = file("$params.map")
gi = file("$params.gi")
pheno = file("$params.pheno")
params.out = "."

process readGWAS {

  publishDir "$params.out", overwrite: true, mode: "copy"

  input:
    file ped
    file map
    file gi
    file pheno
  output:
    file "gwas.*.RData" into gwas_rdata

  """
  #!/usr/bin/env Rscript
  library(martini)
  library(tidyverse)

  # make consistent names for ped and map
  file.rename("$ped", "genotype.ped")
  file.rename("$map", "genotype.map")

  # read data
  gwas <- readGWAS("genotype", "$pheno", "$gi", 0, 0.05)

  # generate random id for the experiment
  id <- runif(1, 1, 10000000)

  save(gwas, id, file = paste0("gwas.", id, ".RData"))
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
      file "causal.*.RData" into causal_rdata

      """
      #!/usr/bin/env Rscript
      library(tidyverse)

      load("$gwas_rdata")

      causal <- read_tsv("$causal", col_types = "cd") %>%
        # reorder based of scones reading and remove filtered out snps
        .[match(gwas\$ids, .\$snp),] %>%
        .\$causalSnp %>%
        as.logical

      save(causal, file = paste0("causal.", id, ".RData"))
      """

  }
}

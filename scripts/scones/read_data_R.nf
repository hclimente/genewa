ped = file("$params.ped")
map = file("$params.map")
gi = file("$params.gi")
pheno = file("$params.pheno")
truth = file("$params.truth")
params.out = "."

process readData {

  publishDir "$params.out", overwrite: true, mode: "copy"

  input:
    file ped
    file map
    file gi
    file pheno
    file truth
  output:
    file "gwas.*.RData" into gwas_rdata

  """
  #!/usr/bin/env Rscript
  library(rscones2)
  library(tidyverse)

  # make consistent names for ped and map
  file.rename("$ped", "genotype.ped")
  file.rename("$map", "genotype.map")

  # read data
  gwas <- readBio("genotype", "$pheno", "$gi", 0, 0.5)

  # generate random id for the experiment
  id <- runif(1, 1, 10000000)

  truth <- read_tsv("$truth", col_types = "iddccdd") %>%
    .\$causalSnp %>%
    as.integer

  save(gwas, truth, id, file = paste0("gwas.", id, ".RData"))
  """

}

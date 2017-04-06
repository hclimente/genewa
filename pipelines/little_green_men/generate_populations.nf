params.rr = 1.2
params.wd = "."

wd = params.wd
lgmScript = file("/Users/hclimente/projects/genewa/pipelines/little_green_men/little_green_men.nf")
snpNetworkscript = file("/Users/hclimente/projects/genewa/pipelines/scones/get_snp_networks.nf")
getPhenotypesScript = file("/Users/hclimente/projects/genewa/pipelines/scones/get_phenotypes.nf")
analyzePopulationScript = file("/Users/hclimente/projects/genewa/pipelines/little_green_men/analyze_population.nf")
getQualityMeasuresScript = file("/Users/hclimente/projects/genewa/pipelines/scones/getQualityMeasures.R")

process generatePopulation {

  input:
    file lgmScript
    each i from 1..10

  output:
    file "genotypes.ped" into ped, ped2
    file "genotypes.map" into map, map2
    file "gene2snp.tsv" into gene2snp
    file "truth.tsv" into truth
    file "ppi.tab" into ppi
    file "snp_list.csv" into snp_list

  """
  nextflow run $lgmScript --rr $params.rr --population random --N 100
  """
}

process getSconesFiles {

  input:
    file snpNetworkscript
    file getPhenotypesScript
    file gene2snp
    file ped2
    file map2
    file ppi

  output:
    file "gs.txt" into gs
    file "gm.txt" into gm
    file "gi.txt" into gi
    file "phenotype.txt" into pheno

  """
  nextflow run $snpNetworkscript
  nextflow run $getPhenotypesScript
  """

}

process readData{

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
  gwas <- readBio("${ped.baseName}", "$pheno", "$gi", 0, 0.5)
  id <- runif(1, 1, 10000000)

  truth <- read_tsv("$truth", col_types = "iddccdd") %>%
    .\$causalSnp %>%
    as.integer

  save(gwas, truth, id, file = paste0("gwas.", id, ".RData"))
  """
}

process analyzePopulation {

  input:
    file gwas_rdata
    file analyzePopulationScript
    file getQualityMeasuresScript

  output:
    file "*.RData" into analyses

  """
  nextflow run $analyzePopulationScript --rdata $gwas_rdata
  """
}

process joinResults {

  input:
    file '*.RData' from analyses.collect()

  output:
    file "summary.RData" into summary
  """
  #!/usr/bin/env Rscript
  library(magrittr)

  quality <- lapply(list.files(pattern = "*.RData"), function(f){
    load(f)
    qual
    }) %>% do.call("rbind", .)

  save(quality, file = "summary.RData")
  """

}

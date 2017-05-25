#!/usr/bin/env nextflow
// files and working directories
params.wd = "."
params.genewawd = "/Users/hclimente/projects/genewa"

wd = params.wd
genewawd = params.genewawd
readDataRScript = file("$genewawd/scripts/scones/r_read_gwas.nf")
simulatePhenoScript = file("$genewawd/scripts/gwas_general/simulate_phenotypes_from_ped.nf")
snpNetworkscript = file("$genewawd/scripts/scones/get_snp_networks.nf")
getPhenotypesScript = file("$genewawd/scripts/scones/get_phenotypes.nf")
analyzeGWASScript = file("$genewawd/scripts/scones/analyze_gwas.nf")
getQualityMeasuresScript = file("$genewawd/scripts/scones/getQualityMeasures.R")

// real GWAS files
ped = file("$genewawd/${params.geno}.ped")
map = file("$genewawd/${params.geno}.map")
gene2snp = file("$genewawd/$params.gene2snp")
tab = file("$genewawd/$params.tab")

// simulation parameters
h2 = params.h2

process getSconesFiles {

  input:
    file snpNetworkscript
    file getPhenotypesScript
    file gene2snp
    file ped
    file map
    file tab

  output:
    file "gi.txt" into net
    file "phenotype.txt" into pheno

  """
  nextflow run $snpNetworkscript -profile bigmem --tab $tab --map $map --snp2gene $gene2snp
  nextflow run $getPhenotypesScript -profile cluster --ped $ped
  """

}

process readData {

  input:
    file readDataRScript
    file ped
    file map
    file net
    file pheno

  output:
    file "gwas.*.RData" into gwas_rdata

  """
  nextflow run $readDataRScript --ped $ped --map $map --gi $net --pheno $pheno -profile bigmem
  """

}

process simulatePhenotype {

  input:
    each i from 1..10
    file simulatePhenoScript
    file gwas_rdata

  output:
    file "simu*RData" into sgwas_rdata
    file "causal*RData" into scausal_rdata

  """
  nextflow run $simulatePhenoScript --h2 $h2 --n 1283 --gwas $gwas_rdata --i $i -profile bigmem --nAssociatedSnps 20
  """

}

process analyzePopulation {

  input:
    file sgwas_rdata
    file scausal_rdata
    file analyzeGWASScript
    file getQualityMeasuresScript

  output:
    file "*.RData" into analyses

  """
  nextflow run $analyzeGWASScript --rgwas $sgwas_rdata --rcausal $scausal_rdata -profile bigmem -resume
  """

}

process joinResults {

  publishDir "$params.wd", overwrite: true, mode: "copy"

  input:
    file '*.RData' from analyses.collect()

  output:
    file "summary.RData" into summary

  """
  #!/usr/bin/env Rscript
  library(magrittr)
  library(tidyverse)

  quality <- lapply(list.files(pattern = "*.RData"), function(f){
    load(f)
    qual %>%
      mutate(h2 = $h2)
    }) %>% do.call("rbind", .)

  save(quality, file = "summary.RData")
  """

}

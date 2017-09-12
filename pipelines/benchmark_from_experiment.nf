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
snp2gene = file("$genewawd/$params.snp2gene")
tab = file("$genewawd/$params.tab")
net = "$params.net"

// simulation parameters
params.permutations = 10
params.n = 1283

h2 = params.h2
n = params.n
permutations = params.permutations

process readData {

  input:
    file readDataRScript
    file ped
    file map
    file snp2gene
    file tab

  output:
    file "gwas.*.RData" into gwas_rdata
    file "net.*.RData" into net_rdata

  """
  nextflow run $readDataRScript --ped $ped --map $map --net $net --snp2gene $snp2gene --tab $tab -profile bigmem
  """

}

process simulatePhenotype {

  input:
    each k from 1..permutations
    file simulatePhenoScript
    file gwas_rdata
    file net_rdata

  output:
    file "simu*RData" into simgwas_rdata
    file "causal*RData" into simcausal_rdata

  """
  nextflow run $simulatePhenoScript --h2 $h2 --n $n --gwas $gwas_rdata --net $net_rdata --k $k -profile bigmem --nAssociatedSnps 100
  """

}

process analyzePopulation {

  input:
    file simgwas_rdata
    file simcausal_rdata
    file analyzeGWASScript
    file getQualityMeasuresScript

  output:
    file "*.RData" into analyses

  """
  nextflow run $analyzeGWASScript --rgwas $simgwas_rdata --rcausal $simcausal_rdata -profile bigmem -resume
  """

}

process joinResults {

  publishDir "$params.wd", overwrite: true, mode: "copy"

  input:
    file '*.RData' from analyses.collect()

  output:
    file "qualityMeasures.RData" into qualityMeasures

  """
  #!/usr/bin/env Rscript
  library(magrittr)
  library(tidyverse)

  results <- list.files(pattern = "*.RData")

  qMeasures <- lapply(results, function(f){
    load(f)
    qual %>%
      mutate(h2 = $h2,
             time = time.taken)
    }) %>% do.call("rbind", .)

  cones <- lapply(results, function(f){
      load(f)
      cones\$selected
  }) %>% do.call("rbind", .)

  model <- gsub(".RData", "", results)

  save(qMeasures, cones, model, file = "qualityMeasures.RData")
  """

}

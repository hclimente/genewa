params.wd = "."
params.genewawd = "~/projects/genewa"

ped = file("$params.ped")
real_map = file("$params.map")
gene2snp = file("$params.gene2snp")
tab = file("$params.tab")
snp_list = file("$params.snp_list")

wd = params.wd
snpNetworkscript = file("$params.genewawd/pipelines/scones/get_snp_networks.nf")
getPhenotypesScript = file("$params.genewawd/pipelines/scones/get_phenotypes.nf")
readDataRScript = file("$params.genewawd/pipelines/scones/read_data_R.nf")
analyzePopulationScript = file("$params.genewawd/pipelines/little_green_men/analyze_population.nf")
getQualityMeasuresScript = file("$params.genewawd/pipelines/scones/getQualityMeasures.R")

process generateBED {

  container 'gauravkaushik/plink'

  input:
    file ped
    file real_map

  output:
    file "${ped.baseName}.bed" into bed
    file "${ped.baseName}.bim" into bim
    file "${ped.baseName}.fam" into fam

  '''
  plink --file ${ped.baseName} --make-bed --out ${ped.baseName}
  '''

}

process selectCausalSNPs {

  input:
    each i from 1..10

  output:
    file "causal_snps.txt" into causalSnps

  """
  """

}

process generatePopulation {

  container 'biodckrdev/gcta'

  input:
    file bed
    file bim
    file fam
    file causalSnps

  output:

    file "genotypes.map" into map, map2

  """
  gcta64 --bfile ${bed.baseName} --simu-cc 500 500 --simu-causal-loci $causalSnps --simu-hsq 0.5 --simu-k 0.1 --simu-rep 3 --out genotypes.map
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
  nextflow run $snpNetworkscript -profile cluster
  nextflow run $getPhenotypesScript -profile cluster
  """

}

process readData {

  input:
    file readDataRScript
    file ped
    file map
    file gi
    file pheno
    file truth

  output:
    file "gwas.*.RData" into gwas_rdata

  """
  nextflow run $readDataRScript --ped $ped --map $map --gi $gi --pheno $pheno --truth $truth -profile bigmem
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
  nextflow run $analyzePopulationScript --rdata $gwas_rdata -profile bigmem -resume
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
	mutate(rr = $params.rr)
    }) %>% do.call("rbind", .)

  save(quality, file = "summary.RData")
  """

}

params.wd = "."
params.genewawd = "/Users/hclimente/projects/genewa"

ped = file("${params.geno}.ped")
map = file("${params.geno}.map")
gene2snp = file("$params.gene2snp")
tab = file("$params.tab")

wd = params.wd
ped2bedScript = file("$params.genewawd/scripts/gwas_general/ped2bed.nf")
simulatePhenoScript = file("$params.genewawd/scripts/gwas_general/simulate_phenotypes_from_ped.nf")
snpNetworkscript = file("$params.genewawd/scripts/scones/get_snp_networks.nf")
getPhenotypesScript = file("$params.genewawd/scripts/scones/get_phenotypes.nf")
readDataRScript = file("$params.genewawd/scripts/scones/read_data_R.nf")
analyzeGWASScript = file("$params.genewawd/scripts/scones/analyze_gwas.nf")
getQualityMeasuresScript = file("$params.genewawd/scripts/scones/getQualityMeasures.R")

process ped2bed {

  input:
    file ped
    file map
    file ped2bedScript

  output:
    file "${ped.baseName}.bed" into bed
    file "${ped.baseName}.bim" into bim
    file "${ped.baseName}.fam" into fam

  """
  nextflow run $ped2bedScript --gen $ped.baseName
  """

}

process simulatePhenotype {

  input:
    each i from 1..1
    file simulatePhenoScript
    file bed
    file bim
    file fam
    file ped
    file gene2snp
    file tab
    file map

  output:
    file "${ped.baseName}.simupheno.ped" into simuPed, simuPed2
    file "truth.tsv" into truth

  """
  nextflow run $simulatePhenoScript \
    --gen $ped.baseName \
    --gene2snp $gene2snp \
    --tab $tab \
    --h2 0.8
  """

}

process getSconesFiles {

  input:
    file snpNetworkscript
    file getPhenotypesScript
    file gene2snp
    file simuPed
    file map
    file tab

  output:
    file "gi.txt" into gi
    file "phenotype.txt" into pheno

  """
  nextflow run $snpNetworkscript -profile cluster --tab $tab --map $map --snp2gene $gene2snp
  nextflow run $getPhenotypesScript -profile cluster --ped $ped
  """

}

process readData {

  input:
    file readDataRScript
    file simuPed2
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
    file analyzeGWASScript
    file getQualityMeasuresScript

  output:
    file "*.RData" into analyses

  """
  nextflow run $analyzeGWASScript --rdata $gwas_rdata -profile bigmem -resume
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

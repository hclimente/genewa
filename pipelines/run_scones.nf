#! /bin/env nextflow

// scones parameters
association_score = params.association_score
model_selection = params.model_selection

// files and paths
ped = file("${params.gen}.ped")
map = file("${params.gen}.map")
gene2snp = file("$params.gene2snp")
tab = file("$params.tab")

params.wd = "."
params.genewawd = "/Users/hclimente/projects/genewa"

wd = params.wd

// scripts
snpNetworkscript = file("$params.genewawd/scripts/scones/get_snp_networks.nf")
getPhenotypesScript = file("$params.genewawd/scripts/scones/get_phenotypes.nf")
runSconesScript = file("$params.genewawd/scripts/scones/run_scones.nf")

process getSconesFiles {

  input:
    file snpNetworkscript
    file getPhenotypesScript
    file gene2snp
    file ped
    file map
    file tab

  output:
    file "gi.txt" into gi
    file "phenotype.txt" into phen

  """
  nextflow run $snpNetworkscript -profile bigmem --tab $tab --map $map --snp2gene $gene2snp
  nextflow run $getPhenotypesScript -profile cluster --ped $ped
  """

}

process runScones {

  input:
    file runSconesScript
    file ped
    file map
    file gi
    file phen

  output:
    file "pmatrix.txt" into pmatrix
    file "selected_snps.txt" into selected_snps

  """
  nextflow run $runSconesScript -profile bigmem \
    --gen ${ped.baseName}\
    --phen $phen \
    --net $gi \
    --association_score $association_score \
    --model_selection $model_selection \
    --depth 3 \
    --maf 0.05 \
    --lambda -1 \
    --eta -1 \
    --outdir . \
    --encoding additive \
    --pc 0 \
    --seed 0
  """

}

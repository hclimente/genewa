#!/usr/bin/env nextflow

////////////////////////
// GENOME PARAMETERS //
///////////////////////
params.population = "random"
params.numCausalGenes = 5
params.numGenesPerPathway = 6
params.numPathways = 40
params.numSnpsPerGene = 4
params.numChromosomes = 6
params.freqCrosstalk = 0.3

// unused
params.freqCausalSnp = 0.5

///////////////////////////
// POPULATION PARAMETERS //
///////////////////////////
params.N = 1000
params.rr = 3

///////////////////////////
//   SCONES PARAMETERS   //
///////////////////////////
params.stat = "chisq"
params.lambda_base = 10
params.eta_base = 10

//////////////
// SCRIPTS //
/////////////
population_generation_script = file("little_green_men.nf")
run_scones_script = file("../scones/run_scones.nf")
analyze_grid_script = file("analyze_grid.R")
nextflow_cfg = file("nextflow.config")

outdir = "gridsearch/${params.stat}_lambdaBase${params.lambda_base}_etaBase${params.eta_base}_rr${params.rr}_N${params.N}_numCausalGenes${params.numCausalGenes}_numGenesPerPathway${params.numGenesPerPathway}_numPathways${params.numPathways}_numSnpsPerGene${params.numSnpsPerGene}_numChromosomes${params.numChromosomes}_freqCrosstalk${params.freqCrosstalk}_freqCausalSnp${params.freqCausalSnp}"
num_causal_snps = params.numCausalGenes * params.numSnpsPerGene
num_genes = params.numSnpsPerGene * params.numGenesPerPathway * params.numPathways

process generate_population {

  publishDir "$outdir", mode: 'copy'

  input:
    file population_generation_script
    file nextflow_cfg
  output:
    file "genotypes.ped" into ped
    file "genotypes.map" into map
    file "gene2snp.tsv" into gene2snp
    file "snp_list.csv" into snp_list
    file "ppi.tab" into ppi

  """
  nextflow run little_green_men.nf -resume --numCausalGenes $params.numCausalGenes --numGenesPerPathway $params.numGenesPerPathway --numPathways $params.numPathways --numSnpsPerGene $params.numSnpsPerGene --numChromosomes $params.numChromosomes --freqCrosstalk $params.freqCrosstalk --N $params.N --rr $params.rr --population $params.population --freqCausalSnp $params.freqCausalSnp
  """
}

process run_scones {

  publishDir "$outdir/grid", mode: 'copy'

  input:
    file run_scones_script
    file nextflow_cfg
    file ped
    file map
    file gene2snp
    file snp_list
    file ppi
  output:
    file "gridsearch/*.out.ext.txt" into snps
    file "gridsearch/*.terms.txt" into terms
    file "gridsearch/*.L.txt" into laplacian

  """
  nextflow run run_scones.nf -profile curie -resume --ped genotypes.ped --map  genotypes.map --snp2gene  gene2snp.tsv --icogsSnps snp_list.csv --tab ppi.tab --gridsearch custom --stat $params.stat --lambda_base $params.lambda_base --eta_base $params.eta_base
  """
}

process analyze_grid {

  publishDir "$outdir", mode: 'copy'

  input:
    file analyze_grid_script
    file snps
    file terms
    file laplacian
  output:
    file "gs.png" into gs
    file "gs.interactions.png" into gs_intx
    file "gm.png" into gm
    file "gm.interactions.png" into gm_intx
    file "gi.png" into gi
    file "gi.interactions.png" into gi_intx

  """
  Rscript $analyze_grid_script $num_causal_snps $num_genes gs $params.lambda_base $params.eta_base
  Rscript $analyze_grid_script $num_causal_snps $num_genes gm $params.lambda_base $params.eta_base
  Rscript $analyze_grid_script $num_causal_snps $num_genes gi $params.lambda_base $params.eta_base
  """
}

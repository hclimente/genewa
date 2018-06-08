#!/usr/bin/env nextflow

// files and working directories
params.wd = "."
params.genewawd = "/Users/hclimente/projects/genewa"

wd = params.wd
genewawd = params.genewawd
srcReadPed = file("$genewawd/martiniflow/io/read_ped.nf")
srcGetNetwork = file("$genewawd/martiniflow/io/get_network.nf")
srcSimuGWAS = file("$genewawd/martiniflow/simulations/simulate_phenotypes.nf")
srcAnalyzeGWAS = file("$genewawd/martiniflow/analyze/analyze_gwas.nf")
srcBenchmark = file("$genewawd/martiniflow/qc/stats_simulations.nf")
srcQualityMeasures = file("$genewawd/martiniflow/qc/getQualityMeasures.R")

// real GWAS files
ped = file("$genewawd/${params.geno}.ped")
map = file("$genewawd/${params.geno}.map")
snp2gene = file("$genewawd/$params.snp2gene")
tab2 = file("$genewawd/$params.tab2")

// simulation parameters
params.ngenes = 20
params.psnps = 1
params.permutations = 10
params.cases = 636
params.controls = 637
params.prevalence = 0.5
params.rld = "None"

params.heritabilities = [0.25, 1]
params.nets = ["gs", "gm", "gi"]

ngenes = params.ngenes
psnps = params.psnps
permutations = params.permutations
cases = params.cases
controls = params.controls
prevalence = params.prevalence

heritabilities = params.heritabilities
nets = params.nets

process readGWAS {

  input:
    file srcReadPed
    file ped
    file map

  output:
    file "gwas*.RData" into rgwas

  """
  nextflow run $srcReadPed --ped $ped --map $map -profile bigmem
  """

}

rgwas.into { rgwas_getNetwork; rgwas_simulate_net; rgwas_simulate }

process getNetwork {

  input:
    file srcGetNetwork
    val net from nets
    file rgwas_getNetwork
    file snp2gene
    file tab2

  output:
    file "net*.RData" into rnet

  """
  nextflow run $srcGetNetwork --gwas $rgwas_getNetwork --net $net --snp2gene $snp2gene --tab2 $tab2 -profile bigmem
  """

}

if (params.rld != "None") {

  rld = file("$genewawd/$params.rld")

  process getLDNetwork {

    input:
      file srcGetNetwork
      val net from nets
      file rgwas_getNetwork
      file snp2gene
      file tab2
      file rld

    output:
      file "net*.RData" into rldnet

    """
    nextflow run $srcGetNetwork --gwas $rgwas_getNetwork --net $net --snp2gene $snp2gene --tab2 $tab2 --rld $rld -profile bigmem
    """

  }

  rnet.mix(rldnet).set{rnet}

}

process getGI4Simulations {

  input:
    file srcGetNetwork
    file snp2gene
    file tab2
    file rgwas_simulate_net

  output:
    file "net*.RData" into gi

  """
  nextflow run $srcGetNetwork --gwas $rgwas_simulate_net --net gi --snp2gene $snp2gene --tab2 $tab2 -profile bigmem
  """

}

process simulatePhenotype {

  input:
    file srcSimuGWAS
    file gi
    file rgwas_simulate
    each k from 1..permutations
    each h2 from heritabilities

  output:
    set file("gwas.simulated.RData"), file("causal.RData") into simulations

  """
  nextflow run $srcSimuGWAS --rgwas $rgwas_simulate --rnet $gi  --h2 $h2 --cases $cases --controls $controls --ngenes $ngenes --psnps $psnps --prevalence $prevalence -profile bigmem
  """
}

process benchmarkSimulation {

  input:

    file srcAnalyzeGWAS
    file srcBenchmark
    file srcQualityMeasures
    each net from rnet
    set file(rgwas), file(rcausal) from simulations

  output:
    file "cones.RData" into simulationCones
    file "benchmark.RData" into simulationBenchmarks

  """
  nextflow run $srcAnalyzeGWAS --rgwas $rgwas --rnet $net -profile bigmem
  nextflow run $srcBenchmark --rcausal $rcausal -profile cluster
  """

}

process joinBenchmarks {

  publishDir "$params.wd", overwrite: true, mode: "copy"

  input:
    file '*.RData' from simulationBenchmarks.collect()

  output:
    file "benchmark.tsv" into benchmark

  """
  #!/usr/bin/env Rscript
  library(tidyverse)

  results <- list.files(pattern = "*.RData")

  benchmark <- lapply(results, function(f){
    load(f)
    benchmark
  }) %>% do.call("rbind", .)

  write_tsv(benchmark, "benchmark.tsv")
  """

}

process joinCones {

  publishDir "$params.wd", overwrite: true, mode: "copy"

  input:
    file "*.RData" from simulationCones.collect()

  output:
    file "cones.tsv" into cones

  """
  #!/usr/bin/env Rscript
  library(magrittr)
  library(tidyverse)

  results <- list.files(pattern = "*.RData")

  cones <- lapply(results, function(f){
    load(f)
    colnames(cones) <- tests
    cones
  }) %>% do.call("cbind", .) %>%
	as.data.frame

  write_tsv(cones, "cones.tsv")
  """
}

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
ld = file("$params.ld")
tab = file("$genewawd/$params.tab")

// simulation parameters
params.p = 2
params.permutations = 10
params.n = 1283

heritabilities = [0.25, 0.5, 0.75, 1]
nets = ["gs", "gm", "gi"]
n = params.n
p = params.p
permutations = params.permutations

process readGWAS {

  input:
    file srcReadPed
    file ped
    file map

  output:
    file "gwas.RData" into rgwas

  """
  nextflow run $srcReadPed --ped $ped --map $map -profile bigmem
  """

}

rgwas.into { rgwas_getNetwork; rgwas_simulate }

if (params.hasProperty('rld') && params.rld) {

  rld = file("$params.rld")

  process getNetwork {

    input:
      file srcGetNetwork
      val net from nets
      file rgwas_getNetwork
      file snp2gene
      file tab
      file rld

    output:
      file "net.RData" into rnet

    """
    nextflow run $srcGetNetwork --gwas $rgwas_getNetwork --net $net --snp2gene $snp2gene --tab $tab --ld $rld -profile bigmem
    """

  }

} else {

  process getNetwork {

    input:
      file srcGetNetwork
      val net from nets
      file rgwas_getNetwork
      file snp2gene
      file tab

    output:
      file "net.RData" into rnet

    """
    nextflow run $srcGetNetwork --gwas $rgwas_getNetwork --net $net --snp2gene $snp2gene --tab $tab -profile bigmem
    """

  }

}

process benchmarkSimulation {

  input:
    file srcSimuGWAS
    file srcAnalyzeGWAS
    file srcBenchmark
    file srcQualityMeasures
    each k from 1..permutations
    each h2 from heritabilities
    file net from rnet
    file rgwas_simulate

  output:
    file "cones.RData" into cones
    file "benchmark.RData" into simulationBenchmarks

  """
  nextflow run $srcSimuGWAS --rgwas $rgwas_simulate --rnet $net  --h2 $h2 --n $n --nAssociatedSnps $p -profile bigmem
  nextflow run $srcAnalyzeGWAS --rgwas simGwas.RData --rnet $net -profile bigmem
  nextflow run $srcBenchmark --rcausal causal.RData -profile cluster
  """

}

process joinBenchmarks {

  publishDir "$params.wd", overwrite: true, mode: "copy"

  input:
    file '*.RData' from simulationBenchmarks.collect()

  output:
    file "benchmark.RData" into benchmark

  """
  #!/usr/bin/env Rscript
  library(magrittr)
  library(tidyverse)

  results <- list.files(pattern = "*.RData")

  benchmark <- lapply(results, function(f){
    load(f)
    benchmark
  }) %>% do.call("rbind", .)

  save(benchmark, file = "benchmark.RData")
  """

}

#!/usr/bin/env nextflow

// Help message
helpMessage = """
Launch several analysis on a GWAS.

Usage:

  nextflow martiniflow/qc/stats_simulations.nf

(causal.RData, cones.*.RData) -> (benchmark.RData, benchmark.tsv, cones.RData)

INPUT

- causal.RData         Contains the logical vector with the true information per SNP.
- cones*.RData              RData files containing at least test, id, cones and time.taken.

OUTPUT

- benchmark.RData      Contains tibble benchmark with accuracy, elapsed time, etc.
- benchmark.tsv        benchmark in TSV format.
- cones.RData          Contains a logical matrix representing which SNP was picked by each method.

PARAMETERS

- Required
  --rcausal            Logical vector indicating causal SNPs.
- Optional
  --out                Path where the results to be saved [Default: '.'].
"""

// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

params.out = "."
rcausal = file("$params.rcausal")
Channel
  .fromPath("cones*.RData")
  .into { rcones_benchmark; rcones_conesData }

srcQualityMeasures = file("getQualityMeasures.R")

process benchmark {

  publishDir "$params.out", overwrite: true

  input:
    file "cones*.RData" from rcones_benchmark.collect()
    file rcausal
    file srcQualityMeasures

  output:
    file "benchmark.RData" into rbenchmark
    file "benchmark.tsv" into benchmark

  """
  #!/usr/bin/env Rscript
  source("$srcQualityMeasures")
  library(tidyverse)

  load("$rcausal")

  benchmark <- lapply(list.files(pattern = "cones.+RData"), function(f){
    load(f)

    getQualityMeasures(as.numeric(cones\$selected), as.numeric(causal), test) %>%
      mutate(time = time.taken,
             test = test,
             id = as.numeric(info\$id),
             solutionSize = sum(cones\$selected),
             causalSnps = as.numeric(info\$nsnps),
             h2 = as.numeric(info\$h2),
             net = info\$network)
    }) %>% do.call("rbind", .)

  save(benchmark, file = "benchmark.RData")
  write_tsv(benchmark, path = "benchmark.tsv")
  """

}

process getCones {

  publishDir "$params.out", overwrite: true

  input:
    file "cones*.RData" from rcones_conesData.collect()

  output:
    file "cones.RData" into rsumcones

  """
  #!/usr/bin/env Rscript
  library(tidyverse)

  results <- list.files(pattern = "cones.+RData")

  cones <- lapply(results, function(f) {
      load(f)
      cones\$selected
  }) %>% do.call("cbind", .)

  tests <- lapply(results, function(f) {
      load(f)
      test
  }) %>% do.call("c", .)

  save(cones, tests, file = "cones.RData")
  """

}

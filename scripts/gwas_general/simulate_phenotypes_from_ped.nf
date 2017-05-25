// Causal SNP information
nAssociatedSnps = params.nAssociatedSnps

// Phenotype information
h2 = params.h2
n = params.n


// File info
params.i = 1

i = params.i
gwas_rdata = file("$params.gwas")
params.out = "."

process simulatePhenotypes {

  publishDir "$params.out", overwrite: true, mode: "copy"

  input:
    file gwas_rdata

  output:

    file "simu*RData" into sgwas_rdata
    file "causal*RData" into scausal_rdata

  """
  #!/usr/bin/env Rscript
  library(martini)

  load("$gwas_rdata")

  causal <- simulateCausalSNPs(gwas\$net, $nAssociatedSnps)
  simulNum <- $i

  # get their effect sizes from a normal distribution and simulate the phenotype
  effectSizes <- rnorm(sum(causal))
  gwas\$Y <- simulatePhenotype(gwas\$X, causal,
                               h2 = 1,
                               effectSize = effectSizes,
                               qualitative = TRUE, ncases = $n, ncontrols = $n)

  save(gwas, id, file = paste("simu", id, simulNum, "RData", sep = "."))
  save(causal, id, file = paste("causal", id, simulNum, "RData", sep = "."))
  """
}

// Causal SNP information
nAssociatedSnps = params.nAssociatedSnps

// Phenotype information
h2 = params.h2
n = params.n

// File info
params.k = 1

k = params.k
gwas_rdata = file("$params.gwas")
net_rdata = file("$params.net")
params.out = "."

process simulatePhenotypes {

  publishDir "$params.out", overwrite: true, mode: "copy"

  input:
    file gwas_rdata
    file net_rdata

  output:

    file "simu*RData" into sgwas_rdata
    file "causal*RData" into scausal_rdata

  """
  #!/usr/bin/env Rscript
  library(martini)

  load("$gwas_rdata")
  load("$net_rdata")

  causalSnps <- simulate_causal_snps(gwas, net, $nAssociatedSnps)
  k <- $k

  # get their effect sizes from a normal distribution and simulate the phenotype
  effectSizes <- rnorm(length(causalSnps))
  gwas\$fam\$affected <- simulate_phenotype(gwas, causalSnps,
                                            h2 = $h2,
                                            effectSize = effectSizes,
                                            qualitative = TRUE,
                                            ncases = $n, ncontrols = $n)

  causal <- gwas\$map\$snp.names %in% causalSnps
  save(gwas, net, id, file = paste("simu", id, k, "RData", sep = "."))
  save(causal, id, file = paste("causal", id, k, "RData", sep = "."))
  """
}

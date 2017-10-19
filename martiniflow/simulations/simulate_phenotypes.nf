#!/usr/bin/env nextflow

// Help message
helpMessage = """
Simulate a set of causal SNPs interconnected in the network and phenotypes.

Usage:

  nextflow martiniflow/simulations/simulate_phenotypes.nf --rgwas gwas.RData --rnet net.RData --ngenes 20 --h2 1 --n 1283

(rgwas, rnet, ngenes, h2, n) -> (simu*RData, causal*RData)

INPUT

- rgwas                snpMatrix with a GWAS experiment and an info tag.
- rnet                 An igraph network and a netType.

OUTPUT

- causal.RData           Logical vector with the SNPs selected as causal. Contains a vector with the unique info of the simulation.
- simGwas.RData          snpMatrix identical to rgwas, with the new simulated phenotypes. Contains a vector with the unique info of the simulation.

PARAMETERS

- Required
  --rgwas                snpMatrix with a GWAS experiment.
  --rnet                 An igraph network and a netType.
  --ngenes               Number of causal genes.
  --psnps                Proportion of causal SNPs in the genes.
  --h2                   Heritability for the simulation.
  --cases                Number of cases.
  --controls             Number of controls.
	--prevalence           Prevalence of the disease in the population.
- Optional:
  --out                  Path where the results to be saved [Default: '.'].
"""

// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

// Causal SNP information
ngenes = params.ngenes
psnps = params.psnps

// Phenotype information
h2 = params.h2
cases = params.cases
controls = params.controls
prevalence = params.prevalence

// File info
rgwas = file("$params.rgwas")
rnet = file("$params.rnet")
params.out = "."

process simulatePhenotypes {

  publishDir "$params.out", overwrite: true, mode: "copy"

  input:
    file rgwas
    file rnet

  output:

    file "simGwas.RData" into srgwas
    file "causal.RData" into scausal_rdata

  """
  #!/usr/bin/env Rscript
  library(martini)
  library(igraph)

  load("$rgwas")
  load("$rnet")

  causalSnps <- simulate_causal_snps(net, $ngenes, $psnps)

  # get their effect sizes from a normal distribution and simulate the phenotype
  effectSizes <- rnorm(length(causalSnps))
  gwas <- simulate_phenotype(gwas, causalSnps,
                             h2 = $h2,
                             effectSize = effectSizes,
                             qualitative = TRUE,
                             ncases = $cases,
														 ncontrols = $controls,
														 prevalence = $prevalence)

  # generate random info for the simulation
  info <- list(origin = "simulation",
               id = floor(runif(1, 1, 10000000)),
               numCausalGenes = $ngenes,
               proportionCausalSnps = $psnps,
               realSolutionSize = length(causalSnps),
               h2 = $h2)

  causal <- gwas\$map\$snp.names %in% names(causalSnps)
  causalGenes <- V(net)\$gene[causalSnps] %>% unique
  save(gwas, info, netType, file = "simGwas.RData")
  save(causal, causalGenes, info, file = "causal.RData")
  """
}

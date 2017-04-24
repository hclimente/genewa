// Causal SNP information
params.idx = 1
params.nAssociatedSnps = 20

idx = params.idx
nAssociatedSnps = params.nAssociatedSnps

// Phenotype information
h2 = params.h2
n = 1283
// balanced datasets
prevalence = 0.5

// File info
bed = file("${params.gen}.bed")
bim = file("${params.gen}.bim")
fam = file("${params.gen}.fam")
ped = file("${params.gen}.ped")
map = file("${params.gen}.map")

gene2snp = file("$params.gene2snp")
tab = file("$params.tab")
params.out = "."

process selectCausalSNPs {

  input:
    file tab
    file gene2snp

  output:
    file "causal_snps.txt" into causalSnps

  """
  #!/usr/bin/env Rscript
  library(magrittr)
  library(tidyverse)
  library(igraph)

  gene2snp <- read_tsv("$gene2snp")

  # select a group of interconnected genes
  ppi <- read_tsv("$tab") %>%
    select(OFFICIAL_SYMBOL_FOR_A, OFFICIAL_SYMBOL_FOR_B) %>%
    filter(OFFICIAL_SYMBOL_FOR_A %in% gene2snp\$GENE & OFFICIAL_SYMBOL_FOR_B %in% gene2snp\$GENE) %>%
    graph_from_data_frame(directed = FALSE)

  cliques <- largest_cliques(ppi)
  clique <- cliques[[$idx]]

  # randomly select their snps
  causalSnps <- gene2snp %>%
    filter(GENE %in% names(clique)) %>%
    sample_n($nAssociatedSnps)

  causalSnps %>%
    select(SNP) %>%
    write_tsv("causal_snps.txt")
  """

}

process simulatePhenotypes {

  input:
    file bed
    file bim
    file fam
    file causalSnps

  output:

    file "simu.phen" into pheno
    file "simu.par" into snpInfo

  """
  gcta64 --bfile ${bed.baseName} \
         --simu-cc $n $n \
         --simu-causal-loci $causalSnps \
         --simu-hsq $h2 \
         --simu-k $prevalence \
         --out simu
  """
}

process getSimuPed {

  publishDir "$params.out", overwrite: true, mode: "copy"

  input:
    file pheno
    file ped

  output:
    file "${ped.baseName}.simupheno.ped" into simuPed

  """
  cut -d' ' -f1-5 $ped >sample.info
  cut -d' ' -f3 $pheno >pheno.info
  cut -d' ' -f7- $ped >geno.info

  paste -d ' ' sample.info pheno.info geno.info >${ped.baseName}.simupheno.ped
  """

}

process getSnpInfo {

  publishDir "$params.out", overwrite: true, mode: "copy"

  input:
    file map
    file snpInfo

  output:
    file "truth.tsv" into truth

  """
  #!/usr/bin/env Rscript
  library(magrittr)
  library(tidyverse)

  causal <- read_tsv("$snpInfo") %>% .\$QTL

  read_tsv("$map", col_names = FALSE) %>%
    set_colnames(c("chr","snp","cm","pos")) %>%
    mutate(causalSnp = as.numeric(snp %in% causal)) %>%
    select(snp, causalSnp) %>%
    write_tsv("truth.tsv")
  """

}

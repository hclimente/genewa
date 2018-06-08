#!/usr/bin/env nextflow

// Help message
helpMessage = """
Create a GS, GM or GI network.

Usage:
  nextflow martiniflow/io/get_network.nf --gwas gwas.RData --net gs
  nextflow martiniflow/io/get_network.nf --gwas gwas.RData --net gm --snp2gene snp2gene.txt
  nextflow martiniflow/io/get_network.nf --gwas gwas.RData --net gi --snp2gene snp2gene.txt --tab2 BIOGRID-ORGANISM-Homo_sapiens-3.4.138.tab2.txt

(gwas,snp2gene,tab2,net) -> (net.Rdata)

OUTPUT:

- net.RData. Contains:
  - net         An igraph with the network.
  - netType     A character variable containing the network type.

PARAMETERS:

- Required
  --gwas        Path to a GWAS experiment in a snpMatrix.
  --snp2gene    Path to a file containing snp2gene mapping.
  --tab2        Path to a PPI information table in TAB2 format.
  --net         Type of network to generate (gs, gm, gi).
  --prune       Remove SNPs not mapped to a gene in a PPI?
- Optional
  --rld         Path to a file with LD information.
  --out         Path where the results to be saved [Default: '.']. Outputs:

"""

// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

params.rld = "None"
params.prune = "FALSE"

rgwas = file("$params.gwas")
snp2gene = file("$params.snp2gene")
tab = file("$params.tab")
net = "$params.net"
prune = "$params.prune"
params.out = "."

process getNetwork {

  publishDir "$params.out", overwrite: true, mode: "copy"

  input:
    file rgwas
    file tab2
    file snp2gene
  output:
    file "net.RData" into rnet

  """
  #!/usr/bin/env Rscript
  library(martini)
  library(tidyverse)
  library(igraph)

  load("$rgwas")
  netType <- "$net"
  LD <- FALSE

  if (netType == "gs") {
    net <- get_GS_network(gwas)
  } else if (netType == "gm") {
    snp2gene <- read_tsv("$snp2gene")  %>%
        rename(snp = SNP, gene = GENE)
    net <- get_GM_network(gwas, snpMapping = snp2gene)
  } else if (netType %in% c('gi', 'gi2', 'ppi')) {
      snp2gene <- read_tsv("$snp2gene") %>%
        rename(snp = SNP, gene = GENE)
      tab2 <- read_tsv("$tab2") %>%
        rename(gene1 = `Official Symbol Interactor A`, gene2 = `Official Symbol Interactor B`) %>%
        select(gene1, gene2)

      if (netType == "gi") {
        net <- get_GI_network(gwas, snpMapping = snp2gene, ppi = tab2)
      } else if (netType == "gi2") {
        gi <- get_GI_network(gwas, snpMapping = snp2gene, ppi = tab2)
        gm <- get_GM_network(gwas, snpMapping = snp2gene)
        gs <- get_GS_network(gwas)

        net <- gi - gm + gs
        net <- set_edge_attr(net, "weight", value = 1)
      } else if (netType == "ppi") {
        gi <- get_GI_network(gwas, snpMapping = snp2gene, ppi = tab2)
        gm <- get_GM_network(gwas, snpMapping = snp2gene)

        net <- gi - gm
        net <- set_edge_attr(net, "weight", value = 1)
      }

    if ("true" == "$prune") {
        ppiGenes <- unique(c(tab2\$gene1, tab2\$gene2))
        net <- martini:::subnet(net, "gene", ppiGenes)
    }

  } else {
    stop("network type not recognized.")
  }

  save(net, netType, LD, file = "net.RData")
  """

}

if (params.rld != "None") {

  rld = file("$params.rld")

  process ld_weight {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
      file rld
      file rnet

    output:
      file "net*.RData" into rldnet

    """
    #!/usr/bin/env Rscript
    library(martini)
    library(tidyverse)

    load("$rnet")
    load("$rld")

    net <- ldweight_edges(net, ld, method = "sigmoid")
    LD <- TRUE

    save(net, netType, LD, file = paste('net', netType, 'RData', sep = '.'))
    """

  }

}

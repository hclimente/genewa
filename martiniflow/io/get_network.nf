#!/usr/bin/env nextflow

// Help message
helpMessage = """
Create a GS, GM or GI network.

Usage:
  nextflow martiniflow/io/get_network.nf --gwas gwas.RData --net gs
  nextflow martiniflow/io/get_network.nf --gwas gwas.RData --net gm --snp2gene snp2gene.txt
  nextflow martiniflow/io/get_network.nf --gwas gwas.RData --net gi --snp2gene snp2gene.txt --tab BIOGRID-ORGANISM-Homo_sapiens-3.4.138.tab.txt

(gwas,snp2gene,tab,net) -> (net.Rdata)

OUTPUT:

- net.RData. Contains:
  - net         An igraph with the network.
  - netType     A character variable containing the network type.

PARAMETERS:

- Required
  --gwas        Path to a GWAS experiment in a snpMatrix.
  --snp2gene    Path to a file containing snp2gene mapping.
  --tab         Path to a PPI information table in TAB format.
  --net         Type of network to generate (gs, gm, gi).
- Optional
  --out         Path where the results to be saved [Default: '.']. Outputs:

"""

// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

gwas_rdata = file("$params.gwas")
snp2gene = file("$params.snp2gene")
tab = file("$params.tab")
net = "$params.net"
params.out = "."

process getNetwork {

  publishDir "$params.out", overwrite: true, mode: "copy"

  input:
    file gwas_rdata
    file tab
    file snp2gene
  output:
    file "net.RData" into net_rdata

  """
  #!/usr/bin/env Rscript
  library(martini)
  library(tidyverse)

  load("$gwas_rdata")
  netType <- "$net"

  if (netType == "gs") {
    net <- get_GS_network(gwas)
  } else if (netType == "gm") {
    snp2gene <- read_tsv("$snp2gene")
    net <- get_GM_network(gwas, snpMapping = snp2gene)
  } else if (netType == "gi") {
    snp2gene <- read_tsv("$snp2gene")
    tab <- read_tsv("$tab") %>%
      select(OFFICIAL_SYMBOL_FOR_A, OFFICIAL_SYMBOL_FOR_B)
    net <- get_GI_network(gwas, snpMapping = snp2gene, ppi = tab)
  } else {
    stop("network type not recognized.")
  }

  save(net, netType, file = "net.RData")
  """

}

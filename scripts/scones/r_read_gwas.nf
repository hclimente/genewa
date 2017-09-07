ped = file("$params.ped")
map = file("$params.map")
snp2gene = file("$params.snp2gene")
tab = file("$params.tab")
net = "$params.net"
params.out = "."

process readGWAS {

  publishDir "$params.out", overwrite: true, mode: "copy"

  input:
    file ped
    file map
  output:
    file "gwas.*.RData" into gwas_rdata

  """
  #!/usr/bin/env Rscript
  library(martini)
  library(snpStats)
  library(tidyverse)

  gwas <- read.pedfile("$ped", snps = "$map")

  # generate random id for the experiment
  id <- floor(runif(1, 1, 10000000))

  save(gwas, id, file = paste0("gwas.", id, ".RData"))
  """

}

process getNetwork {

  publishDir "$params.out", overwrite: true, mode: "copy"

  input:
    file gwas_rdata
    file tab
    file snp2gene
  output:
    file "net.*.RData" into net_rdata

  """
  #!/usr/bin/env Rscript
  library(martini)
  library(tidyverse)

  load("$gwas_rdata")

  if ("$net" == "gs") {
    net <- get_GS_network(gwas)
  } else if ("$net" == "gm") {
    snp2gene <- read_tsv("$snp2gene")
    net <- get_GM_network(gwas, snpMapping = snp2gene)
  } else if ("$net" == "gi") {
    snp2gene <- read_tsv("$snp2gene")
    tab <- read_tsv("$tab") %>%
      select(OFFICIAL_SYMBOL_FOR_A, OFFICIAL_SYMBOL_FOR_B)
    net <- get_GI_network(gwas, snpMapping = snp2gene, ppi = tab)
  } else {
    error("Option not recognized.")
  }

  save(net, id, file = paste0("net.", id, ".RData"))
  """

}

if (params.hasProperty('causal') && params.causal) {

  causal = file("$params.causal")

  process getCausal {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
      file gwas_rdata
      file causal
    output:
      file "causal.*.RData" into causal_rdata

      """
      #!/usr/bin/env Rscript
      library(tidyverse)

      load("$gwas_rdata")

      causal <- read_tsv("$causal", col_types = "cd") %>%
        # reorder based of scones reading and remove filtered out snps
        .[match(gwas\$ids, .\$snp),] %>%
        .\$causalSnp %>%
        as.logical

      save(causal, file = paste0("causal.", id, ".RData"))
      """

  }
}

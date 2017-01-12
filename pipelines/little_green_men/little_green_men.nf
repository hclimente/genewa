#!/usr/bin/env nextflow

params.N = 100
params.rr = 1.5

snpInvolved = 'chromosomes[["involved"]] <- list(list(1,3,c(2)),list(2,4,c(1,5)))'
snpNeutral = 'chromosomes[["neutral"]] <- list(list(1,2,c(4)),list(1,2,c(3,5)),list(2,4,c(4,2)))'

process get_ped {

  publishDir "../../data/little_green_men", overwrite: true

  output:
    file "genotypes.ped" into ped

  """
  #!/usr/bin/env Rscript
  library(readr)
  library(magrittr)
  library(dplyr)

  rr = $params.rr
  N = $params.N
  maf = 0.4

  chromosomes <- list()
  $snpInvolved
  $snpNeutral

  # involved snps
  # A = low risk allele
  # C = high risk allele
  involvedSnps <- list()
  i <- 1
  involvedSnps <- list()
  for (gene in chromosomes[['involved']]){
    for (n in 1:gene[[2]]){
      cases <- rep("A",N)
      cases[runif(N) < (rr/(rr+1))] <- "C"
      controls <- rep("A",N)
      controls[runif(N) > (rr/(rr+1))] <- "C"

      involvedSnps[[i]] <- c(cases, controls)
      i <- i + 1
    }
  }

  involvedSnps <- do.call("cbind", involvedSnps)

  # neutral snps
  # T,G = alleles without impact
  i <- 1
  neutralSnps <- list()
  for (gene in chromosomes[['neutral']]){
    for (n in 1:gene[[2]]){
      cases <- rep("G",N)
      cases[runif(N) < (rr/(rr+1))] <- "T"
      controls <- rep("G",N)
      controls[runif(N) > (rr/(rr+1))] <- "T"

      neutralSnps[[i]] <- c(cases, controls)
      i <- i + 1
    }
  }

  neutralSnps <- do.call("cbind", neutralSnps)

  snps <- cbind(involvedSnps,neutralSnps)
  ids <- 1:(2*N)
  data.frame(family = ids, indiv = ids,
             pater = as.integer(0), mater = as.integer(0),
             pheno = as.integer(c(rep(2, N), rep(1, N)))) %>%
    cbind(snps) %>%
    write_delim("genotypes.ped", col_names=FALSE)

  """

}

process get_map {

  publishDir "../../data/little_green_men", overwrite: true

  output:
    file "genotypes.map" into map

  """
  #!/usr/bin/env Rscript
  library(readr)
  library(magrittr)
  library(dplyr)

  chromosomes <- list()
  $snpInvolved
  $snpNeutral

  perchr <- list("1" = 0, "2" = 0)
  map <- list()
  step <- 5000
  i <- 1
  for (snpType in chromosomes){
    for (gene in snpType){
      chr <- as.character(gene[[1]])
      nsnps <- gene[[2]]
      for (snp in 1:nsnps){
        map[[i]] <- data.frame(chr = chr, position = step + perchr[[chr]]*step, name = paste0("rs", i))
        i <- i + 1
        perchr[[chr]] <- perchr[[chr]] + 1
      }
    }
  }

  do.call("rbind", map) %>%
    mutate(geneticDistance = 0) %>%
    select(chr, name, geneticDistance, position) %>%
    write.table("genotypes.map", sep = "\t", row.names = F, col.names = F, quote = F)

  """

}

process get_phenotypes {

  publishDir "../../data/little_green_men", overwrite: true

  input:
    file ped
  output:
    file "phenotype.txt" into pheno

  """
  echo FID IID TOY > phenotype.txt
  cut -d' ' -f1,2,6 $ped >> phenotype.txt
  """
}

process get_gene2snp {

  publishDir "../../data/little_green_men", overwrite: true

  output:
    file "gene2snp.tsv" into gene2snp

  """
  #!/usr/bin/env Rscript

  library(readr)
  library(magrittr)
  library(dplyr)

  chromosomes <- list()
  $snpInvolved
  $snpNeutral

  snp2gene <- list()
  i <- 1
  g <- 1
  for (snpType in chromosomes){
    for (gene in snpType){
      for (ind in 1:gene[[2]]){
        snp2gene[[i]] <- data.frame(SNP = paste0("rs", i), GENE = paste0("g", g))
        i <- i + 1
      }
      g <- g + 1
    }
  }

  do.call("rbind", snp2gene) %>%
    write_tsv("gene2snp.tsv")

  """

}

process get_ppi {

  publishDir "../../data/little_green_men", overwrite: true

  output:
    file "ppi.tab" into ppi

  """
  #!/usr/bin/env Rscript
  library(readr)
  library(magrittr)
  library(dplyr)

  chromosomes <- list()
  $snpInvolved
  $snpNeutral

  ppi <- list()
  g <- 1
  i <- 1
  for (snpType in chromosomes){
    for (gene in snpType){
      for (intx in gene[[3]]){
        ppi[[i]] <- data.frame(OFFICIAL_SYMBOL_FOR_A = paste0("g", g), OFFICIAL_SYMBOL_FOR_B = paste0("g", intx))
        i <- i + 1
      }
      g <- g + 1
    }
  }

  do.call("rbind", ppi) %>%
    mutate(INTERACTOR_A = "", INTERACTOR_B = "", ALIASES_FOR_A = "",
           ALIASES_FOR_B = "", EXPERIMENTAL_SYSTEM = "", SOURCE = "",
           PUBMED_ID = "", ORGANISM_A_ID = "", ORGANISM_B_ID = "") %>%
    write_tsv("ppi.tab")

  """

}

process get_snp_list {

  publishDir "../../data/little_green_men", overwrite: true

  input:
    file map
  output:
    file "snp_list.csv" into snp_list

  """
  #!/usr/bin/env Rscript
  library(readr)
  library(magrittr)
  library(dplyr)

  read_tsv("$map", col_names=FALSE) %>%
    set_colnames(c("Chromosome","Illumina_SNP_Name","geneDist","Build37_Position")) %>%
    select(Illumina_SNP_Name,Chromosome,Build37_Position) %>%
    write_csv("snp_list.csv")

  """

}

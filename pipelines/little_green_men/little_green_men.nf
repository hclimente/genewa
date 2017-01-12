#!/usr/bin/env nextflow

params.N = 100
params.rr = 1.5

process get_ped {

  publishDir "$HOME/genewa/data/little_green_men", overwrite: true

  output:
    file "genotypes.ped" into ped

  """
  #!/usr/bin/env Rscript
  library(readr)
  library(magrittr)
  library(dplyr)

  rr = $params.rr
  N = $params.N
  nInvolvedSnps = 7
  nNeutralSnps = 8
  maf = 0.4

  # involved snps
  # A = low risk allele
  # C = high risk allele
  involvedSnps <- list()
  for (i in 1:nInvolvedSnps){
    cases <- rep("A",N)
    cases[runif(N) < (rr/(rr+1))] <- "C"
    controls <- rep("A",N)
    controls[runif(N) > (rr/(rr+1))] <- "C"

    involvedSnps[[i]] <- c(cases, controls)
  }

  involvedSnps <- do.call("cbind", involvedSnps)

  # neutral snps
  # T,G = alleles without impact
  neutralSnps <- list()
  for (i in 1:nNeutralSnps){
    cases <- rep("G",N)
    cases[runif(N) < maf] <- "T"
    controls <- rep("G",N)
    controls[runif(N) < maf] <- "T"

    neutralSnps[[i]] <- c(cases, controls)
  }

  neutralSnps <- do.call("cbind", neutralSnps)

  snps <- cbind(involvedSnps,neutralSnps)
  ids <- 1:(2*N)
  data.frame(family = ids, indiv = ids, pater = 0, mater = 0,
             pheno = c(rep(2, N), rep(1, N))) %>%
    cbind(snps) %>%
    write_delim("genotypes.ped", col_names=FALSE)

  """

}

process get_map {

  publishDir "$HOME/genewa/data/little_green_men", overwrite: true

  output:
    file "genotypes.map" into map

  """
  #!/usr/bin/env Rscript
  library(readr)
  library(magrittr)
  library(dplyr)

  chromosomes <- list()
  chromosomes[["involved"]] <- list(c(1,3),c(2,4))
  chromosomes[["neutral"]] <- list(c(1,2),c(1,2),c(2,4))

  perchr <- list("1" = 0, "2" = 0)
  map <- list()
  step <- 5000
  i <- 1
  for (snpType in chromosomes){
    for (gene in snpType){
      chr <- as.character(gene[1])
      nsnps <- gene[2]
      start <- step + perchr[[chr]]*step
      end <- start + nsnps*step
      map[[i]] <- data.frame(chr = chr, position = seq(start, end, step), name = paste0("rs", i))
      perchr[[chr]] <- perchr[[chr]] + nsnps
      i <- i + 1
    }
  }

  do.call("rbind", map) %>%
    mutate(geneticDistance = 0) %>%
    select(chr, name, geneticDistance, position) %>%
    write_delim("genotypes.map", col_names=FALSE)
  """

}

process get_phenotypes {

  publishDir "$HOME/genewa/data/little_green_men", overwrite: true

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

  publishDir "$HOME/genewa/data/little_green_men", overwrite: true

  output:
    file "gene2snp.tsv" into gene2snp

  """
  #!/usr/bin/env Rscript

  library(readr)
  library(magrittr)
  library(dplyr)

  chromosomes <- list()
  chromosomes[["involved"]] <- list(list(1,3,c(2),list(2,4,c(1,5))
  chromosomes[["neutral"]] <- list(list(1,2,c(3),list(1,2,c(2,5),list(2,4,c(4,1))

  snp2gene <- list()
  i <- 1
  g <- 1
  for (snpType in chromosomes){
    for (gene in snpType){
      for (ind in 1:gene[2]){
        snp2gene[[i]] <- data.frame(SNP = paste0("rs", i), GENE = paste0("g", g))
        i <- i + 1
      }
    }
    g <- g + 1
  }

  do.call("rbind", snp2gene) %>%
    write_delim("gene2snp.tsv")

  """

}

process get_ppi {

  publishDir "$HOME/genewa/data/little_green_men", overwrite: true

  output:
    file "ppi.tab" into ppi

  """
  #!/usr/bin/env Rscript
  library(readr)
  library(magrittr)
  library(dplyr)

  chromosomes <- list()
  chromosomes[["involved"]] <- list(list(1,3,c(2),list(2,4,c(1,5))
  chromosomes[["neutral"]] <- list(list(1,2,c(3),list(1,2,c(2,5),list(2,4,c(4,1))

  ppi <- list()
  g <- 1
  for (snpType in chromosomes){
    for (gene in snpType){
      for (intx in gene[3]){
        ppi[[i]] <- data.frame(OFFICIAL_SYMBOL_FOR_A = paste0("g", g), OFFICIAL_SYMBOL_FOR_B = paste0("g", intx))
      }
      g <- g + 1
    }
  }

  do.call("rbind", ppi) %>%
    mutate(INTERACTOR_A = "", INTERACTOR_B = "", ALIASES_FOR_A = "",
           ALIASES_FOR_B = "", EXPERIMENTAL_SYSTEM = "", SOURCE = "",
           PUBMED_ID = "", ORGANISM_A_ID = "", ORGANISM_B_ID = "") %>%
    write_delim("ppi.tab")

  """

}

process get_ppi {

  publishDir "$HOME/genewa/data/little_green_men", overwrite: true

  input:
    file map
  output:
    file "snp_list.csv" into snp_list

  """
  #!/usr/bin/env Rscript
  library(readr)
  library(magrittr)
  library(dplyr)

  read_delim("$map", " ", col_names=FALSE) %>%
    set_colnames(c("Chromosome","Illumina_SNP_Name","geneDist","Build37_Position")) %>%
    select(Illumina_SNP_Name,Chromosome,Build37_Position) %>%
    write_csv("snp_list.csv")

  """

}

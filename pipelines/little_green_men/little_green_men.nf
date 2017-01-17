#!/usr/bin/env nextflow

// params for genotype generation
params.N = 3000
params.rr = 30

// params for genome generation
params.numCausalGenes = 3
params.numGenesPerPathway = 6
params.numPathways = 30
params.numSnpsPerGene = 4
params.numChromosomes = 5
params.crosstalkFreq = 0.3

if (params.population == "custom"){

  process getCustomGenome {

    output:
      file "genome.Rdata" into genome_ped, genome_map, genome_ppi

    """
    #!/usr/bin/env Rscript
    genome <- data.frame(id = 1:5, chr = c(1,1,1,2,2),
                         causalGene = c(F,T,F,F,T),
                         snps = c(2,3,2,4,4))
    genome\$ppi <- list(c(2), c(1,5), c(4), c(3,5), c(4,2))

    save(genome, file = "genome.Rdata")
    """
  }
} else if (params.population == "random") {

  process generateRandomGenome {

    output:
      file "genome.Rdata" into genome_ped, genome_map, genome_ppi

    """
    #!/usr/bin/env Rscript
    library(magrittr)
    library(dplyr)

    # define parameters
    numCausalGenes <- $params.numCausalGenes
    numGenesPerPathway <- $params.numGenesPerPathway
    numPathways <- $params.numPathways
    numSnpsPerGene <- $params.numSnpsPerGene
    numChromosomes <- $params.numChromosomes
    crosstalkFreq <- $params.crosstalkFreq

    # create data frame
    numGenes <- numPathways * numGenesPerPathway
    causalGenes <- sample(1:numGenes, numCausalGenes)
    genome <- data.frame(id = 1:numGenes,
                        chr = sample(numChromosomes, numGenes, replace = TRUE)) %>%
      mutate(causalGene = id %in% causalGenes,
             snps = numSnpsPerGene)

    #calculate ppi
    ppiPathway.model <- expand.grid(1:numGenesPerPathway, 1:numGenesPerPathway)

    ppiPathway <- lapply(0:(numPathways - 1), function(i){
      ppiPathway.model + i * numGenesPerPathway
    }) %>% do.call("rbind", .) %>%
      set_colnames(c("Gene1", "Gene2"))

    ppiCrosstalk <- lapply(1:numPathways, function(i){
      first <- 1 + numGenesPerPathway * (i - 1)
      last <- numGenesPerPathway * i

      pwGenes <- genome\$id[first:last]
      otherGenes <- genome\$id[-(first:last)]

      possibleCrosstalk <- expand.grid(pwGenes, otherGenes)

      sample_n(possibleCrosstalk, floor(numGenesPerPathway * crosstalkFreq))
    }) %>% do.call("rbind", .) %>%
      set_colnames(c("Gene1", "Gene2"))

    ppiCrosstalk <- data.frame(Gene1 = c(ppiCrosstalk\$Gene1, ppiCrosstalk\$Gene2),
               Gene2 = c(ppiCrosstalk\$Gene2, ppiCrosstalk\$Gene1))

    ppi <- rbind(ppiPathway, ppiCrosstalk)

    genome\$ppi <- apply(genome, 1, function(x){
      ppi\$Gene2[ppi\$Gene1 == x[[1]]] %>% sort %>% unique
    })

    # save results
    save(genome, file = "genome.Rdata")
    """
  }

}

process get_ped {

  publishDir "../../data/little_green_men", overwrite: true

  input:
    file genome_ped
  output:
    file "genotypes.ped" into ped

  """
  #!/usr/bin/env Rscript
  library(readr)
  library(magrittr)
  library(dplyr)

  getGenotypes <- function(i, wt, alt, maf, N){
    cases.chr1 <- rep(wt,N)
    cases.chr1[runif(N) < maf] <- alt
    cases.chr2 <- rep(wt,N)
    cases.chr2[runif(N) < maf] <- alt

    controls.chr1 <- rep(wt,N)
    controls.chr1[runif(N) > maf] <- alt
    controls.chr2 <- rep(wt,N)
    controls.chr2[runif(N) > maf] <- alt

    data.frame(chr1 = c(cases.chr1, controls.chr1),
               chr2 = c(cases.chr2, controls.chr2))
  }

  rr = $params.rr
  N = $params.N
  maf = 0.4

  load("$genome_ped")

  # involved snps
  # A = low risk allele
  # C = high risk allele
  involvedSnps <- lapply(1:sum(genome\$snps[genome\$causalGene]), getGenotypes, "A", "C", rr/(rr+1), N) %>%
    do.call("cbind", .)

  # neutral snps
  # T,G = alleles without impact
  neutralSnps <- lapply(1:(sum(genome\$snps) - sum(genome\$snps[genome\$causalGene])), getGenotypes, "T", "G", maf, N) %>%
    do.call("cbind", .)

  snps <- cbind(involvedSnps,neutralSnps)
  ids <- 1:(2*N)
  data.frame(family = ids, indiv = ids,
             pater = as.integer(0), mater = as.integer(0),
             sex = as.integer(2),
             pheno = as.integer(c(rep(2, N), rep(1, N)))) %>%
    cbind(snps) %>%
    write_delim("genotypes.ped", col_names=FALSE)

  """

}

process get_map {

  publishDir "../../data/little_green_men", overwrite: true

  input:
    file genome_map
  output:
    file "genotypes.map" into map
    file "gene2snp.tsv" into gene2snp

  """
  #!/usr/bin/env Rscript
  library(readr)
  library(magrittr)
  library(dplyr)

  load("$genome_map")
  step <- 5000

  df <- genome %>%
    arrange(desc(causalGene), chr, id)

  map <- df[rep(row.names(df), df\$snps), ] %>%
    select(chr, id) %>%
    mutate(name = paste0("rs", 1:n()),
           geneticDistance = 0) %>%
    group_by(chr) %>%
    mutate(position = step + 1:n() * step)

  map %>%
    select(chr, name, geneticDistance, position) %>%
    write.table("genotypes.map", sep = "\\t", row.names = F, col.names = F, quote = F)

  map %>%
    select(name, id) %>%
    rename(SNP = name, GENE = id) %>%
    write_tsv("gene2snp.tsv")

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

process get_ppi {

  publishDir "../../data/little_green_men", overwrite: true

  input:
    file genome_ppi
  output:
    file "ppi.tab" into ppi

  """
  #!/usr/bin/env Rscript
  library(readr)
  library(magrittr)
  library(dplyr)

  load("$genome_ppi")

  genome %>%
    apply(1, function(gene){
      lapply(gene[[5]], function(intx, id){
        data.frame(OFFICIAL_SYMBOL_FOR_A = id, OFFICIAL_SYMBOL_FOR_B = as.integer(intx))
        }, gene[[1]]) %>% do.call("rbind", .)
      }) %>% do.call("rbind", .) %>%
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

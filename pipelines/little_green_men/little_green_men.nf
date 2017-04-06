#!/usr/bin/env nextflow

// params for genotype generation
params.N = 3000
params.rr = 30

// params for genome generation
params.numCausalGenes = 4
params.numGenesPerPathway = 6
params.numPathways = 30
params.numSnpsPerGene = 10
params.numChromosomes = 6
params.freqCrosstalk = 0.3
params.freqCausalSnp = 0.4

// params for system generation
params.out = "."

if (params.population == "custom"){

  process getCustomGenome {

    output:
      file "genome.Rdata" into genome_ped, genome_map, genome_ppi

    """
    #!/usr/bin/env Rscript
    genome <- data.frame(geneId = paste0("g", 1:5), chr = c(1,1,1,2,2),
                         causalGene = c(F,T,F,F,T),
                         snps = c(2,3,2,4,4))
    genome\$ppi <- list(c("g2"), c("g1","g5"), c("g4"), c("g3","g5"), c("g4","g2))

    save(genome, file = "genome.Rdata")
    """
  }
} else if (params.population == "random") {

  process generateGenesChromosomesAndPPIs {

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
    freqCrosstalk <- $params.freqCrosstalk

    # create data frame
    numGenes <- numPathways * numGenesPerPathway
    causalGenes <- paste0("g", 1:numCausalGenes)

    genome <- data.frame(geneId = paste0("g", 1:numGenes),
                         chr = sample(numChromosomes, numGenes, replace = TRUE)) %>%
      mutate(causalGene = geneId %in% causalGenes,
             snps = numSnpsPerGene)

    # calculate ppi
    # clique of causal genes (all combinations are possible)
    ppiCausal <- expand.grid(causalGenes, causalGenes, stringsAsFactors = FALSE) %>%
      set_colnames(c("Gene1", "Gene2"))

    # sequential
    ppiPathway.model <- data.frame( Gene1 = (1:numGenesPerPathway)[-numGenesPerPathway],
                                    Gene2 = (1:numGenesPerPathway)[-1])

    ppiPathway <- lapply(0:(numPathways - 1), function(i){
      ppiPathway.model + i * numGenesPerPathway
    }) %>% do.call("rbind", .) %>%
      set_colnames(c("Gene1", "Gene2")) %>%
      mutate(Gene1 = paste0("g", Gene1),
             Gene2 = paste0("g", Gene2))

    ppiCrosstalk <- lapply(1:numPathways, function(i){
      first <- 1 + numGenesPerPathway * (i - 1)
      last <- numGenesPerPathway * i

      pwGenes <- genome\$geneId[first:last]
      otherGenes <- genome\$geneId[-(first:last)]

      possibleCrosstalk <- expand.grid(pwGenes, otherGenes)

      sample_n(possibleCrosstalk, floor(numGenesPerPathway * freqCrosstalk))
    }) %>% do.call("rbind", .) %>%
      set_colnames(c("Gene1", "Gene2")) %>%
      mutate(Gene1 = as.character(Gene1), Gene2 = as.character(Gene2))

    ppi <- rbind(ppiCausal, ppiPathway, ppiCrosstalk)
    # make it symmetrical
    ppi <- data.frame(Gene1 = c(ppi\$Gene1, ppi\$Gene2),
                      Gene2 = c(ppi\$Gene2, ppi\$Gene1))

    genome\$ppi <- apply(genome, 1, function(x){
      ppi\$Gene2[ppi\$Gene1 == x[[1]]] %>% sort %>% unique
    })

    # save results
    save(genome, file = "genome.Rdata")
    """
  }

}

process getGenotypes{

  publishDir "$params.out", overwrite: true, mode: "copy"

  input:
    file genome_ped
  output:
    file "genotypes.ped" into ped

  """
  #!/usr/bin/env Rscript
  library(readr)
  library(magrittr)
  library(dplyr)

  sampleCausalSNP <- function(i, ref, causal, rr, N){
    L <- rr/(rr+1)

    gt.cases <- rep(ref, 2 * N)
    gt.cases[runif(2 * N) < L] <- causal

    gt.controls <- rep(ref, 2 * N)
    gt.controls[runif(2 * N) > L] <- causal

    data.frame(chr1 = c(gt.cases[1:N], gt.controls[1:N]),
               chr2 = c(gt.cases[(N + 1):(2 * N)], gt.controls[(N + 1):(2 * N)]))
  }

  sampleNeutralSNP <- function(i, ref, alt, maf, N){

    gt <- rep(ref, 4 * N)
    gt[runif(4 * N) < maf] <- alt

    data.frame(chr1 = gt[1:(2 * N)], chr2 = gt[(2 * N + 1):(4 * N)])
  }

  rr = $params.rr
  N = $params.N
  freqCausalSnp <- $params.freqCausalSnp
  maf = 0.25

  load("$genome_ped")

  # involved snps
  numCausalSnps <- genome %>% filter(causalGene) %>% .\$snps %>% sum(.) * freqCausalSnp
  numCausalSnps <- floor(numCausalSnps)

  # A = low risk allele
  # C = high risk allele
  causalSNPs <- lapply(1:numCausalSnps, sampleCausalSNP, "A", "C", rr, N) %>%
    do.call("cbind", .)

  # neutral snps
  # T,G = alleles without impact
  neutralSNPs <- lapply(1:(sum(genome\$snps) - numCausalSnps), sampleNeutralSNP, "T", "G", maf, N) %>%
    do.call("cbind", .)

  snps <- cbind(causalSNPs, neutralSNPs)
  ids <- 1:(2*N)
  data.frame(family = ids, indiv = ids,
             pater = as.integer(0), mater = as.integer(0),
             sex = as.integer(2),
             pheno = as.integer(c(rep(2, N), rep(1, N)))) %>%
    cbind(snps) %>%
    write_delim("genotypes.ped", col_names=FALSE)

  """

}

process getSNPs {

  publishDir "$params.out", overwrite: true, mode: "copy"

  input:
    file genome_map
  output:
    file "genotypes.map" into map
    file "gene2snp.tsv" into gene2snp
    file "truth.tsv" into truth

  """
  #!/usr/bin/env Rscript
  library(readr)
  library(magrittr)
  library(dplyr)

  # disable scientific notation
  options(scipen=999)

  load("$genome_map")
  step <- 5000

  causal <- genome %>%
    .[rep(row.names(.), .\$snps), ] %>%
    group_by(geneId) %>%
    mutate(causalSnp = causalGene * c(
                                    rep(1, floor($params.freqCausalSnp * n())),
                                    rep(0, n() - floor($params.freqCausalSnp * n()))),
           causalSnp = sample(causalSnp)) %>%
    ungroup %>%
    arrange(desc(causalSnp)) %>%
    mutate(rs = paste0("rs", 1:n())) %>%
    arrange(chr, geneId) %>%
    select(geneId, rs, causalSnp)

  pos <- genome %>%
    .[rep(row.names(.), .\$snps), ] %>%
    arrange(chr, geneId) %>%
    group_by(chr) %>%
    mutate(geneId = sample(geneId),
           position = step + 1:n() * step,
           geneticDistance = 0) %>%
    ungroup %>%
    select(chr, position, geneticDistance)

  snps <- cbind(pos, causal) %>%
    mutate(snpNum = gsub("rs", "", rs) %>% as.numeric) %>%
    arrange(snpNum)

  snps %>%
    write_tsv("truth.tsv")

  snps %>%
    select(chr, rs, geneticDistance, position) %>%
    write.table("genotypes.map", sep = "\\t", row.names = F, col.names = F, quote = F)

  snps %>%
    select(rs, geneId) %>%
    unique %>%
    rename(SNP = rs, GENE = geneId) %>%
    write_tsv("gene2snp.tsv")

  """

}

process getTAB {

  publishDir "$params.out", overwrite: true, mode: "copy"

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
      lapply(gene[[5]], function(intx, geneId){
        data.frame(OFFICIAL_SYMBOL_FOR_A = geneId, OFFICIAL_SYMBOL_FOR_B = intx)
        }, gene[[1]]) %>% do.call("rbind", .)
      }) %>% do.call("rbind", .) %>%
    mutate(INTERACTOR_A = "", INTERACTOR_B = "", ALIASES_FOR_A = "",
           ALIASES_FOR_B = "", EXPERIMENTAL_SYSTEM = "", SOURCE = "",
           PUBMED_ID = "", ORGANISM_A_ID = "", ORGANISM_B_ID = "") %>%
    write_tsv("ppi.tab")

  """

}

process getSnp_list {

  publishDir "$params.out", overwrite: true, mode: "copy"

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

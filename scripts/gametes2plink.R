#!/usr/bin/env Rscript
source("~/wisdom/r/data_analysis_environment.R")

setwd(paste0("~/genewa"))
wd <- "populations/gametes/"

args <- commandArgs(trailingOnly = TRUE)
replicate <- paste0(wd,args[1])

# generate MAP file
snpNames <- c(paste0("N",0:97),c("M0P0","M0P1"))

data_frame(chromosome = 0, 
           snp = snpNames,
           morgans = 0,
           pos = 0) %>%
  write_tsv(paste0(wd, "map.txt"), col_names = FALSE)

# generate PED file
genotypes <- read_tsv(replicate) %>%
  ## create individual identifiers
  mutate(family_id = 1:2000, 
         individual_id = 1:2000, 
         paternal_id = 1:2000, 
         maternal_id = 1:2000,
         sex = "unknown",
         phenotype = Class + 1) %>%
  select(-Class) %>%
  select(family_id:phenotype, N0:M0P1)

## convert from 0,1,2 to phased phenotypes
snps <- genotypes %>% select(N0:M0P1) %>% as.matrix
snps[snps == 0] <- "1_1"
snps[snps == 1] <- "1_2"
snps[snps == 2] <- "2_2"

## split into columns and create the matrix
snps.chr <- as.vector(snps) %>% data_frame(genotypes = .) %>% separate(genotypes, c("chr1","chr2"), "_") %>% as.matrix
dim(snps.chr) <- c(dim(snps)[1],2 * dim(snps)[2])

genotypes %>%
  select(family_id:phenotype) %>%
  cbind(snps.chr) %>%
  write_tsv(paste0(replicate,".ped"), col_names = FALSE)

system2("plink", 
        c("--noweb",
          "--epistasis",
          "--ped", paste0(replicate,".ped"),
          "--map", paste0(wd, "map.txt"),
          "--out", paste0(replicate, ".plink"),
          "--allow-no-sex",
          "--epi1 1", "--epi2 1"))
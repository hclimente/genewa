#!/usr/bin/env Rscript
library(tidyr)
library(magrittr)
library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
n <- as.integer(args[2])
p <- as.integer(args[3])
map <- args[4]
ped <- args[5]

# generate MAP file
snpNames <- c(paste0("N",0:(p-3)),c("M0P0","M0P1"))

if (!file.exists(map)){
  data_frame(chromosome = 0,
             snp = snpNames,
             morgans = 0,
             pos = 0) %>%
    write.table(map, col.names = F, sep = "\t", quote = F, row.names = F)
}

# generate PED file
genotypes <- read_tsv(sample) %>%
  ## create individual identifiers
  mutate(family_id = 1:(2*n),
         individual_id = 1:(2*n),
         paternal_id = 1:(2*n),
         maternal_id = 1:(2*n),
         sex = "unknown",
         phenotype = as.integer(Class + 1)) %>%
  select(-Class) %>%
  select(family_id:phenotype, N0:M0P1)

## convert from 0,1,2 to phased genotypes
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
  write_tsv(ped, col_names = FALSE)

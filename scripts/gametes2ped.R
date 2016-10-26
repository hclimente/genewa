#!/usr/bin/env Rscript
source("~/libs/wisdom/r/data_analysis_environment.R")

setwd("~/genewa/populations/gametes")

args <- commandArgs(trailingOnly = TRUE)
h <- args[1]
maf <- args[2]
N <- args[3]
modelNo <- args[4]
repNo <- args[5]
  
# generate MAP file
snpNames <- c(paste0("N",0:(as.integer(N)-3)),c("M0P0","M0P1"))

if (!file.exists(paste0("map_",N,".txt"))){
  data_frame(chromosome = 0, 
             snp = snpNames,
             morgans = 0,
             pos = 0) %>%
    #write_tsv("map.txt", col_names = FALSE)
    write.table(paste0("map_",N,".txt"), col.names = F, sep = "\t", quote = F, row.names = F)
}

model <- paste0("h",h,"_maf",maf,"_N",N,"_EDM-",modelNo)
replicate <- paste0(model,"_",repNo,".txt")
          
# generate PED file
genotypes <- read_tsv(paste(model,replicate, sep = "/")) %>%
  ## create individual identifiers
  mutate(family_id = 1:2000, 
         individual_id = 1:2000, 
         paternal_id = 1:2000, 
         maternal_id = 1:2000,
         sex = "unknown",
         phenotype = Class + 1) %>%
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
  write_tsv(paste(model,paste0(replicate,".ped"), sep = "/"), col_names = FALSE)
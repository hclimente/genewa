#!/usr/bin/env Rscript
source("~/libs/wisdom/r/data_analysis_environment.R")

setwd("~/genewa/populations/gametes")

# generate MAP file
snpNames <- c(paste0("N",0:97),c("M0P0","M0P1"))

data_frame(chromosome = 0, 
           snp = snpNames,
           morgans = 0,
           pos = 0) %>%
  write_tsv("map.txt", col_names = FALSE)

for (h in c("005","01","025","05","1","2")){
  for (maf in c("2","4")){
    for (model in formatC(1:10, width = 2, flag = "0")){
      model <- paste0("h",h,"_maf",maf,"_EDM-",model)
      for (repNo in formatC(1:100, width = 3, flag = "0")){
        replicate <- paste0(model,"_",repNo,".txt")
        
        # generate PED file
        genotypes <- read_tsv(paste(model,replicate, sep = "/")) %>%
          ## create individual identifiers
          mutate(family_id = 1:2000, 
                 individual_id = 1:2000, 
                 paternal_id = 1:2000, 
                 maternal_id = 1:2000,
                 sex = "unknown",
                 phenotype = Class) %>%
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
        
      }
    }
  }
}
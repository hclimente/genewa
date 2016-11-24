#!/usr/bin/env Rscript
library(magrittr)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
ped <- args[1]
beam <- args[2]

snps <- read_tsv(ped, col_names = F) %>%
  # remove columns related to family, individual, father, mother and sex
  select(-X1, -X2, -X3, -X4, -X5) %>%
  # cases are 1 and controls are 0
  as.matrix
snps <- snps - 1L

# fetch phenotype information and remove that column
phenotypes <- snps[,1] %>% rep(each=2)
snps <- snps[,-1]

# take separately the two sets of chromosomes
haplo1 <- snps[,c(T,F)] %>% t %>% as.data.frame %>% set_colnames(., paste0("h1",1:dim(.)[2]))
haplo2 <- snps[,c(F,T)] %>% t %>% as.data.frame %>% set_colnames(., paste0("h2",1:dim(.)[2]))

# combine the data to create beam format
intertwinedCols <- paste0(c("h1","h2"),rep(1:dim(haplo2)[2], each=2))
cbind(haplo1, haplo2) %>%
  .[,intertwinedCols] %>%
  rbind(phenotypes, .) %>%
  write.table(beam, col.names = F, sep = " ", quote = F, row.names = F)

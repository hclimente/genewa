#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(tidyverse);
ensembl2hgnc <- read_tsv('../preprocessing/non_alt_loci_set.txt') %>%
    select(symbol, ensembl_gene_id)
vegas <- read_tsv(args[1], col_types = 'iciddddddcd') %>%
    left_join(ensembl2hgnc, by = c('Gene' = 'ensembl_gene_id')) %>%
    mutate(Gene = ifelse(!is.na(symbol), symbol, Gene)) %>%
    select(-symbol) %>%
    write_tsv(args[1])

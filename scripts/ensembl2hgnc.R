#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(tidyverse);
ensembl2hgnc <- read_tsv(args[2], col_types = cols(.default = "c")) %>%
    select(symbol, ensembl_gene_id)
vegas <- read_tsv(args[1], col_types = 'iciddddddcd') %>%
    inner_join(ensembl2hgnc, by = c('Gene' = 'ensembl_gene_id')) %>%
    mutate(Gene = symbol) %>%
    select(-symbol) %>%
    write_tsv(args[1])

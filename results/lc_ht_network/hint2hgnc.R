#!/usr/bin/env Rscript

library(tidyverse)

hgnc <- read_tsv('../preprocessing/non_alt_loci_set.txt') %>%
  separate_rows(alias_symbol, sep = '\\|') %>%
  separate_rows(prev_symbol, sep = '\\|') %>%
  select(symbol, alias_symbol, prev_symbol)

symbols <- select(hgnc, symbol) %>% mutate(matched = TRUE) %>% unique
aliases <- select(hgnc, symbol, alias_symbol) %>% unique
previous <- select(hgnc, symbol, prev_symbol) %>% unique

hint <- bind_rows(read_tsv('HomoSapiens_binary_hq.txt', col_types = 'ccccccccc'),
                  read_tsv('HomoSapiens_cocomp_hq.txt', col_types = 'ccccccccc')) %>%
  mutate(Gene_A = gsub('ORF', 'orf', Gene_A),
         Gene_B = gsub('ORF', 'orf', Gene_B),
         id = 1:n())

correct_gene_A <- left_join(hint, symbols, by = c('Gene_A' = 'symbol')) %>%
  left_join(aliases, by = c('Gene_A' = 'alias_symbol')) %>%
  left_join(previous, by = c('Gene_A' = 'prev_symbol'), suffix = c('_alias', '_prev')) %>%
  mutate(`Official Symbol Interactor A` = Gene_A,
         `Official Symbol Interactor A` = ifelse(is.na(matched) & !is.na(symbol_alias), symbol_alias, `Official Symbol Interactor A`),
         `Official Symbol Interactor A` = ifelse(is.na(matched) & is.na(symbol_alias) & !is.na(symbol_prev), symbol_prev, `Official Symbol Interactor A`)) %>%
  filter(!(is.na(matched) & is.na(symbol_alias) & is.na(symbol_prev))) %>%
  select(-matched, -symbol_alias, -symbol_prev) %>%
  unique

correct_genes <- left_join(correct_gene_A, symbols, by = c('Gene_B' = 'symbol')) %>%
  left_join(aliases, by = c('Gene_B' = 'alias_symbol')) %>%
  left_join(previous, by = c('Gene_B' = 'prev_symbol'), suffix = c('_alias', '_prev')) %>%
  mutate(`Official Symbol Interactor B` = Gene_B,
         `Official Symbol Interactor B` = ifelse(is.na(matched) & !is.na(symbol_alias), symbol_alias, `Official Symbol Interactor B`),
         `Official Symbol Interactor B` = ifelse(is.na(matched) & is.na(symbol_alias) & !is.na(symbol_prev), symbol_prev, `Official Symbol Interactor B`)) %>%
  filter(!(is.na(matched) & is.na(symbol_alias) & is.na(symbol_prev))) %>%
  select(-matched, -symbol_alias, -symbol_prev) %>%
  unique

write_tsv(correct_genes, 'hint.hgnc.pseudo.tab2')

#!/usr/bin/env Rscript
library(tidyverse)

snp2gene <- read_tsv('../preprocessing/snp2hgnc.tsv', col_types = 'cc')

snp_assoc <- read_tsv('../conventional_gwas/univariate_models.no_covars.tsv', 
                      col_types = 'icdccccddd') %>%
    rename(Chr = CHR) %>%
    select(-A1, -F_A, -F_U, -A2)

right_join(snp2gene, snp_assoc, by = c('snp' = 'SNP')) %>%
    group_by(gene) %>%
    top_frac(0.1, -P) %>%
    select(snp, gene) %>%
    write_tsv('snp2hgnc.top10.tsv')

#!/usr/bin/env Rscript

library(tidyverse)

# from ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/non_alt_loci_set.txt
symbols <- read_tsv('non_alt_loci_set.txt') %>%
	select(entrez_id, symbol)

biogrid <- read_tsv("BIOGRID-MV-Physical-3.5.171.tab2.txt") %>%
	filter(`Organism Interactor A` == 9606 & `Organism Interactor B` == 9606) %>%
	inner_join(symbols, by = c('Entrez Gene Interactor A' = 'entrez_id')) %>%
	inner_join(symbols, by = c('Entrez Gene Interactor B' = 'entrez_id')) %>%
	mutate(`Official Symbol Interactor A` = symbol.x, 
		   `Official Symbol Interactor B` = symbol.y) %>%
	select(-symbol.x, -symbol.y) %>%
	write_tsv('BIOGRID-MV-Physical-3.5.168.tab2.hgnc.tsv')

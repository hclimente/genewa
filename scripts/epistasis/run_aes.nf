ped = file("$params.ped")
map = file("$params.map")

process ped2beam {

	input:
		file ped
		file map

	output:
		file "${ped.baseName}.aes" into aesIn

	"""
	#!/usr/bin/env Rscript

	library(tidyverse)
	library(dplyr)
	library(magrittr)

	# read files
	map <- read_tsv("$map", col_names = F) %>%
	  set_colnames(c("chr","snp","cm","pos")) %>%
	  select(snp)
	ped <- read_delim("$ped", " ", col_names = F)

	class <- as.integer(ped\$X6 - 1)

	# artificially phase the genotypes
	chr1 <- ped %>%
	  select(-num_range("X", 1:6)) %>%
	  .[,c(T,F)] %>%
	  set_colnames(map\$snp)

	chr2 <- ped %>%
	  select(-num_range("X", 1:6)) %>%
	  .[,c(F,T)] %>%
	  set_colnames(map\$snp)

	# get mode
	mode <- rbind(chr1, chr2) %>%
	  apply(2, function(x) {
	    ux <- unique(x)
	    ux[which.max(tabulate(match(x, ux)))]})

	# get # different alleles per snp
	nalleles <- rbind(chr1, chr2) %>%
	  apply(2, function(x) {length(unique(x))})

	aes <- (chr1 == mode) + (chr2 == mode) %>%
	  as.data.frame %>%
	  # remove those snps with more than 2 alleles (cannot be encoded as 0,1,2)
	  subset(select=(nalleles <= 2))
	aes\$class = class

	write_csv(aes, "${ped.baseName}.aes")
	"""

}

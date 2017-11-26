#!/usr/bin/env nextflow

rdata = file("$params.rdata")
map = file("$params.map")
params.out = "."

p = Channel
        .fromPath("$params.map")
        .splitText()
        .count()

process r2boost {

	input:
		file rdata

	output:
		file "genotypes.boost" into boostIn

	"""
	#!/usr/bin/env Rscript

	library(snpStats)
	library(tidyverse)
	load("$rdata")

	X <- as(gwas\$genotypes, "numeric")
	Y <- gwas\$fam\$affected - 1

	cbind(Y,X) %>% as.data.frame %>% write_delim("genotypes.boost", col_names = F)
	"""

}

process runBOOST {

	input:
		file boostIn

	output:
		file "tempInteractionRecords.txt" into interactions

	"""
	echo $boostIn >listFile
	GBOOST -i listFile -wm GPU
	"""

}

process calculatePValues {

    publishDir "$params.out", overwrite: true, mode: "copy"

	input:
		file map
		file interactions
		val p

	output:
		file "${ped.baseName}.boost.txt" into epistasis

	"""
	#!/usr/bin/env Rscript

	library(tidyverse)
	library(magrittr)

	map <- read_tsv("$map", col_names = F) %>%
		set_colnames(c("chr", "snp", "pos", "gpos"))

	read_tsv("$interactions", col_names = F) %>%
		set_colnames(c("index","SNP1","SNP2","singlelocusAssoc1","singlelocusAssoc2","InteractionBOOST","InteractionPLINK")) %>%
		mutate(SNP1 = map\$snp[SNP1],
		       SNP2 = map\$snp[SNP2],
		       p_singlelocusAssoc1 = pchisq(singlelocusAssoc1, 2, lower.tail = F),
		       p_singlelocusAssoc2 = pchisq(singlelocusAssoc2, 2, lower.tail = F),
		       p_InteractionBOOST = pchisq(InteractionBOOST, 4, lower.tail = F),
		       padj_singlelocusAssoc1 = p.adjust(p_singlelocusAssoc1, n = $p),
                       padj_singlelocusAssoc2 = p.adjust(p_singlelocusAssoc2, n = $p),
                       padj_InteractionBOOST = p.adjust(p_InteractionBOOST, n = $p*($p - 1)/2)) %>%
		write_tsv("${ped.baseName}.boost.txt")
	"""

}

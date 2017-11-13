#!/usr/bin/env nextflow

ped = file("$params.ped")
map = file("$params.map")
params.out = "."

srcReadPed = file("$genewawd/martiniflow/io/read_ped.nf")

process ped2boost {

	input:
		file ped
		file map

	output:
		file "genotypes.boost" into boostIn

	beforeScript: "nextflow run $srcReadPed --ped $ped --map $map"

	script:
	"""
	#!/usr/bin/env Rscript

	library(snpStats)
	library(tidyverse)
	load("gwas.RData")

	X <- as(gwas\$genotypes, "numeric")
	Y <- ifelse(gwas\$fam\$affected == 2, 1, 2)

	cbind(Y,X) %>% write_delim("genotypes.boost")
	"""

}

process runBOOST {

	input:
		file boostIn

	"""
	echo $boostIn >listFile
	GBOOST -i listFile -wm GPU -o gboost.
	"""

}

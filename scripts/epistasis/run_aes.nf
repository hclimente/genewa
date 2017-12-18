rdata = file("$params.rdata")
map = file("$params.map")
params.out = "."

p = Channel
        .fromPath("$params.map")
        .splitText()
        .count()

process r2aes {

	input:
		file rdata

	output:
		file "genotypes.aes" into aesIn

	"""
	#!/usr/bin/env Rscript

	library(tidyverse)
	library(snpStats)
	load("$rdata")

	X <- as(gwas\$genotypes, "numeric")
	Y <- gwas\$fam\$affected - 1

	cbind(X,Y) %>% as.data.frame %>% write_csv("genotypes.aes")
	"""

}

process makeParameters {

	input:
		val p

	output:
		file "parameters.txt" into parameters

	"""
	cat << EOF >parameters.txt
	iAntCount	1000    		//number of ants
	iItCountLarge	${p * 0.1}  //number of iterations for the large haplotypes
	iItCountSmall	${p * 0.05}	//number of iterations for the small haplotypes
	alpha	1               	//weight given to pheromone deposited by ants
	iTopModel	1000     		//number of top ranking haplotypes in the first stage
	iTopLoci	200      		//number of loci with top ranking pheromone in the first stage
	rou	0.05       				//evaporation rate in Ant Colony Optimizaion
	phe	100        				//initial pheromone level for each locus
	largehapsize	6 			//size of the large haplotypes
	smallhapsize	3 			//size of the small haplotypes
	iEpiModel	2     			//number of SNPs in an epistatic interaction
	pvalue	0.01   				//p value threshold (after Bonferroni correction)
	INPFILE	genotypes.aes		//input file name for case-control genotype data
	OUTFILE result.txt 			//output file name for detected epistatic interactions
	EOF
	"""

}

process runAntEpiSeeker {

	validExitStatus 0,1

	input:
		file aesIn
		file parameters

	output:
		file "result.txt" into aesOut

	"""
	AntEpiSeeker
	"""

}

process processResults {

	publishDir "$params.out", overwrite: true, mode: "copy"

	input:
		file aesOut

	output:
		file "${map.baseName}.aes.txt" into epistasis

	"""
	echo -e "snp1\\tsnp2\\tchi2\\tp" >${map.baseName}.aes.txt
	tail -n+3 $aesOut | sed -E 's/[0-9]+\\(//g' | sed 's/)//g' | sed -E 's/[ \\t]+/\\t/g' >>${map.baseName}.aes.txt
	"""

}

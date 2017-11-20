rdata = file("$params.rdata")
params.out = "."

p = Channel
	.fromPath("$params.map")
	.splitText()
	.count()

process r2beam {

	input:
		file rdata

	output:
		file "genotypes.beam" into beamIn

	"""
	#!/usr/bin/env Rscript
	library(tidyverse)
	library(snpStats)
	load("$rdata")

	X <- as(gwas\$genotypes, "numeric")
	Y <- gwas\$fam\$affected - 1
	snps <- gwas\$map[,c(2,1,4)]
	snps[,1] <- as.character(snps[,1])
	snps <- rbind(c("ID","Chr", "Pos"), snps)

	cbind(Y, X) %>%
		as.data.frame %>%
		arrange(-Y) %>%
		t %>%
		as.data.frame %>%
		cbind(snps, .) %>%
		write_delim("genotypes.beam", col_names = F)
	"""

}

process get_parameters {

	input:
		val p

	output:
		file "parameters.txt" into beamParams

	"""
	cat << EOF >parameters.txt
	[RUNNING PARAMETER]
	CHAIN				100		// number of MCMC chains
	BURNIN				0		// number of burnin updates
	MCMC				0		// number of updates after burnin
	THIN				$p		// number of updates between two posterior samples; Recommended: = number of markers
	PRIOR1				0.01	// prior probability for each marker to belong to group 1 (marginal group)
	PRIOR2				0.01	// prior probability for each marker to belong to group 2 (interacting group)
	SINGLE_ONLY			0		// 1: test for marginal associations only; 0: test for both marginal and interaction associatoin
	MINDISTANCE			0		// let BEAM to detect interaction between SNPs that are at least MINDISTANCE bps apart
	INITIALTRYS			20		// if set to 0, the chains will run from random starting points, with annealing;
								// if set to a positive integer (x), we first search for x good starting points, and then run chains from these points, without annealing.
	TRY_LENGTH 			0		// this is a multiplier (x) that determines how long each trial will run. For L markers, each trial runs for xL iterations.
	AUTORESTART			5		// during MCMC, whenever max(LogP) is larger than initial(LogP) + AUTORESTART, then MCMC restarts from this "better" mode.
								// this should be a positive number. to disable autorestart, set this value high (e.g., 1000000)

	[INPUT FORMAT]
	INPFILE				genotypes.beam	// input file name for case-control data
	INC_SNP_ID			1				// input file includes SNP ID, e.g., rs329040
	INC_SNP_POS			1				// input file includes SNP locations (in bytes), e.g., Chr10  1042329

	[OUTPUT FORMAT]
	OUTFILE				${ped.baseName}.beam.txt 	// output file name for posterior distributions and detected associations
	P_THRESHOLD 		1
	EOF
	"""

}

process run_beam {

	publishDir "$params.out", overwrite: true, mode: "copy"

	input:
		file beamIn
		file beamParams

	output:
		file "${ped.baseName}.beam.txt" into beamOut

	"""
	BEAM
	"""

}

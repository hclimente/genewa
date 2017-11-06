ped = file("$params.ped")
map = file("$params.map")

p = Channel
	.fromPath("$params.map")
	.splitText()
	.count()

process ped2beam {

	input:
		file ped
		file map

	output:
		file "${ped.baseName}.beam" into beamIn

	"""
	#!/usr/bin/env Rscript
	library(tidyverse)
	library(snpStats)

	gwas <- read.pedfile("$ped", snps = "$map")
	X <- as(gwas\$genotypes, "numeric") %>% t
	X1 <- matrix(0, nrow(X), ncol(X))
	colnames(X1) <- paste0("h1",1:ncol(X1))
	X2 <- X1
	colnames(X2) <- paste0("h2",1:ncol(X2))

	X1[X > 0] <- 1
	X2[X == 2] <- 1
	
	phenotypes <- (gwas\$fam\$affected - 1)  %>% rep(each=2)

	# combine the data to create beam format
	intertwinedCols <- paste0(c("h1","h2"),rep(1:ncol(X1), each=2))
	cbind(X1, X2) %>%
	  .[,intertwinedCols] %>%
	  rbind(phenotypes, .) %>%
	  write.table("${ped.baseName}.beam", col.names = F, sep = " ", quote = F, row.names = F)
	"""

}

process get_parameters {

	output:
		file "parameters.txt" into beamParams

	"""
	cat << EOF >parameters.txt
	[RUNNING PARAMETER]
	CHAIN					100						// number of MCMC chains
	BURNIN				0							// number of burnin updates
	MCMC					0							// number of updates after burnin
	THIN					$p						// number of updates between two posterior samples; Recommended: = number of markers
	PRIOR1				0.01					// prior probability for each marker to belong to group 1 (marginal group)
	PRIOR2				0.01					// prior probability for each marker to belong to group 2 (interacting group)
	SINGLE_ONLY		0							// 1: test for marginal associations only; 0: test for both marginal and interaction associatoin
	MINDISTANCE		0							// let BEAM to detect interaction between SNPs that are at least MINDISTANCE bps apart
	INITIALTRYS		20						// if set to 0, the chains will run from random starting points, with annealing;
															// if set to a positive integer (x), we first search for x good starting points, and then run chains from these points, without annealing.
	TRY_LENGTH		0							// this is a multiplier (x) that determines how long each trial will run. For L markers, each trial runs for xL iterations.
	AUTORESTART		5							// during MCMC, whenever max(LogP) is larger than initial(LogP) + AUTORESTART, then MCMC restarts from this "better" mode.
															// this should be a positive number. to disable autorestart, set this value high (e.g., 1000000)

	[INPUT FORMAT]
	INPFILE				${ped.baseName}.beam				// input file name for case-control data
	INC_SNP_ID		0							// input file includes SNP ID, e.g., rs329040
	INC_SNP_POS		0							// input file includes SNP locations (in bytes), e.g., Chr10  1042329

	[OUTPUT FORMAT]
	OUTFILE				${ped.baseName}.beam.txt	// output file name for posterior distributions and detected associations
	P_THRESHOLD 1
	EOF
	"""

}

process run_beam {

	input:
		file beamIn
		file beamParams

	output:
		file "${ped.baseName}.beam.txt" into beamOut

	"""
	BEAM $beamIn
	"""

}

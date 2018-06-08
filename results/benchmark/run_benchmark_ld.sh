#!/bin/bash

benchmark_from_experiment.nf --geno data/genesis/genesis.processed --rld data/genesis/ld_snpStats.RData --snp2gene data/genesis/snp2hgnc.tsv --tab2 data/genesis/BIOGRID-ORGANISM-Homo_sapiens-3.4.160.tab2.hgnc.tsv --cases 636 --controls 637 --permutations 10  --genewawd /data/users/hcliment/projects/genewa --ngenes 20 --nets gm --psnps 0.5 --heritabilities 1 -resume

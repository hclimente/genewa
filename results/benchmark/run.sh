#!/bin/bash

benchmark_from_experiment.nf --geno data/genesis/genesis.processed --rld data/genesis/ld_snpStats.RData --snp2gene data/genesis/snp2hgnc.tsv --tab data/genesis/BIOGRID-ORGANISM-Homo_sapiens-3.4.160.tab.txt --cases 636 --controls 637 --permutations 10  --genewawd /data/users/hcliment/projects/genewa --ngenes 20 -resume

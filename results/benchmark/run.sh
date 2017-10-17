#!/bin/bash

benchmark_from_experiment.nf --geno data/genesis/genesis.processed --snp2gene data/genesis/gene2snp.hg19 --tab data/genesis/BIOGRID-ORGANISM-Homo_sapiens-3.4.138.tab.txt --cases 636 --controls 637 --permutations 10 --rld data/genesis/ld.RData --genewawd /data/users/hcliment/projects/genewa --ngenes 20 -resume

#!/bin/bash

benchmark_from_experiment.nf --geno datasets/genesis/genesis.processed --snp2gene datasets/genesis/gene2snp.hg19 --tab datasets/genesis/BIOGRID-ORGANISM-Homo_sapiens-3.4.138.tab.txt --cases 636 --controls 637 --permutations 10 --rld datasets/genesis/ld.RData --genewawd /data/users/hcliment/projects/genewa --ngenes 20 -resume

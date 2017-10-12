#!/bin/bash

benchmark_from_experiment.nf --geno datasets/lgm/genotypes --snp2gene datasets/lgm/gene2snp.tsv --tab datasets/lgm/ppi.tab --cases 1500 --controls 1500 --permutations 2 --rld ld.RData --ngenes 3 -resume

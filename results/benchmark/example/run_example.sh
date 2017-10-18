#!/bin/bash

benchmark_from_experiment.nf --geno data/lgm/genotypes --snp2gene data/lgm/gene2snp.tsv --tab data/lgm/ppi.tab --cases 1500 --controls 1500 --permutations 1 --rld data/lgm/ld.RData --ngenes 3 -resume

#!/bin/bash

generate_gwas_from_experiment.nf --geno datasets/lgm/genotypes --snp2gene datasets/lgm/gene2snp.tsv --tab datasets/lgm/ppi.tab --h2 1 --net gi -resume --n 3000 --permutations 2

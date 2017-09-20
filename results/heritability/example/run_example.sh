#!/bin/bash

benchmark_from_experiment.nf  --geno datasets/lgm/genotypes --snp2gene datasets/lgm/gene2snp.tsv --tab datasets/lgm/ppi.tab --n 3000 --permutations 2 --p 2 -resume

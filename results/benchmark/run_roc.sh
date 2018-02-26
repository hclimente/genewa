#!/bin/bash

roc_curve.nf --geno data/genesis/genesis.processed.pruned --rld data/genesis/ld_snpStats.RData --snp2gene data/genesis/gene2snp.hg19 --tab data/genesis/BIOGRID-ORGANISM-Homo_sapiens-3.4.138.tab.txt --cases 636 --controls 637 --permutations 10 --nets gs,gm,gi,gi2 --etas 0,1,2,3,-1 --lambdas 0,1,2,3,-1 --genewawd /data/users/hcliment/projects/genewa --ngenes 20 --psnps 0.5 -resume

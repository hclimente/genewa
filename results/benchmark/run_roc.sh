#!/bin/bash

roc_curve.nf --geno data/genesis/genesis.processed --rld data/genesis/ld_snpStats.RData --snp2gene data/genesis/gene2snp.hg19 --tab data/genesis/BIOGRID-ORGANISM-Homo_sapiens-3.4.138.tab.txt --cases 636 --controls 637 --permutations 10 --net gs,gm,gi,gi2 --etas 0,0.5,1,1.5,2,2.5,3 --lambdas 4,3.5,3,2.5,2,1.5,1,0.5,0,-0.5,-1 --genewawd /data/users/hcliment/projects/genewa --ngenes 20 --psnps 0.5 -resume

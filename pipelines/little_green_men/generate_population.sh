nextflow run ../../pipelines/little_green_men/little_green_men.nf --population random --rr $1
nextflow run ../../pipelines/scones/get_snp_networks.nf --snp2gene gene2snp.tsv
nextflow run ../../pipelines/scones/get_phenotypes.nf --ped genotypes.ped

vegas2.nf --bfile genesis_2019 --gencode 31 --genome 37 --buffer 50000 --vegas_params '-top 10' -resume -profile bigmem
../../scripts/ensembl2hgnc.R scored_genes.vegas.txt non_alt_loci_set.txt

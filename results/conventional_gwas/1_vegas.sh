# ld pruned data 
vegas2.nf --bfile genesis_2019_pruned --genome 37 --buffer 50000 --vegas_params '-top 10 -chr 23' -resume -profile bigmem
mv scored_genes.vegas.txt scored_genes.X.vegas.txt

vegas2.nf --bfile genesis_2019_pruned --genome 37 --buffer 50000 --vegas_params '-top 10' -resume -profile bigmem
tail -n +2 scored_genes.X.vegas.txt >>scored_genes.vegas.txt

rm scored_genes.X.vegas.txt
mv scored_genes.vegas.txt scored_genes.vegas.pruned.txt
../../scripts/ensembl2hgnc.R scored_genes.vegas.pruned.txt ../preprocessing/non_alt_loci_set.txt

# with covariate correction 
vegas2.nf --bfile genesis_2019 --genome 37 --buffer 50000 --vegas_params '-top 10 -chr 23' --covar CT_age_cens_tronq.cov -resume -profile bigmem
mv scored_genes.vegas.txt scored_genes.X.vegas.txt

vegas2.nf --bfile genesis_2019 --genome 37 --buffer 50000 --vegas_params '-top 10' --covar CT_age_cens_tronq.cov -resume -profile bigmem
tail -n +2 scored_genes.X.vegas.txt >>scored_genes.vegas.txt

rm scored_genes.X.vegas.txt
mv scored_genes.vegas.txt scored_genes.vegas.age_adjust.txt
../../scripts/ensembl2hgnc.R scored_genes.vegas.age_adjust.txt ../preprocessing/non_alt_loci_set.txt

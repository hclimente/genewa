plink --bfile genesis_2019 -indep-pairwise 50 5 0.75 -out pruned
plink --bfile genesis_2019 -extract pruned.prune.in -make-bed -out genesis_2019_pruned

# run vegas
run_vegas --bfile genesis_2019_pruned --genome GRCh37 --vegas_params '\-top 10 -upper 50000 -lower 50000 -chr 23' -resume -profile bigmem
mv scored_genes.vegas.txt scored_genes.X.vegas.txt

run_vegas --bfile genesis_2019_pruned --genome GRCh37 --vegas_params '\-top 10 -upper 50000 -lower 50000' -resume -profile bigmem
tail -n +2 scored_genes.X.vegas.txt >>scored_genes.vegas.txt

rm scored_genes.X.vegas.txt
mv scored_genes.vegas.txt scored_genes.vegas.pruned.txt

# run vegas
run_vegas --bfile genesis_2019 --genome GRCh37 --vegas_params '\-top 10 -upper 50000 -lower 50000 -chr 23' --covar CT_age_cens_tronq.cov -resume -profile bigmem
mv scored_genes.vegas.txt scored_genes.X.vegas.txt

run_vegas --bfile genesis_2019 --genome GRCh37 --vegas_params '\-top 10 -upper 50000 -lower 50000' --covar covar CT_age_cens_tronq.cov -resume -profile bigmem
tail -n +2 scored_genes.X.vegas.txt >>scored_genes.vegas.txt

rm scored_genes.X.vegas.txt
mv scored_genes.vegas.txt scored_genes.vegas.age_adjust.txt

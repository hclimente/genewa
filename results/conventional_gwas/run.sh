# SNP LEVEL STATISTICS
plink -bfile genesis_2019 --assoc --adjust
sed 's/ \+/\t/g' plink.assoc | sed 's/^\t//' | sed 's/\t$//' >univariate_models.no_covars.tsv
rm plink.assoc

plink -bfile genesis_2019 --logistic --covar CT_age_cens_tronq.cov
sed 's/ \+/\t/g' plink.assoc.logistic | sed 's/^\t//' | sed 's/\t$//' >univariate_models.covars.tsv
rm plink.assoc.logistic

plink -bfile genesis_2019 --logistic
sed 's/ \+/\t/g' plink.assoc.logistic | sed 's/^\t//' | sed 's/\t$//' >univariate_models.no_covars.logistic.tsv
rm plink.assoc.logistic

# GENE LEVEL STATISTICS
plink --bfile genesis_2019 -indep-pairwise 50 5 0.75 -out pruned
plink --bfile genesis_2019 -extract pruned.prune.in -make-bed -out genesis_2019_pruned

# run vegas
run_vegas --bfile genesis_2019_pruned --genome 37 --vegas_params '-top 10 -upper 50000 -lower 50000 -chr 23' -resume -profile bigmem
mv scored_genes.vegas.txt scored_genes.X.vegas.txt

run_vegas --bfile genesis_2019_pruned --genome 37 --vegas_params '-top 10 -upper 50000 -lower 50000' -resume -profile bigmem
tail -n +2 scored_genes.X.vegas.txt >>scored_genes.vegas.txt

rm scored_genes.X.vegas.txt
mv scored_genes.vegas.txt scored_genes.vegas.pruned.txt

# run vegas
run_vegas --bfile genesis_2019 --genome 37 --vegas_params '-top 10 -upper 50000 -lower 50000 -chr 23' --covar CT_age_cens_tronq.cov -resume -profile bigmem
mv scored_genes.vegas.txt scored_genes.X.vegas.txt

run_vegas --bfile genesis_2019 --genome 37 --vegas_params '-top 10 -upper 50000 -lower 50000' --covar CT_age_cens_tronq.cov -resume -profile bigmem
tail -n +2 scored_genes.X.vegas.txt >>scored_genes.vegas.txt

rm scored_genes.X.vegas.txt
mv scored_genes.vegas.txt scored_genes.vegas.age_adjust.txt

# ld prune 0.75
plink -bfile genesis_2019 -indep-pairwise 50 5 0.75 -out snps_75

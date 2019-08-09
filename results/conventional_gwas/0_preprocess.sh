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
plink -bfile genesis_2019 -indep-pairwise 50 5 0.75 -out snps_75
plink --bfile genesis_2019 -extract snps_75.prune.in -make-bed -out genesis_2019_pruned

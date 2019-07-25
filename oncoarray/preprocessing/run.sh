# compute pcs
plink -bfile dbgap_bcac_oncoarray_c8 --maf 0.05 -indep-pairwise 50 5 0.1 -out indep_snps
plink -bfile dbgap_bcac_oncoarray_c8 -extract indep_snps.prune.in -make-bed
plink -bfile plink --pca

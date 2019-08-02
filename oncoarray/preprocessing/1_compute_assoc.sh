cut -d' ' -f1-12 bcac_maf0.05_ld0.1.eigenvec >covars.txt
plink -bfile dbgap_bcac_oncoarray_c8 -logistic -covar covars.txt

sed 's/ \+/\t/g' plink.assoc.logistic | sed 's/^\t//' | sed 's/\t$//' >univariate.pc_corrected.tsv

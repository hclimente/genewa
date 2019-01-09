impute --ped Genesis.ped --map Genesis.map --reference 1000GP_Phase3 --strand_info strand_info --population EUR -resume -profile bigmem

# exclude 24 samples with bad genotyping
awk '{print $2,$2}' OFS='\t' Genesis.irem >excluded_samples

# exclude 11 cases (non-familiar, BRCA mutations)
awk '{print $2,$2}' OFS='\t' listNEW_IND_to_suppress.lst >>excluded_samples

# exclude samples and and 20 FGFR2 SNPs
plink -bfile out --remove excluded_samples --exclude FGFR2_SNPs_to_exclude.lst --make-bed --out filtered

# run vegas
run_vegas --bed filtered.bed --bim filtered.bim --fam filtered.fam --vegas_opts '\-top 10 -upper 50000 -lower 50000' -profile cluster -resume

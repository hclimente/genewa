# remove snps with unknown positions
awk '$4 == -1 {print $2}' genesis_raw.bim >unknown_pos
plink --bfile genesis_raw --exclude unknown_pos --make-bed --out genesis_raw.pos_ok

# impute
impute --bfile genesis_raw.pos_ok --genome GRCh37 --reference 1000GP_Phase3 --strand_info strand_info --population EUR -resume -profile bigmem

# exclude 11 cases (non-familiar, BRCA mutations)
awk '{print $2,$2}' OFS='\t' listNEW_IND_to_suppress.lst >excluded_samples

# exclude samples and and 20 FGFR2 SNPs and baly genotyped
plink -bfile out --remove excluded_samples --exclude FGFR2_SNPs_to_exclude.lst --geno 0.1 --make-bed --out filtered

# run vegas
run_vegas --bed filtered.bed --bim filtered.bim --fam filtered.fam --genome GRCh37 --vegas_params '\-top 10 -upper 50000 -lower 50000  -chr 23' --covar CT_age_cens_tronq.cov -resume -profile bigmem
mv scored_genes.vegas.txt scored_genes.X.vegas.txt

run_vegas --bed filtered.bed --bim filtered.bim --fam filtered.fam --genome GRCh37 --vegas_params '\-top 10 -upper 50000 -lower 50000' --covar CT_age_cens_tronq.cov -resume -profile bigmem
tail -n +2 scored_genes.X.vegas.txt >>scored_genes.vegas.txt

rm scored_genes.X.vegas.txt

# map snps to genes
snp2gene --bim filtered.bim --genome GRCh37 -profile cluster

# map BIOGRID to HGNC based on Entrez Ids
./biogrid2hgnc.R

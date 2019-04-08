# remove snps with unknown positions
awk '$4 == -1 {print $2}' genesis_raw.bim >unknown_pos
plink --bfile genesis_raw --exclude unknown_pos --make-bed --out genesis_raw_pos_ok

# impute
impute --bfile genesis_raw_pos_ok --genome GRCh37 --reference 1000GP_Phase3 --strand_info strand_info --population EUR -resume -profile bigmem

# exclude 11 cases
awk '{print $2,$2}' OFS='\t' listNEW_IND_to_suppress.lst >excluded_samples

# exclude samples and and 20 FGFR2 SNPs and baly genotyped
plink -bfile out --remove excluded_samples --exclude FGFR2_SNPs_to_exclude.lst --maf 0.001 --mind 0.10 --geno 0.1 --hwe 0.001 --make-bed --out tmp
plink -bfile tmp -indep-pairwise 50 5 0.75 -out pruned
plink -bfile tmp -extract pruned.prune.in -make-bed -out filtered

# run vegas
run_vegas --bfile filtered --genome GRCh37 --vegas_params '\-top 10 -upper 50000 -lower 50000  -chr 23' --covar CT_age_cens_tronq.cov -resume -profile bigmem
mv scored_genes.vegas.txt scored_genes.X.vegas.txt

run_vegas --bfile filtered --genome GRCh37 --vegas_params '\-top 10 -upper 50000 -lower 50000' --covar CT_age_cens_tronq.cov -resume -profile bigmem
tail -n +2 scored_genes.X.vegas.txt >>scored_genes.vegas.txt

rm scored_genes.X.vegas.txt

# map snps to genes
snp2gene --bim filtered.bim --genome GRCh37 -profile cluster

# map BIOGRID to HGNC based on Entrez Ids
./biogrid2hgnc.R

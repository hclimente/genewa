# remove snps with unknown positions
awk '$4 == -1 {print $2}' genesis_raw.bim >unknown_pos
plink --bfile genesis_raw --exclude unknown_pos --make-bed --out genesis_raw_pos_ok

# exclude 11 cases
awk '{print $2,$2}' OFS='\t' listNEW_IND_to_suppress.lst >excluded_samples

# exclude samples and and 20 FGFR2 SNPs and badly genotyped
plink -bfile genesis_raw_pos_ok --remove excluded_samples --exclude FGFR2_SNPs_to_exclude.lst --maf 0.001 --mind 0.10 --geno 0.1 --hwe 0.001 --make-bed --out genesis_2019 

# run vegas
run_vegas --bfile genesis_2019 --genome GRCh37 --vegas_params '\-top 10 -upper 50000 -lower 50000 -chr 23' -resume -profile bigmem
mv scored_genes.vegas.txt scored_genes.X.vegas.txt

run_vegas --bfile genesis_2019 --genome GRCh37 --vegas_params '\-top 10 -upper 50000 -lower 50000' -resume -profile bigmem
tail -n +2 scored_genes.X.vegas.txt >>scored_genes.vegas.txt

rm scored_genes.X.vegas.txt

# map snps to genes
snp2gene --bim genesis_2019.bim --genome GRCh37 -profile cluster

# map BIOGRID to HGNC based on Entrez Ids
./biogrid2hgnc.R

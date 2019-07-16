# remove snps with unknown positions
awk '$4 == -1 {print $2}' genesis_raw.bim >unknown_pos
plink --bfile genesis_raw --exclude unknown_pos --make-bed --out genesis_raw_pos_ok

# exclude 11 cases
awk '{print $2,$2}' OFS='\t' listNEW_IND_to_suppress.lst >excluded_samples

# exclude samples and and 20 FGFR2 SNPs and badly genotyped
plink -bfile genesis_raw_pos_ok --remove excluded_samples --exclude FGFR2_SNPs_to_exclude.lst --maf 0.001 --mind 0.10 --geno 0.1 --hwe 0.001 --make-bed --out genesis_2019

# BIOGRID to HGNC based on Entrez Ids
wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/non_alt_loci_set.txt
./hint2hgnc.R

# map snps to genes
snp2gene --bim genesis_2019.bim --genome GRCh37 --gencode_version 31 -profile cluster
sed 's/symbol/gene/' snp2hgnc.tsv >tmp

mv tmp snp2hgnc.tsv

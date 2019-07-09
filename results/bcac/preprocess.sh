# make plink output
echo -e 'CHR\tSNP\tBP\tA1\tF_A\tF_U\tA2\tCHISQ\tP\tOR' >icogs_bcac_public_results_euro.assoc
tail -n +2 ~/data/bcac/icogs_bcac_public_results_euro.txt | awk 'BEGIN { OFS = "\t" }{print $2, $6, $3, $4, $9, $9, $5, $14, $15, $13}' | grep -v NULL | sed 's/:[^\t]\+//' >>icogs_bcac_public_results_euro.assoc

# select only genesis snps
R -e 'library(tidyverse); snps <- read_tsv("~/data/genesis/genesis_2019.bim", col_names = FALSE) %>% rename(CHR = X1, SNP = X2, BP = X4) %>% select(CHR, SNP, BP); read_tsv("icogs_bcac_public_results_euro.assoc") %>% select(-SNP) %>% inner_join(snps) %>% select(CHR,SNP,BP,A1,F_A,F_U,A2,CHISQ,P,OR) %>% write_tsv("icogs_bcac_public_results_euro.genesis.assoc")'

cut -f2,9 icogs_bcac_public_results_euro.genesis.assoc | tail -n +2 >icogs_bcac_public_results_euro.genesis.assoc.vegas

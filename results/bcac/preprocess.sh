# make plink output
echo 'CHR\tSNP\tBP\tA1\tF_A\tF_U\tA2\tCHISQ\tP\tOR' >icogs_bcac_public_results_euro.assoc
tail -n +2 ~/data/bcac/icogs_bcac_public_results_euro.txt | awk 'BEGIN { OFS = "\t" }{print $2, $6, $3, $4, $9, $9, $5, $14, $15, $13}' | grep -v NULL | sed 's/:[^\t]\+//' >>icogs_bcac_public_results_euro.assoc

# select only genesis snps
R -e 'library(tidyverse); snps <- read_tsv("~/data/genesis/genesis_2019.bim", col_names = FALSE)$X2; read_tsv("icogs_bcac_public_results_euro.assoc") %>% filter(SNP %in% snps) %>% write_tsv("icogs_bcac_public_results_euro.genesis.assoc")'

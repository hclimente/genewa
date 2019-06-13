echo 'CHR\tSNP\tBP\tA1\tF_A\tF_U\tA2\tCHISQ\tP\tOR' >icogs_bcac_public_results_euro.assoc
tail -n +2 ~/data/bcac/icogs_bcac_public_results_euro.txt | awk 'BEGIN { OFS = "\t" }{print $2, $6, $3, $4, $9, $9, $5, $14, $15, $13}' | grep -v NULL | sed 's/:[^\t]\+//' >>icogs_bcac_public_results_euro.assoc

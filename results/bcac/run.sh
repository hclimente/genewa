# run vegas
run_vegas --snp_association icogs_bcac_public_results_euro.assoc --ld_controls g1000p3_EUR --genome GRCh37 --vegas_params '\-top 10 -upper 50000 -lower 50000 -chr 23' -resume -profile bigmem
mv scored_genes.vegas.txt scored_genes.X.vegas.txt

run_vegas --snp_association icogs_bcac_public_results_euro.assoc --ld_controls g1000p3_EUR --genome GRCh37 --vegas_params '\-top 10 -upper 50000 -lower 50000' -resume -profile bigmem
tail -n +2 scored_genes.X.vegas.txt >>scored_genes.vegas.txt

rm scored_genes.X.vegas.txt

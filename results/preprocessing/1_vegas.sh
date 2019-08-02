# run vegas
run_vegas --bfile genesis_2019 --gencode 31 --genome 37 --vegas_params '-top 10 -upper 50000 -lower 50000 -chr 23' -resume -profile bigmem
mv scored_genes.vegas.txt scored_genes.X.vegas.txt

run_vegas --bfile genesis_2019 --gencode 31 --genome 37 --vegas_params '-top 10 -upper 50000 -lower 50000' -resume -profile bigmem
tail -n +2 scored_genes.X.vegas.txt >>scored_genes.vegas.txt

rm scored_genes.X.vegas.txt

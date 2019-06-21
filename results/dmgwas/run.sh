cut -f2,9 ../preprocessing/scored_genes.vegas.txt | sed 's/Top-0.1-pvalue/Pvalue/' >scored_genes.top10.txt
run_dmGWAS --vegas scored_genes.top10.txt --tab2 ../preprocessing/BIOGRID-ORGANISM-Homo_sapiens-3.5.172.tab2.hgnc.txt -profile bigmem -resume

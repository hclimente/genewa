cut -f2,9 ../preprocessing/scored_genes.vegas.txt | sed 's/Top-0.1-pvalue/Pvalue/' >scored_genes.top10.txt
lean.nf --vegas scored_genes.top10.txt --tab2 ../preprocessing/hint.ht_complex.hgnc.pseudo.tab2 -resume -profile cluster

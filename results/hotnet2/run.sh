#git clone git@github.com:raphael-group/hotnet2.git
#cut -f2,9 ../preprocessing/scored_genes.vegas.txt | sed 's/Top-0.1-pvalue/Pvalue/' >scored_genes.top10.txt
hotnet2.nf --scores scored_genes.top10.txt --tab2 ../preprocessing/hint.ht_complex.hgnc.pseudo.tab2 --hotnet2_path hotnet2 -profile cluster -resume

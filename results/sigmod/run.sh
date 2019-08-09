#wget https://github.com/YuanlongLiu/SigMod/raw/master/SigMod_v2.zip && unzip SigMod_v2.zip
cut -f2,9 ../preprocessing/scored_genes.vegas.txt | sed 's/Top-0.1-pvalue/Pvalue/' >scored_genes.top10.txt
sigmod.nf --sigmod SigMod_v2 --vegas scored_genes.top10.txt --tab2 ../preprocessing/hint.ht_complex.hgnc.pseudo.tab2 -profile bigmem -resume 

#wget https://github.com/YuanlongLiu/SigMod/raw/master/SigMod_v2.zip && unzip SigMod_v2.zip
cut -f2,9 ../preprocessing/scored_genes.vegas.txt | sed 's/Top-0.1-pvalue/Pvalue/' >scored_genes.top10.txt
run_sigmod --sigmod SigMod_v2 --vegas scored_genes.top10.txt --tab2 ../preprocessing/BIOGRID-ORGANISM-Homo_sapiens-3.5.172.tab2.hgnc.txt -profile bigmem -resume 

#wget https://github.com/YuanlongLiu/SigMod/raw/master/SigMod_v2.zip && unzip SigMod_v2.zip
run_sigmod --sigmod SigMod_v2 --vegas ../preprocessing/scored_genes.vegas.txt --tab2 ../preprocessing/BIOGRID-ORGANISM-Homo_sapiens-3.5.172.tab2.hgnc.txt -profile bigmem -resume 

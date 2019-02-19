wget https://github.com/YuanlongLiu/SigMod/raw/master/SigMod_v2.zip && unzip SigMod_v2.zip
run_sigmod --sigmod SigMod_v2 --vegas ../preprocessing/scored_genes.vegas.txt --tab2 ../preprocessing/BIOGRID-MV-Physical-3.5.168.tab2.hgnc.tsv -profile bigmem -resume 

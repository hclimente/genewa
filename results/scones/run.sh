run_scones --bfile filtered --network gs -resume -profile cluster && mv cones.tsv cones_gs.tsv
run_scones --bfile filtered --network gm --snp2gene ../preprocessing/snp2hgnc.tsv -profile bigmem -resume && mv cones.tsv cones_gm.tsv
run_scones --bfile filtered --network gi --snp2gene ../preprocessing/snp2hgnc.tsv --tab2 ../preprocessing/BIOGRID-ORGANISM-Homo_sapiens-3.5.172.tab2.hgnc.txt -profile bigmem -resume && mv cones.tsv cones_gi.tsv

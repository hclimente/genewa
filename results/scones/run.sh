run_old_scones --bfile ../preprocessing/genesis_2019 --network gs -resume -profile bigmem && mv cones.tsv cones_gs.tsv
run_old_scones --bfile ../preprocessing/genesis_2019 --network gm --snp2gene ../preprocessing/snp2hgnc.tsv -profile bigmem -resume && mv cones.tsv cones_gm.tsv
run_old_scones --bfile ../preprocessing/genesis_2019 --network gi --snp2gene ../preprocessing/snp2hgnc.tsv --tab2 ../preprocessing/BIOGRID-ORGANISM-Homo_sapiens-3.5.172.tab2.hgnc.txt -profile bigmem -resume -with-trace && mv cones.tsv cones_gi.tsv

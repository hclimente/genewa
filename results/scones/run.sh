run_scones --bfile genesis_2019 --network gs -resume && mv cones.tsv cones_gs.tsv
run_scones --bfile genesis_2019 --network gm --snp2gene snp2hgnc.tsv -resume && mv cones.tsv cones_gm.tsv
run_scones --bfile genesis_2019 --network gi --snp2gene snp2hgnc.tsv --tab2 BIOGRID-ORGANISM-Homo_sapiens-3.4.160.tab.txt -resume && mv cones.tsv cones_gi.tsv

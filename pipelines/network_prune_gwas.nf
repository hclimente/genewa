#!/usr/bin/env nextflow
// files and working directories
params.wd = "."
params.out = "."
params.genewawd = "/Users/hclimente/projects/genewa"

wd = params.wd
genewawd = params.genewawd
srcReadPed = file("$genewawd/martiniflow/io/read_ped.nf")
srcGetNetwork = file("$genewawd/martiniflow/io/get_network.nf")
srcPruneGwas = file("$genewawd/martiniflow/io/prune_gwas.nf")

// GWAS files
ped = file("$genewawd/${params.geno}.ped")
map = file("$genewawd/${params.geno}.map")
snp2gene = file("$genewawd/$params.snp2gene")
tab = file("$genewawd/$params.tab")

params.nets = "ppi"
nets = params.nets.split(",")

process readData {

    input:
        file srcReadPed
        file ped
        file map

    output:
        file "gwas*.RData" into rgwas

    """
    nextflow run $srcReadPed --ped $ped --map $map -profile bigmem
    """

}

process getNetwork {

    input:
        file srcGetNetwork
        val net from nets
        file rgwas
        file snp2gene
        file tab

    output:
        file "net*.RData" into rnet

    """
    nextflow run $srcGetNetwork --gwas $rgwas --net $net --snp2gene $snp2gene --tab $tab -profile bigmem
    """

}

process pruneGwas {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file ped
        file map
        file rnet
        file tab

    output:
        file "${ped.baseName}.pruned.ped" into ped_pruned
        file "${map.baseName}.pruned.map" into map_pruned

    """
    nextflow run $srcPruneGwas --ped $ped --map $map --rnet $rnet --tab $tab -profile cluster
    """

}

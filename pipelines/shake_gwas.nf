#!/usr/bin/env nextflow
// files and working directories
params.wd = "."
params.out = "."
params.genewawd = "/Users/hclimente/projects/genewa"

wd = params.wd
genewawd = params.genewawd
srcReadPed = file("$genewawd/martiniflow/io/read_ped.nf")
srcGetNetwork = file("$genewawd/martiniflow/io/get_network.nf")
srcRunEvo = file("$genewawd/martiniflow/analyze/quick_evo.nf")

// GWAS files
ped = file("$genewawd/${params.geno}.ped")
map = file("$genewawd/${params.geno}.map")
snp2gene = file("$genewawd/$params.snp2gene")
tab = file("$genewawd/$params.tab")

// evo params
params.encoding = "additive"

nets = ["gs","gm","gi"]
associationScore = "chi2"
modelScore = "consistency"
encoding = params.encoding

process readData {

    input:
        file srcReadPed
        file ped
        file map

    output:
        file "gwas.RData" into rgwas

    """
    nextflow run $srcReadPed --ped $ped --map $map -profile bigmem
    """

}

rgwas.into { rgwas_getNetwork; rgwas_evo }

process getNetwork {

    input:
        file srcGetNetwork
        val net from nets
        file rgwas_getNetwork
        file snp2gene
        file tab

    output:
        file "net.RData" into rnets

    """
    nextflow run $srcGetNetwork --gwas $rgwas_getNetwork --net $net --snp2gene $snp2gene --tab $tab -profile bigmem
    """

}

process run_evo {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file rgwas_evo
        file rnet from rnets
        file srcRunEvo

    output:
        file "cones.*.RData" into rcones
        file "cones.*.tsv" into cones

    """
	nextflow run $srcRunEvo --rgwas $rgwas_evo --rnet $rnet --associationScore $associationScore --modelScore $modelScore --encoding $encoding -profile bigmem
    """

}

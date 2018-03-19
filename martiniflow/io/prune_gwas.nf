#!/usr/bin/env nextflow

params.out = "."

map = file("$params.map")
ped = file("$params.ped")
tab = file("$params.tab")
rnet = file("$params.rnet")

process get_snps {

    input:
        file rnet
        file tab

    output:
        file "snps.tsv" into snps

    """
    #!/usr/bin/env Rscript
    library(martini)
    library(tidyverse)

    load("$rnet")

    tab <- read_tsv("$tab")
    genes <- unique(c(tab\$OFFICIAL_SYMBOL_FOR_A, tab\$OFFICIAL_SYMBOL_FOR_B))

    snps <- martini:::subvert(net, 'gene', genes)
    snps <- data.frame(snps = names(snps))

    write_tsv(snps, path = 'snps.tsv')
    """

}

process filter_snps {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file map
        file ped
        file snps

    output:
        file "${ped.baseName}.pruned.ped" into pruned_ped
        file "${map.baseName}.pruned.map" into pruned_map

    """
    plink --file ${ped.baseName} --extract $snps --recode --out ${ped.baseName}.pruned
    """

}

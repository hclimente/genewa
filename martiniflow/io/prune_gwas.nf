#!/usr/bin/env nextflow

params.out = "."

map = file("$params.map")
ped = file("$params.ped")
tab2 = file("$params.tab2")
rnet = file("$params.rnet")

process get_snps {

    input:
        file rnet
        file tab2

    output:
        file "snps.tsv" into snps

    """
    #!/usr/bin/env Rscript
    library(martini)
    library(tidyverse)

    load("$rnet")

    tab2 <- read_tsv("$tab2")
    genes <- unique(c(tab2\$`Official Symbol Interactor A`, tab2\$`Official Symbol Interactor B`))

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

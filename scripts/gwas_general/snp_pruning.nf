#!/usr/bin/env nextflow

ped = file("${params.input}.ped")
map = file("${params.input}.map")

process prune {

    publishDir ".", overwrite: true, mode: 'copy'

    input:
        file ped
        file map

    output:
        file "plink.prune.in" into prunedIn
        file "plink.prune.out" into prunedOut
        file "${ped.baseName}.pruned.ped" into prunedPed
        file "${ped.baseName}.pruned.map" into prunedMap

    """
    plink --file ${ped.baseName} --filter-controls --indep-pairwise 50 1 0.8
    plink --file ${ped.baseName} --extract plink.prune.in --recode --out ${ped.baseName}.pruned
    """

}

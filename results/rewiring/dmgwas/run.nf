#!/usr/bin/env nextflow

tabs = file('../nets/rewired_*.tab2')
vegas = file('scored_genes.top10.txt')

process dmgwas {

    publishDir ".", overwrite: true, mode: "copy"

    input:
        file tab2 from tabs
        file vegas

    output:
        file "dmgwas_*.txt"

    """
    dmgwas.nf --vegas $vegas --tab2 $tab2 -resume
    mv selected_genes.dmgwas.txt dmgwas_`echo $tab2 | sed 's/.tab2//' | sed 's/rewired_//'`.txt
    """

}

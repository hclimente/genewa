#!/usr/bin/env nextflow

tabs = file('../nets/rewired_*.tab2')
vegas = file('scored_genes.top10.txt')
sigmod = file('SigMod_v2')

process sigmod {

    publishDir ".", overwrite: true, mode: "copy"

    input:
        file tab2 from tabs
        file vegas
        file sigmod

    output:
        file "sigmod_*.txt"

    """
    sigmod.nf --sigmod $sigmod --vegas $vegas --tab2 $tab2 -resume
    mv selected_genes.sigmod.txt sigmod_`echo $tab2 | sed 's/.tab2//' | sed 's/rewired_//'`.txt 
    """

}

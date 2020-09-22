#!/usr/bin/env nextflow

tabs = file('../nets/rewired_*.tab2')
vegas = file('scored_genes.top10.txt')

process lean {

    publishDir ".", overwrite: true, mode: "copy"

    input:
        file tab2 from tabs
        file vegas

    output:
        file "lean_*.txt"

    """
    lean.nf --vegas $vegas --tab2 $tab2 -resume
    mv scored_genes.lean.txt lean_`echo $tab2 | sed 's/.tab2//' | sed 's/rewired_//'`.txt
    """

}

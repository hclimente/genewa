#!/usr/bin/env nextflow

tab2 = file('../../preprocessing/hint.ht_complex.hgnc.pseudo.tab2')
vegas = file('scored_genes.top10.txt')

rs = Channel. from ( 0.01, 0.05, 0.1, 0.25, 0.5, 1 )
ds = Channel. from ( 1, 2 )

process dmgwas {

    publishDir ".", overwrite: true, mode: "copy"

    input:
        each r from rs
        each d from ds
        file tab2
        file vegas

    output:
        file "dmgwas_d_${d}_r_${r}.txt"

    """
    dmgwas.nf --vegas $vegas --tab2 $tab2 --d $d --r $r -profile bigmem -resume
    mv selected_genes.dmgwas.txt dmgwas_d_${d}_r_${r}.txt
    """

}

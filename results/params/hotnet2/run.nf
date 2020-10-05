#!/usr/bin/env nextflow

tab2 = file('../../preprocessing/hint.ht_complex.hgnc.pseudo.tab2')
vegas = file('scored_genes.top10.txt')
hotnet2 = file('hotnet2')

cutoffs = Channel. from ( 0.001, 0.01, 0.05, 0.125, 0.25, 0.5 )

process hotnet2 {

    publishDir ".", overwrite: true, mode: "copy"

    input:
        each cutoff from cutoffs
        file tab2
        file vegas
        file hotnet2

    output:
        file "hotnet2_${cutoff}.tsv"

    """
    hotnet2.nf --scores $vegas --tab2 $tab2 --hotnet2_path $hotnet2 --lfdr_cutoff $cutoff -resume -profile cluster
    mv selected_genes.hotnet2.tsv hotnet2_${cutoff}.tsv
    """

}

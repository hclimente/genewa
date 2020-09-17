#!/usr/bin/env nextflow

snp2gene = file('../../preprocessing/snp2hgnc.tsv')
tabs = file('../nets/rewired_*.tab2')
bed = file('../../preprocessing/genesis_2019.bed')
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")

process scones {
    
    publishDir ".", overwrite: true, mode: "copy"

    input:
        file bed
        file bim
        file fam
        file tab2 from tabs
        file snp2gene

    output:
        file "cones_*.tsv"

    """
    old_scones.nf --bfile $bed.baseName --network gi --snp2gene $snp2gene --tab2 $tab2 -resume
    mv cones.tsv cones_`echo $tab2 | sed 's/.tab2//' | sed 's/rewired_//'`.tsv
    """

}

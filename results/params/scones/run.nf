#!/usr/bin/env nextflow

snp2gene = file('../../preprocessing/snp2hgnc.tsv')
tab2 = file('../../preprocessing/hint.ht_complex.hgnc.pseudo.tab2')
bed = file('../../preprocessing/genesis_2019.bed')
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")

values = file(params.values)

etas = values. readLines()
lambdas = values. readLines() 

process scones {
    
    tag { "e $eta l $lambda" }
    publishDir ".", overwrite: true, mode: "copy"

    input:
        each eta from etas
        each lambda from lambdas
        file bed
        file bim
        file fam
        file tab2
        file snp2gene

    output:
        file "cones_eta_${eta}_lambda_${lambda}.tsv"

    """
    old_scones.nf --bfile $bed.baseName --network gi --snp2gene $snp2gene --tab2 $tab2 --eta $eta --lambda $lambda -resume
    mv cones.tsv cones_eta_${eta}_lambda_${lambda}.tsv
    """

}

#!/usr/bin/env nextflow

tab2 = file('../../preprocessing/hint.ht_complex.hgnc.pseudo.tab2')
vegas = file('scored_genes.top10.txt')
sigmod = file('SigMod_v2')

max_lambdas = Channel. from ( 1 )
n_maxs = Channel. from ( 10, 50, 100, 300, 700, 1000, 10000 )
max_jumps = Channel. from ( 5, 10, 20, 30, 50 )

process sigmod {

    publishDir ".", overwrite: true, mode: "copy"
    errorStrategy 'ignore'

    input:
        each nmax from n_maxs
        each maxjump from max_jumps
        each lambdamax from max_lambdas
        file tab2
        file vegas
        file sigmod

    output:
        file "sigmod_lambdamax_${lambdamax}_nmax_${nmax}_maxjump_${maxjump}.txt"

    """
    sigmod.nf --sigmod $sigmod --vegas $vegas --tab2 $tab2 --lambdamax $lambdamax --nmax $nmax --maxjump $maxjump
    mv selected_genes.sigmod.txt sigmod_lambdamax_${lambdamax}_nmax_${nmax}_maxjump_${maxjump}.txt
    """

}

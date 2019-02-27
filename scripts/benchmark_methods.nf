#!/usr/bin/env nextflow

params.out = '.'
params.k = 5
params.covar = ''

// gwas
bed = file("${params.bfile}.bed")
bim = file("${bed.baseName}.bim")
fam = file("${bed.baseName}.fam")
covar = file(params.covar)

scones_nets = ['gs','gm','gi']

// annotation
tab2 = file(params.tab2)
snp2gene = file(params.snp2gene)

//  SPLIT PREPARATION
/////////////////////////////////////
process make_splits {

    input:
        file FAM from fam
        val I from 1..params.k
        val K from params.k

    output:
        file 'samples_*.txt' into splits

    script:
    template 'genotypes/train_test_split.sh'

}

splits .into {splits_vegas; splits_scones}

process run_vegas {

    input:
        file BED from bed
        file FAM from fam
        file BIM from bim
        file SPLIT from splits_vegas
        file COVAR from covar

    output:
        file 'scored_genes.vegas.txt' into vegas_split

    """
    plink --bfile ${BED.baseName} --keep ${SPLIT} --make-bed --out input
    run_vegas --bfile input --genome GRCh37 --covar ${COVAR}
    """

}

//  BIOMARKER SELECTION
/////////////////////////////////////
process run_scones {

    tag { "${NET} (${SPLIT})" }

    input:
        file BED from bed
        file FAM from fam
        file BIM from bim
        file SPLIT from splits_scones
        file COVAR from covar
        val NET from scones_nets
        file TAB2 from tab2
        file SNP2GENE from snp2gene

    output:
        set ${SPLIT}, 'cones.tsv' into scones_biomarkers

    """
    plink --bfile ${BED.baseName} --keep ${SPLIT} --make-bed --out input
    run_scones --bfile input --covar ${COVAR} --network ${NET} --snp2gene ${SNP2GENE} --tab2 ${TAB2} -profile bigmem
    """

}

//  RISK COMPUTATION
/////////////////////////////////////

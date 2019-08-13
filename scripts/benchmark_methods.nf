#!/usr/bin/env nextflow

params.out = '.'
params.k = 5

// gwas
bed = file("${params.bfile}.bed")
bim = file("${bed.baseName}.bim")
fam = file("${bed.baseName}.fam")

scones_nets = ['gs','gm','gi']
// annotation
tab2 = file(params.tab2)
snp2gene = file(params.snp2gene)
r_ensg2hgnc = file(params.r_ensg2hgnc)
ensg2hgnc = file(params.ensg2hgnc)

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

splits .into {splits_vegas; splits_nothing; splits_scones}

process vegas {

    input:
        file BED from bed
        file FAM from fam
        file BIM from bim
        file SPLIT from splits_vegas
        file R_ENSG2HGNC from r_ensg2hgnc
        file ENSG2HGNC from ensg2hgnc

    output:
        set file(SPLIT), 'scored_genes.vegas.txt' into vegas

    """
    plink --bfile ${BED.baseName} --keep ${SPLIT} --make-bed --out input
    vegas2.nf --bfile input --genome 37 --gencode 31 --buffer 50000 --vegas_params '-top 10' -profile bigmem
    ./${R_ENSG2HGNC} scored_genes.vegas.txt ${ENSG2HGNC}
    """

}

vegas .into {vegas_sigmod; vegas_lean; vegas_heinz; vegas_dmgwas}

//  BIOMARKER SELECTION
/////////////////////////////////////
process do_nothing {

    tag { "${SPLIT}" }

    input:
        file BIM from bim
        file FAM from fam
        file SPLIT from splits_nothing

    output:
        set val("all_snps"), file(SPLIT), 'snps' into all_snps

    """
    echo snp >snps
    cut -f2 $BIM >>snps
    """

}

process scones {

    tag { "${NET}, ${SPLIT}" }

    input:
        file BED from bed
        file FAM from fam
        file BIM from bim
        file SPLIT from splits_scones
        each NET from scones_nets
        file TAB2 from tab2
        file SNP2GENE from snp2gene

    output:
        set val("scones_${NET}"), file(SPLIT), 'snps' into scones_snps, scones_snps_stability

    """
    plink --bfile ${BED.baseName} --keep ${SPLIT} --make-bed --out input
    old_scones.nf --bfile input --network ${NET} --snp2gene ${SNP2GENE} --tab2 ${TAB2} -profile bigmem
    echo snp >snps
    grep TRUE cones.tsv | cut -f1 >>snps
    """

}

process sigmod {

    tag { "${SPLIT}" }

    input:
        set file(SPLIT), file(VEGAS) from vegas_sigmod
        file TAB2 from tab2
        file SNP2GENE from snp2gene

    output:
        set val('sigmod'), file(SPLIT), 'snps' into sigmod_snps, sigmod_snps_stability
    
    """
    wget https://github.com/YuanlongLiu/SigMod/raw/master/SigMod_v2.zip && unzip SigMod_v2.zip
    cut -f2,9 ${VEGAS} | sed 's/Top-0.1-pvalue/Pvalue/' >scored_genes.top10.txt
    sigmod.nf --bfile input --sigmod SigMod_v2 --vegas scored_genes.top10.txt --tab2 ${TAB2} -profile cluster 
    R -e 'library(tidyverse); snp2gene <- read_tsv("${SNP2GENE}"); read_tsv("selected_genes.sigmod.txt") %>% inner_join(snp2gene, by = "gene") %>% select(snp) %>% write_tsv("snps")'
    """

}

process lean {

    tag { "${SPLIT}" }

    input:
        set file(SPLIT), file(VEGAS) from vegas_lean
        file TAB2 from tab2
        file SNP2GENE from snp2gene

    output:
        set val('lean'), file(SPLIT), 'snps' into lean_snps, lean_snps_stability
    
    """
    cut -f2,9 ${VEGAS} | sed 's/Top-0.1-pvalue/Pvalue/' >scored_genes.top10.txt
    lean.nf --vegas scored_genes.top10.txt --tab2 ${TAB2} -profile cluster
    R -e 'library(tidyverse); snp2gene <- read_tsv("${SNP2GENE}"); read_tsv("scored_genes.lean.txt") %>% filter(PLEAN < 0.05) %>% inner_join(snp2gene, by = c("Gene" = "gene")) %>% select(snp) %>% write_tsv("snps")'
    """

}

process heinz {

    tag { "${SPLIT}" }

    input:
        set file(SPLIT), file(VEGAS) from vegas_heinz
        file TAB2 from tab2
        file SNP2GENE from snp2gene

    output:
        set val('heinz'), file(SPLIT), 'snps' into heinz_snps, heinz_snps_stability

    """
    cut -f2,9 ${VEGAS} | sed 's/Top-0.1-pvalue/Pvalue/' >scored_genes.top10.txt
    heinz.nf --vegas scored_genes.top10.txt --tab2 ${TAB2} --fdr 0.5 -profile cluster -resume 
    R -e 'library(tidyverse); snp2gene <- read_tsv("${SNP2GENE}"); read_tsv("selected_genes.heinz.txt") %>% inner_join(snp2gene, by = "gene") %>% select(snp) %>% write_tsv("snps")'
    """

}

process dmgwas {

    tag { "${SPLIT}" }

    input:
        set file(SPLIT), file(VEGAS) from vegas_dmgwas
        file TAB2 from tab2
        file SNP2GENE from snp2gene

    output:
        set val('dmgwas'), file(SPLIT), 'snps' into dmgwas_snps, dmgwas_snps_stability
    
    """
    cut -f2,9 ${VEGAS} | sed 's/Top-0.1-pvalue/Pvalue/' >scored_genes.top10.txt
    dmgwas.nf --vegas scored_genes.top10.txt --tab2 ${TAB2} -profile cluster
    R -e 'library(tidyverse); snp2gene <- read_tsv("${SNP2GENE}"); read_tsv("selected_genes.dmgwas.txt") %>% inner_join(snp2gene, by = "gene") %>% select(snp) %>% write_tsv("snps")'
    """

}
/*
process hierarchichal_hotnet {

    tag { "${SPLIT}" }

    input:
        set file(SPLIT), file(VEGAS) from vegas_hotnet
        file TAB2 from tab2
        file SNP2GENE from snp2gene

    output:
        set val('hotnet'), file(SPLIT), 'snps' into hotnet_snps, hotnet_snps_stability
   , dmgwas_snps 
    """
    cut -f2,9 ${VEGAS} | sed 's/Top-0.1-pvalue/Pvalue/' >scored_genes.top10.txt
    git clone https://github.com/raphael-group/hierarchical-hotnet.git
    run_hhotnet --scores scored_genes.top10.txt --tab2 ${TAB2} --hhnet_path hierarchical-hotnet/src -profile bigmem
    R -e 'library(tidyverse); snp2gene <- read_tsv("${SNP2GENE}"); read_tsv("selected_genes.hotnet.txt") %>% inner_join(snp2gene, by = "gene") %>% select(snp) %>% write_tsv("snps")'
    """

}
*/
//  RISK COMPUTATION
/////////////////////////////////////
snps = scones_snps .mix( sigmod_snps, lean_snps, heinz_snps, dmgwas_snps, all_snps )

process lasso {

    errorStrategy 'ignore'
    tag { "${METHOD}, ${SPLIT}" }
    memory '40 GB'
    executor 'slurm'

    input:
        file BED from bed
        file FAM from fam
        file BIM from bim
        set val(METHOD), file(SPLIT), file(SNPS) from snps

    output:
        file 'bm' into bm

    """
    #!/usr/bin/env Rscript
    library(biglasso)
    library(snpStats)
    library(tidyverse)
    library(caret)

    # read dataset
    gwas <- read.plink("${BED}", "${BIM}", "${FAM}")
    X <- as(gwas[['genotypes']], 'numeric')
    ## remove NAs
    X[is.na(X)] <- 0	
    y <- gwas[['fam']][['affected']] - 1
    names(y) <- gwas[['fam']][['pedigree']] %>% as.character

    # read selected and splits
    selected <- read_tsv('${SNPS}', col_types = 'c')\$snp
    train <- read_delim('${SPLIT}', delim = ' ', col_names = FALSE, col_types = 'cc')\$X1
    test <- setdiff(as.character(gwas[['fam']][['pedigree']]), train)

    X_train <- X[train, selected] %>% as.big.matrix
    X_test <- X[test, selected] %>% as.big.matrix
    y_train <- y[train]
    y_test <- y[test] %>% as.factor
    rm(gwas, X, y)

    # train and evaluate classifier
    cvfit <- cv.biglasso(X_train, y_train, penalty = 'lasso', family = "binomial")
    y_pred <- predict(cvfit, X_test, type = "class") %>% as.numeric %>% as.factor

    tibble(method = "${METHOD}",
           n_selected = length(selected),
           n_active_set = sum(cvfit\$fit\$beta[,cvfit\$lambda == cvfit\$lambda.min][-1] != 0),
           sensitivity = sensitivity(y_pred, y_test),
           specificity = specificity(y_pred, y_test)) %>%
        write_tsv('bm')
    """

}

process join_analyses {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file "bm*" from bm. collect()

    output:
        file "prediction.tsv"

    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    lapply(list.files(pattern = 'bm*'), read_tsv, col_types = 'ciidd') %>% 
        do.call(rbind, .) %>%
        as.tibble %>%
        write_tsv('prediction.tsv')
    """

}

//  JACCARD
/////////////////////////////////////
snps_stability = scones_snps_stability 
    .mix( sigmod_snps_stability, lean_snps_stability, heinz_snps_stability, dmgwas_snps_stability )
    .groupTuple(size: params.k)

process stability {

    input:
        set val(METHOD), file('samples*'), file('snps*') from snps_stability

    output:
        file 'stability_method' into stability_method

    """
    #!/usr/bin/env python

    import csv
    import numpy as np
    from glob import glob

    def read_snps(path):
        FILE = open(path, 'r')
        FILE.readline()
        snps = FILE.read()
        snps = snps.split('\\n')
        snps = set(snps)
        FILE.close()

        return(snps)

    snps_files = glob('snps*')
    selected_snps = []
    jaccards = []
    for i in range(len(snps_files)):

        snps1 = read_snps(snps_files[i])

        for j in range(len(selected_snps)):
            snps2 = read_snps(snps_files[j])
            if len(snps1) and len(snps2):
                intersection = snps1 & snps2
                union = snps1 | snps2
                J = len(intersection)/len(union)
                jaccards.append(('{}-{}'.format(i, j), J))
            else:
                jaccards.append(('{}-{}'.format(i, j), 'nan'))

        selected_snps.append(snps1)

    with open('stability_method', 'w', newline='') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        tsv_output.writerow(['method', 'num_replicates', 'idx', 'jaccard'])
        for idx, J in jaccards:
            row = ['${METHOD}', len(snps_files), idx, J ]
            tsv_output.writerow(row)
    """

}

process join_stability {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file "stability_method*" from stability_method. collect()

    output:
        file "stability.tsv"

    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    lapply(list.files(pattern = 'stability_method*'), read_tsv, col_types = 'cicd') %>% 
        do.call(rbind, .) %>%
        as.tibble %>%
        write_tsv('stability.tsv')
    """

}

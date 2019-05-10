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
        set file(SPLIT), 'scored_genes.vegas.txt' into vegas_split

    """
    plink --bfile ${BED.baseName} --keep ${SPLIT} --make-bed --out input
    run_vegas --bfile input --genome GRCh37 --covar ${COVAR} --vegas_params '\\-top 10 -upper 50000 -lower 50000' -profile bigmem
    """

}

vegas_split .into {vegas_sigmod; vegas_lean}

//  BIOMARKER SELECTION
/////////////////////////////////////
process run_scones {

    tag { "${NET}, ${SPLIT}" }

    input:
        file BED from bed
        file FAM from fam
        file BIM from bim
        file SPLIT from splits_scones
        file COVAR from covar
        each NET from scones_nets
        file TAB2 from tab2
        file SNP2GENE from snp2gene

    output:
        set val("scones_${NET}"), file(SPLIT), 'snps' into scones_biomarkers

    """
    plink --bfile ${BED.baseName} --keep ${SPLIT} --make-bed --out input
    run_scones --bfile input --covar ${COVAR} --network ${NET} --snp2gene ${SNP2GENE} --tab2 ${TAB2} -profile bigmem
    grep TRUE cones.tsv | cut -f1 >snps
    """

}

process run_sigmod {

    tag { "${SPLIT}" }

    input:
        set file(SPLIT), file(VEGAS) from vegas_sigmod
        file TAB2 from tab2
        file SNP2GENE from snp2gene

    output:
        set val('sigmod'), file(SPLIT), 'snps' into sigmod_biomarkers
    
    """
    wget https://github.com/YuanlongLiu/SigMod/raw/master/SigMod_v2.zip && unzip SigMod_v2.zip
    run_sigmod --bfile input --sigmod SigMod_v2 --vegas ${VEGAS} --tab2 ${TAB2} -profile bigmem
    R -e 'library(tidyverse); snp2gene <- read_tsv("${SNP2GENE}"); read_tsv("selected_genes.sigmod.txt") %>% inner_join(snp2gene, by = "gene") %>% select(snp) %>% write_tsv("snps")'
    """

}

process run_lean {

    tag { "${SPLIT}" }

    input:
        set file(SPLIT), file(VEGAS) from vegas_lean
        file TAB2 from tab2
        file SNP2GENE from snp2gene

    output:
        set val('lean'), file(SPLIT), 'snps' into lean_biomarkers
    
    """
    run_lean --vegas ${VEGAS} --tab2 ${TAB2} -profile cluster
    R -e 'library(tidyverse); snp2gene <- read_tsv("${SNP2GENE}"); read_tsv("scored_genes.lean.txt") %>% inner_join(snp2gene, by = c("Gene" = "gene")) %>% select(snp) %>% write_tsv("snps")'
    """

}

//  RISK COMPUTATION
/////////////////////////////////////
biomarkers = scones_biomarkers .mix( sigmod_biomarkers, lean_biomarkers ) 

process benchmark {

    errorStrategy 'ignore'
    tag { "${METHOD}, ${SPLIT}" }

    input:
        file BED from bed
        file FAM from fam
        file BIM from bim
        file COVAR from covar
        set val(METHOD), file(SPLIT), file(SNPS) from biomarkers

    output:
        file 'auc' into auc

    """
    #!/usr/bin/env Rscript
    library(biglasso)
    library(snpStats)
    library(tidyverse)
    library(pROC)

    # read dataset
    gwas <- read.plink("${BED}", "${BIM}", "${FAM}")
    X <- as(gwas[['genotypes']], 'numeric') %>% as.big.matrix
    y <- gwas[['fam']][['phenotype']]

    # read selected and splits
    selected <- read_tsv('${SNPS}', col_names = FALSE)\$X1
    train <- read_delim('${SPLIT}', delim = ' ', col_names = FALSE)$X1
    test <- setdiff(gwas[['fam']][['pedigree']], train)

    X_train <- X[train, selected]
    X_test <- X[test, selected]
    y_train <- y[train]
    y_train <- y[test]
    rm(gwas, X, y)

    # train and evaluate classifier
    cvfit <- cv.biglasso(X_train, y_train, penalty = 'lasso', family = "binomial")
    y_pred <- predict(cvfit, X_test)

    X_test <- X[gwas[['fam']][['phenotype']] %in% test,]
    y_train <- y[gwas[['fam']][['phenotype']] %in% test]
    y_pred <- predict(cvfit, X_test)

    tibble(method = "${METHOD}",
           n_selected = length(snps),
           n_active_set = sum(cvfit\$fit\$beta[,cvfit\$lambda == cvfit\$lambda.min][-1] != 0),
           auc = auc(y_pred, y_train)) %>%
        write_tsv('auc')
    """

}

process join_analyses {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file "auc*" from auc. collect()

    output:
        file "prediction.tsv"

    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    lapply(list.files(pattern = 'auc*'), read_tsv, col_types = 'ciid') %>% 
        do.call(rbind, .) %>%
        as.tibble %>%
        write_tsv('prediction.tsv')
    """

}

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

process vegas {

    input:
        file BED from bed
        file FAM from fam
        file BIM from bim
        file SPLIT from splits_vegas
        file COVAR from covar

    output:
        set file(SPLIT), 'scored_genes.vegas.txt' into vegas

    """
    plink --bfile ${BED.baseName} --keep ${SPLIT} --make-bed --out input
    run_vegas --bfile input --genome GRCh37 --covar ${COVAR} --vegas_params '\\-top 10 -upper 50000 -lower 50000' -profile bigmem
    """

}

//vegas .into {vegas_sigmod; vegas_lean; vegas_heinz; vegas_hotnet}
vegas .into {vegas_sigmod; vegas_lean; vegas_heinz}

//  BIOMARKER SELECTION
/////////////////////////////////////
process scones {

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
        set val('sigmod'), file(SPLIT), 'snps' into sigmod_biomarkers
    
    """
    wget https://github.com/YuanlongLiu/SigMod/raw/master/SigMod_v2.zip && unzip SigMod_v2.zip
    run_sigmod --bfile input --sigmod SigMod_v2 --vegas ${VEGAS} --tab2 ${TAB2} -profile bigmem
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
        set val('lean'), file(SPLIT), 'snps' into lean_biomarkers
    
    """
    run_lean --vegas ${VEGAS} --tab2 ${TAB2} -profile cluster
    R -e 'library(tidyverse); snp2gene <- read_tsv("${SNP2GENE}"); read_tsv("scored_genes.lean.txt") %>% inner_join(snp2gene, by = c("Gene" = "gene")) %>% select(snp) %>% write_tsv("snps")'
    """

}

process heinz {

    tag { "${SPLIT}" }

    input:
        set file(SPLIT), file(VEGAS) from vegas_heinz
        file TAB2 from tab2
        file SNP2GENE from snp2gene

    output:
        set val('heinz'), file(SPLIT), 'snps' into heinz_biomarkers

    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(igraph)
    library(BioNet)

    # search subnetworks
    vegas <- read_tsv('${VEGAS}') %>% 
        select(Gene, Pvalue) %>%
        mutate(score = qnorm(1 - (Pvalue/2)))
    scores <- vegas\$score
    names(scores) <- vegas\$Gene

    net <- read_tsv("${TAB2}") %>%
        rename(gene1 = `Official Symbol Interactor A`, 
               gene2 = `Official Symbol Interactor B`) %>%
        filter(gene1 %in% names(scores) & gene2 %in% names(scores)) %>%
        select(gene1, gene2) %>%
        graph_from_data_frame(directed = FALSE)

    selected <- runFastHeinz(net, scores)
    
    # map selected genes to snps
    snp2gene <- read_tsv("${SNP2GENE}")
    tibble(gene = names(V(selected))) %>% 
        inner_join(snp2gene, by = 'gene') %>% 
        select(snp) %>% 
        write_tsv("snps")
    """


}

/*
process hotnet {

    tag { "${SPLIT}" }

    input:
        set file(SPLIT), file(VEGAS) from vegas_hotnet
        file TAB2 from tab2
        file SNP2GENE from snp2gene

    output:
        set val('hotnet'), file(SPLIT), 'snps' into hotnet_biomarkers
    
    """
    git clone https://github.com/raphael-group/hierarchical-hotnet.git
    run_hhotnet --scores ${VEGAS} --tab2 ${TAB2} --hhnet_path hierarchical-hotnet/src -profile bigmem
    R -e 'library(tidyverse); snp2gene <- read_tsv("${SNP2GENE}"); read_tsv("selected_genes.sigmod.txt") %>% inner_join(snp2gene, by = "gene") %>% select(snp) %>% write_tsv("snps")'
    """

}
*/

//  RISK COMPUTATION
/////////////////////////////////////
// biomarkers = scones_biomarkers .mix( sigmod_biomarkers, lean_biomarkers, heinz_biomarkers, hotnet_biomarkers ) 
biomarkers = scones_biomarkers .mix( sigmod_biomarkers, lean_biomarkers, heinz_biomarkers ) 

process lasso {

    errorStrategy 'ignore'
    tag { "${METHOD}, ${SPLIT}" }
    clusterOptions = '-l mem=20G'

    input:
        file BED from bed
        file FAM from fam
        file BIM from bim
        file COVAR from covar
        set val(METHOD), file(SPLIT), file(SNPS) from biomarkers

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
    covars_df <- read_tsv('${COVAR}')
    covars <- select(covars_df, AGE) %>% as.matrix()
    rownames(covars) <- covars_df\$IID
    X <- as(gwas[['genotypes']], 'numeric')
    y <- gwas[['fam']][['affected']] - 1
    names(y) <- gwas[['fam']][['member']]

    # read selected and splits
    selected <- read_tsv('${SNPS}', col_types = 'c')\$snp
    train <- read_delim('${SPLIT}', delim = ' ', col_names = FALSE, col_types = 'cc')\$X1
    test <- setdiff(gwas[['fam']][['pedigree']], train)

    X_train <- X[train, selected] %>% cbind(covars[train,]) %>% as.big.matrix
    X_test <- X[test, selected] %>% cbind(covars[test,]) %>% as.big.matrix
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

    lapply(list.files(pattern = 'bm*'), read_tsv, col_types = 'ciid') %>% 
        do.call(rbind, .) %>%
        as.tibble %>%
        write_tsv('prediction.tsv')
    """

}

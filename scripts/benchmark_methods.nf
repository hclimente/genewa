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

    output:
        set file(SPLIT), 'scored_genes.vegas.txt' into vegas

    """
    plink --bfile ${BED.baseName} --keep ${SPLIT} --make-bed --out input
    run_vegas --bfile input --genome GRCh37 --gencode 31 --vegas_params '\\-top 10 -upper 50000 -lower 50000' -profile bigmem
    """

}

vegas .into {vegas_sigmod; vegas_lean; vegas_heinz; vegas_dmgwas}

//  BIOMARKER SELECTION
/////////////////////////////////////
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
        set val("scones_${NET}"), file(SPLIT), 'snps' into scones_biomarkers

    """
    plink --bfile ${BED.baseName} --keep ${SPLIT} --make-bed --out input
    run_old_scones --bfile input --network ${NET} --snp2gene ${SNP2GENE} --tab2 ${TAB2} -profile bigmem
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
    cut -f2,9 ${VEGAS} | sed 's/Top-0.1-pvalue/Pvalue/' >scored_genes.top10.txt
    run_sigmod --bfile input --sigmod SigMod_v2 --vegas scored_genes.top10.txt --tab2 ${TAB2} -profile cluster 
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
    cut -f2,9 ${VEGAS} | sed 's/Top-0.1-pvalue/Pvalue/' >scored_genes.top10.txt
    run_lean --vegas scored_genes.top10.txt --tab2 ${TAB2} -profile cluster
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
        set val('heinz'), file(SPLIT), 'snps' into heinz_biomarkers

    """
    cut -f2,9 ${VEGAS} | sed 's/Top-0.1-pvalue/Pvalue/' >scored_genes.top10.txt
    run_heinz --vegas scored_genes.top10.txt --tab2 ${TAB2} --fdr 0.5 -profile cluster -resume 
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
        set val('dmgwas'), file(SPLIT), 'snps' into dmgwas_biomarkers
    
    """
    cut -f2,9 ${VEGAS} | sed 's/Top-0.1-pvalue/Pvalue/' >scored_genes.top10.txt
    run_dmGWAS --vegas scored_genes.top10.txt --tab2 ${TAB2} -profile cluster
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
        set val('hotnet'), file(SPLIT), 'snps' into hotnet_biomarkers
   , dmgwas_biomarkers 
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
biomarkers = scones_biomarkers .mix( sigmod_biomarkers, lean_biomarkers, heinz_biomarkers, dmgwas_biomarkers ) 

process lasso {

    errorStrategy 'ignore'
    tag { "${METHOD}, ${SPLIT}" }
    clusterOptions = '-l mem=20G'

    input:
        file BED from bed
        file FAM from fam
        file BIM from bim
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

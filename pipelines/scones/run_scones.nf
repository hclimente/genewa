#!/usr/bin/env nextflow

ped = file("$params.ped")
map = file("$params.map")

// network information
snp2gene = file("$params.snp2gene")
icogsSnps = file("$params.icogsSnps")
ppi = file("$params.ppi")

process get_phenotypes {

  input:
    file ped
  output:
    file "phenotype.txt" into pheno

  """
  echo FID IID BRCA > phenotype.txt
  cut -d' ' -f1,2,6 $ped >> phenotype.txt
  """
}

process get_GS {

  input:
    file map
  output:
    file "gs.txt" into gs

  """
  #!/usr/bin/env Rscript

  library(readr)
  library(dplyr)
  library(magrittr)

  map <- read_tsv("$map", col_names = FALSE) %>%
    set_colnames(c("chr","gene","cm","pos")) %>%
    select(chr,pos) %>%
    unique %>%
    arrange(chr, pos)

  by(map, map\$chr, function(x){
      chr <- unique(x\$chr)
      pos1 <- head(x\$pos, n = length(x\$pos) - 1)
      pos2 <- tail(x\$pos, n = length(x\$pos) - 1)
      data.frame(chr1 = chr, pos1 = pos1, chr2 = chr, pos2 = pos2)
    }) %>% do.call("rbind", .) %>%
      write_tsv("gs.txt", col_names=FALSE)

  """

}

gs.into { gs; gs_tmp }

process get_GM {

  input:
    file snp2gene
    file map
    file gs from gs_tmp
  output:
    file "gm.txt" into gm

  """
  #!/usr/bin/env Rscript

  library(readr)
  library(dplyr)
  library(magrittr)

  genes <- read_tsv("$snp2gene") %>%
    set_colnames(c("snp","gene"))
  map <- read_tsv("$map", col_names = FALSE) %>%
    set_colnames(c("chr","snp","cm","pos")) %>%
    merge(genes) %>%
    select(chr,gene,pos) %>%
    # in some cases the same position is linked to two different variants
    unique

  gm <- by(map, map\$gene, function(x){
    chr <- unique(x\$chr)
    if (nrow(x) > 1){
      comb <- combn(x\$pos, 2)
      data.frame(chr1 = chr, pos1 = comb[1,], chr2 = chr, pos2 = comb[2,])
    }
  }) %>% do.call("rbind", .)

  read_tsv("$gs", col_names = FALSE) %>%
    set_colnames(c("chr1", "pos1", "chr2", "pos2")) %>%
    rbind(gm) %>%
    write_tsv("gm.txt", col_names=FALSE)

  """

}

gm.into { gm; gm_tmp }

process get_GI {

  clusterOptions = '-l mem=10G'

  input:
    file ppi
    file snp2gene
    file map
    file gm from gm_tmp
  output:
    file "gi.txt" into gi

  """
  #!/usr/bin/env Rscript

  library(readr)
  library(dplyr)
  library(magrittr)

  genes <- read_tsv("$snp2gene") %>%
    set_colnames(c("snp","gene"))

  map <- read_tsv("$map", col_names = FALSE) %>%
      set_colnames(c("chr","snp","cm","pos")) %>%
      select(chr,snp,pos)

  ppi <- read_tsv("$ppi") %>%
    select(OFFICIAL_SYMBOL_FOR_A,OFFICIAL_SYMBOL_FOR_B,ALIASES_FOR_A,ALIASES_FOR_B) %>%
    # CHECK IF THE NAME IS IN THE ALIASES
    select(OFFICIAL_SYMBOL_FOR_A,OFFICIAL_SYMBOL_FOR_B) %>%
    rename(gene1 = OFFICIAL_SYMBOL_FOR_A,
           gene2 = OFFICIAL_SYMBOL_FOR_B) %>%
    # remove self-interactions
    filter(gene1 != gene2) %>%
    unique

  gi <- merge(ppi, genes, by.x = "gene1", by.y = "gene") %>%
    merge(genes, by.x = "gene2", by.y = "gene") %>%
    merge(map, by.x = "snp.x", by.y = "snp") %>%
    merge(map, by.x = "snp.y", by.y = "snp") %>%
    select(chr.x, pos.x, chr.y, pos.y) %>%
    set_colnames(c("chr1", "pos1", "chr2", "pos2")) %>%
    # remove cases of the same snp mapped to interacting genes
    filter(chr1 != chr2 | pos1 != pos2)

  read_tsv("$gm", col_names = FALSE) %>%
    set_colnames(c("chr1", "pos1", "chr2", "pos2")) %>%
    rbind(gi) %>%
    write_tsv("gi.txt", col_names=FALSE)

  """

}

process run_scones {

  clusterOptions = '-l mem=30G'

  input:
    file ped from ped
    file map from map
    file phenotype from pheno.first()
    file net from gs .mix(gm) . mix(gi)
  output:
    file "BRCA_${net.baseName}.scones.out.txt" into selectedSnps
    file "BRCA_${net.baseName}.scones.pmatrix.txt" into pmatrix

  """
  scones2 -p ${ped.baseName} -f $phenotype -n $net
  mv BRCA.scones.pmatrix.txt BRCA_${net.baseName}.scones.pmatrix.txt
  mv BRCA.scones.out.txt BRCA_${net.baseName}.scones.out.txt
  """
}

process get_genes {

  input:
    file selectedSnps
    file snp2gene
    file icogsSnps
  output:
    file "candidateGenes.txt" into selectedGenes
    file "candidateRegions.txt" into selectedRegions

  """
  #!/usr/bin/env Rscript

  library(readr)
  library(dplyr)
  library(magrittr)

  selectedSnps <- read_tsv("$selectedSnps", skip = 2) %>%
    set_colnames(c("snp","chr","pos"))

  snps <- read_csv("$icogsSnps") %>%
      select(Illumina_SNP_Name,Chromosome,Build37_Position) %>%
      set_colnames(c("snp","chr","pos")) %>%
      arrange(chr, pos) %>%
      mutate(selected = ifelse(snp %in% selectedSnps\$snp, TRUE, FALSE))

  rlex <- rle(snps\$selected)
  consecutive <- rep(ifelse(rlex\$lengths > 1, rlex\$values, FALSE), rlex\$lengths)

  read_tsv("$snp2gene") %>%
    set_colnames(c("snp","gene")) %>%
    filter(snp %in% selectedSnps\$snp) %>%
    select(gene) %>%
    unique %>%
    write_tsv("candidateGenes.txt")

  lower <- 1
  regions <- list()
  for (i in 1:length(rlex\$lengths)){
    l <- rlex\$lengths[i]
    v <- rlex\$values[i]

    if (v & l > 1){
      upper <- lower + l -1
      first <- snps[lower,]
      last <- snps[upper,]
      regions[[i]] <- c(unique(c(first\$chr, last\$chr)), first\$pos, last\$pos, l)
    }
    lower <- lower + l
  }

  do.call("rbind",regions) %>%
    as.data.frame %>%
    set_colnames(c("chr","start","end","num")) %>%
    write_tsv("candidateRegions.txt")

  """

}

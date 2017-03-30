#!/usr/bin/env nextflow

params.tab = "ppi.tab"
params.map = "genotypes.map"

map = file("$params.map")
snp2gene = file("$params.snp2gene")
tab = file("$params.tab")

process get_GS {
  publishDir ".", overwrite: true

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

  gs <- by(map, map\$chr, function(x){
      chr <- unique(x\$chr)
      pos1 <- head(x\$pos, n = length(x\$pos) - 1)
      pos2 <- tail(x\$pos, n = length(x\$pos) - 1)
      data.frame(chr1 = chr, pos1 = pos1, chr2 = chr, pos2 = pos2)
    }) %>% do.call("rbind", .)

    gs %>%
      # make it bidirectional
      rename(chrn = chr1, posn = pos1, chr1 = chr2, pos1 = pos2) %>%
      rename(chr2 = chrn, pos2 = posn) %>%
      select(chr1, pos1, chr2, pos2) %>%
      rbind(gs) %>%
      write_tsv("gs.txt", col_names=FALSE)

  """

}

gs.into { gs; gs_tmp }

process get_GM {
  publishDir ".", overwrite: true

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

  # make it bidirectional
  gm <- gm %>%
    rename(chrn = chr1, posn = pos1, chr1 = chr2, pos1 = pos2) %>%
    rename(chr2 = chrn, pos2 = posn) %>%
    select(chr1, pos1, chr2, pos2) %>%
    rbind(gm)

  read_tsv("$gs", col_names = FALSE) %>%
    set_colnames(c("chr1", "pos1", "chr2", "pos2")) %>%
    rbind(gm) %>%
    unique %>%
    write_tsv("gm.txt", col_names=FALSE)

  """

}

gm.into { gm; gm_tmp }

process get_GI {
  publishDir ".", overwrite: true

  input:
    file tab
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

  ppi <- read_tsv("$tab") %>%
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

  # make it bidirectional
  gi <- gi %>%
    rename(chrn = chr1, posn = pos1, chr1 = chr2, pos1 = pos2) %>%
    rename(chr2 = chrn, pos2 = posn) %>%
    select(chr1, pos1, chr2, pos2) %>%
    rbind(gi)

  read_tsv("$gm", col_names = FALSE) %>%
    set_colnames(c("chr1", "pos1", "chr2", "pos2")) %>%
    rbind(gi) %>%
    unique %>%
    write_tsv("gi.txt", col_names=FALSE)

  """

}

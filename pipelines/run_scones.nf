#!/usr/bin/env nextflow

ped = file("$HOME/genewa/data/genesis/genesis.imputed.filtered.ped")
map = file("$HOME/genewa/data/genesis/genesis.imputed.filtered.map")
genes = file("$HOME/genewa/data/genesis/glist-hg19")
snp2gene = file("$HOME/genewa/data/genesis/gene2snp.hg19")
ppi = file("$HOME/genewa/data/genesis/BIOGRID-ORGANISM-Homo_sapiens-3.4.138.tab.txt")

scones = "$HOME/genewa/libs/easyGWASCore/bin/linux2/tools/scones"

process get_GS {
  input:
    file map
  output:
    file "gs.txt" into gs, gs_tmp

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

process get_GM {
  input:
    file snp2gene
    file map
    file gs_tmp
  output:
    file "gm.txt" into gm, gm_tmp

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

  read_tsv("$gs_tmp", col_names = FALSE) %>%
    set_colnames(c("chr1", "pos1", "chr2", "pos2")) %>%
    rbind(gm) %>%
    write_tsv("gm.txt", col_names=FALSE)

  """

}

process get_GI {
  input:
    file ppi
    file snp2gene
    file map
    file gm_tmp
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
    set_colnames(c("chr1", "pos1", "chr2", "pos2"))

  read_tsv("$gm_tmp", col_names = FALSE) %>%
    set_colnames(c("chr1", "pos1", "chr2", "pos2")) %>%
    rbind(gi) %>%
    write_tsv("gi.txt", col_names=FALSE)

  """

}

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

process run_scones {

  input:
    file ped from ped.first()
    file map from map.first()
    file phenotype from pheno.first()
    file net from gs .mix(gm) . mix(gi)
  output:
    file "BRCA_${net.baseName}.scones.out.txt" into out
    file "BRCA_${net.baseName}.scones.pmatrix.txt" into pmatrix

  """
  $scones ${ped.baseName} $phenotype $net 0.05 `pwd` additive 0
  mv BRCA.scones.pmatrix.txt BRCA_${net.baseName}.scones.pmatrix.txt
  mv BRCA.scones.out.txt BRCA_${net.baseName}.scones.out.txt
  """
}

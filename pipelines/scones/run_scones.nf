#!/usr/bin/env nextflow

// GENESIS data
peds = Channel.fromPath("$HOME/genewa/data/genesis/chr*.processed.ped")
maps = Channel.fromPath("$HOME/genewa/data/genesis/chr*.processed.map")
// ped = file("$HOME/genewa/data/genesis/Genesis.ped")
// map = file("$HOME/genewa/data/genesis/Genesis.map")

// network information
genes = file("$HOME/genewa/data/genesis/glist-hg19")
snp2gene = file("$HOME/genewa/data/genesis/gene2snp.hg19")
icogsSnps = file("$HOME/genewa/data/genesis/icogs_snp_list.csv")
ppi = file("$HOME/genewa/data/genesis/BIOGRID-ORGANISM-Homo_sapiens-3.4.138.tab.txt")

process join_files {

  input:
    file peds from peds.toList()
    file maps from maps.toList()
  output:
    file "genesis.processed.ped" into ped
    file "genesis.processed.map" into map

  """
  cat chr20.processed.ped | sed 's/[\\t ]\$//' >genesis.processed.ped
  cat chr20.processed.map >genesis.processed.map

  ###################
  # CORRECT ALL CHR, REPLACE 20 BY 1 #
  ###################
  for i in `seq 21 22`
  do
    paste -d' ' genesis.processed.ped <(cut -d' ' -f7- chr\$i.processed.ped | sed 's/[\\t ]\$//') > tmp.ped && mv tmp.ped genesis.processed.ped
    cat genesis.processed.map chr\$i.processed.map >tmp.map && mv tmp.map genesis.processed.map
  done
  """

}

ped.into { ped; ped_pheno }
map.into { map; map_gs; map_gm; map_gi }

/*
process numeric2acgt {
  input:
    file ped
    file map
  output:
    file "genesis.filtered.ped" into ped_filtered
    file "genesis.filtered.map" into map_filtered

  """
  # remove indels for the moment
  grep '\\[D/I\\]\\|\\[I/D\\]' ~/genewa/data/genesis/icogs_snp_list.csv | cut -d',' -f2 >indels.txt
  plink --file ${ped.baseName} --recode --alleleACGT --exclude indels.txt --out genesis.filtered --noweb
  """
}
*/

process get_phenotypes {

  input:
    file ped from ped_pheno
  output:
    file "phenotype.txt" into pheno

  """
  echo FID IID BRCA > phenotype.txt
  cut -d' ' -f1,2,6 $ped >> phenotype.txt
  """
}

process get_GS {

  input:
    file map from map_gs
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
    file map from map_gm
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
    file map from map_gi
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
    file ped from ped.first()
    file map from map.first()
    file phenotype from pheno.first()
    file net from gs .mix(gm) . mix(gi)
  output:
    file "BRCA_${net.baseName}.scones.out.txt" into selectedSnps
    file "BRCA_${net.baseName}.scones.pmatrix.txt" into pmatrix

  """
  scones ${ped.baseName} $phenotype $net 0.05 `pwd` additive 3
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

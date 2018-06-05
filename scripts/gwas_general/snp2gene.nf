#!/usr/bin/env nextflow

params.out = '.'

gff = file("$params.gff")
map = file("$params.map")
ensembl2hgnc = file("$params.hgnc")
biogrid = file("$params.biogrid")

process extractGenes {

	input:
		file gff

	output:
		file 'genes.gff' into genes_gff

	"""
	awk '{if(\$3=="gene" || \$3=="pseugene"){print \$0}}' $gff >genes.gff
	"""

}

process gff2bed {

	input:
		file genes_gff

	output:
		file 'genes.bed' into genes_bed

	"""
	gff2bed < $genes_gff >genes.bed
	"""

}

process map2bed {

	input:
		file map

	output:
		file 'snps.bed' into snps_bed

	"""
	awk '{print "chr" \$1 "\\t" \$4 "\\t" \$4 "\\t" \$2 "\\t.\\t." }' $map >tmp
	sed 's/chr23/chrX/' tmp >snps.bed
	"""

}

process snp2gene {

	input:
		file snps_bed
		file genes_bed

	output:
		file 'snp2ensembl.tsv' into snp2ensembl

	"""
	bedtools intersect -a $snps_bed -b $genes_bed -wa -wb >tmp
	cut -f4,16 tmp | gsed 's/\\tID=/\\t/' | gsed 's/\\.[0-9]\\+;.\\+//' >snp2ensembl.tsv
	"""

}

process ensembl2hgnc {

	publishDir "$params.out", overwrite: true, mode: "copy"

	input:
		file ensembl2hgnc
		file snp2ensembl

	output:
		file 'snp2hgnc.tsv' into snp2hgnc

	"""
	#!/usr/bin/env Rscript
	library(tidyverse)
	library(magrittr)

	ensembl2hgnc <- read_tsv('$ensembl2hgnc') %>%
		select(symbol, ensembl_gene_id)

	read_tsv('$snp2ensembl', col_names = FALSE) %>%
		set_colnames(c('snp','ensembl_gene_id')) %>%
		inner_join(ensembl2hgnc, by = 'ensembl_gene_id') %>%
		select(snp, symbol) %>%
		write_tsv('snp2hgnc.tsv')
	"""

}

process biogrid2hgnc {

	publishDir "$params.out", overwrite: true, mode: "copy"

	input:
		file ensembl2hgnc
		file biogrid

	output:
		file "${biogrid.baseName}.hgnc.tsv" into hgncBiogrid

	"""
	#!/usr/bin/env Rscript
	library(tidyverse)

	symbols <- read_tsv('$ensembl2hgnc') %>%
		select(entrez_id, symbol)

	read_tsv("$biogrid") %>%
		inner_join(symbols, by = c('Entrez Gene Interactor A' = 'entrez_id')) %>%
		inner_join(symbols, by = c('Entrez Gene Interactor B' = 'entrez_id')) %>%
		rename(`Official Symbol Interactor A` = symbol.x, `Official Symbol Interactor B` = symbol.y) %>%
		write_tsv('${biogrid.baseName}.hgnc.tsv')
	"""

}

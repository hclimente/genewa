snp2gene --bim genesis_raw_fix.bim --grch_version 37 --gencode_version 31 -profile cluster -resume
sed 's/symbol/gene/' snp2hgnc.tsv >snp2gene.preqc.tsv

rm snp2hgnc.tsv

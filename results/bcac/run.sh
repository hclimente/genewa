vegas2.nf --snp_association oncoarray_bcac_public_release_oct17.assoc.vegas --ld_controls g1000p3_EUR --genome 37 --gencode 31 --buffer 50000 --vegas_params '-top 10' -resume -profile bigmem "$@"
../../scripts/ensembl2hgnc.R scored_genes.vegas.txt ../preprocessing/non_alt_loci_set.txt

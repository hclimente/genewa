# cut -f2,9 ../../preprocessing/scored_genes.vegas.txt | sed 's/Top-0.1-pvalue/Pvalue/' >scored_genes.top10.txt

for i in `seq 0 0.05 1`
do
  heinz.nf --vegas scored_genes.top10.txt --tab2 ../../preprocessing/hint.ht_complex.hgnc.pseudo.tab2 --fdr $i -profile bigmem -resume
  mv selected_genes.heinz.txt heinz_$i.txt
done

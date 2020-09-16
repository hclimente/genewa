# cut -f2,9 ../../preprocessing/scored_genes.vegas.txt | sed 's/Top-0.1-pvalue/Pvalue/' >scored_genes.top10.txt

for i in `seq 1 100`
do
  heinz.nf --vegas scored_genes.top10.txt --tab2 ../nets/rewired_${i}.tab2 --fdr 0.5 -profile bigmem -resume
  mv selected_genes.heinz.txt heinz_$i.txt
done

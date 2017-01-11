unalias mv
unalias cp

########################################
#               IMPUTE                 #
########################################

# impute
for i in `seq 1 22`
do
  echo CHROMOSOME $i
  ./prepare_ped.nf -work-dir chr$i --chr $i -profile curie -resume
done

./prepare_ped.nf -work-dir chrX_NONPAR --chr 23 --chr_region X_NONPAR -profile curie -resume
./prepare_ped.nf -work-dir chrX_PAR1 --chr 23 --chr_region X_PAR1 -profile curie -resume
./prepare_ped.nf -work-dir chrX_PAR2 --chr 23 --chr_region X_PAR2 -profile curie -resume

# join chromosomes
## start by the X, as is the one we cannot use in the loop
cp ~/genewa/data/genesis/chrX_NONPAR.processed.ped genesis.processed.ped
cp ~/genewa/data/genesis/chrX_NONPAR.processed.map genesis.processed.map

for i in `seq 1 22`
do
  echo CHROMOSOME $i
  paste -d' ' genesis.processed.ped <(cut -f7- ~/genewa/data/genesis/chr$i.processed.ped -d' ') > tmp.ped
  mv tmp.ped genesis.processed.ped

  cat ~/genewa/data/genesis/chr$i.processed.map >>genesis.processed.map
done

########################################
#     STUDY POPULATION STRUCTURE       #
########################################

smartpca.perl -i genesis.processed.ped -a genesis.processed.map -b genesis.processed.ped -o genesis.processed.pca -p genesis.processed.plot -e genesis.processed.eval -l genesis.processed.pca.log
smarteigenstrat.perl -i genesis.processed.ped -a genesis.processed.map -b genesis.processed.ped -p genesis.processed.pca -l genesis.processed.eigenstrat.log -o genesis.processed.chisq
gc.perl genesis.processed.chisq genesis.processed.lambda

mv genesis.processed.??? ~/genewa/data/genesis

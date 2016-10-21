#!/bin/bash
#$ -cwd
#$ -V
# launch with qsub -t 1-11 scripts/launch_plink.sh

# get pedfile
param_file="populations/gametes/parameters.txt"

h=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -f1`
maf=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -f2`
N=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -f3`
modelNo=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -f4`
repNo=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -f5`

scripts/gametes2plink $h $maf $N $modelNo $repNo

# run plink
ped=populations/gametes/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".txt.ped

plink --noweb --epistasis --ped $ped --map populations/gametes/map_$N.txt --out $ped.plink --epi1 1 --epi2 1 --allow-no-sex

# launch sixpac
java -jar libs/Sixpac.jar --raw $ped --map populations/gametes/map.txt --out $ped.sixpac --mode complete

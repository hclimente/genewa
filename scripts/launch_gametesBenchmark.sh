#!/bin/bash
#$ -N gametes
#$ -cwd
#$ -V
#$ -S /bin/sh
# maximum number of nodes to take
#$ -tc 100

# generate pedfile
param_file="populations/gametes/parameters.txt"

h=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -f1`
maf=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -f2`
N=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -f3`
modelNo=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -f4`
repNo=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -f5`

echo Analyzing population h $h maf $maf N $N modelNo $modelNo repNo $repNo

# scripts/gametes2plink.R $h $maf $N $modelNo $repNo
ped=populations/gametes/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".txt.ped

# run plink
plinkOut=populations/gametes/plink/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".plink.txt
plink --noweb --epistasis --ped $ped --map populations/gametes/map_$N.txt --out $plinkOut --epi1 1 --epi2 1 --allow-no-sex

# run sixpac
# sixpacOut=populations/gametes/sixpac/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".sixpac.txt
# java -jar libs/Sixpac.jar --raw $ped --map populations/gametes/map.txt --out $sixpacOut --mode complete

# run MDR
mdrOut=populations/gametes/mdr/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".mdr.txt
gametesFile=libs/mdr_3.0.2/mdr_3.0.2.jar populations/gametes/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".txt
java -jar libs/mdr_3.0.2/mdr_3.0.2.jar -min=1 -max=2 -cv=10 -table_data=true $gametesFile | \
## crop the Top models table
awk '/### Top Models ###/{f=1} /FINISHED/{f=0;print} f' | \
## remove first and last line and blank lines
grep -v '###' | sed '/^$/d' | sed '$ d' | \
## format header
sed 's/# /No/' | sed 's/ /_/g' >$mdrOut

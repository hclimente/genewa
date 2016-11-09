#!/bin/bash
#$ -N gametesBenchmark
#$ -cwd
#$ -V
#$ -S /bin/sh
# maximum number of nodes to take
#$ -tc 200
# qsub -t 1-36000 -e populations/gametes/logs -o populations/gametes/logs scripts/gametesBenchmark.sh

#####################
#    ENVIRONMENT    #
#####################

SGE_TASK_ID=1

mkdir populations/gametes/aes

# make sure no interactive mode is active
unalias rm

#####################
#     GET PARAMS    #
#####################

param_file="populations/gametes/parameters_long.txt"

h=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -d' ' -f1`
maf=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -d' ' -f2`
N=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -d' ' -f3`
modelNo=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -d' ' -f4`
repNo=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -d' ' -f5`

# AntEpiSeeker
#######################
gametesOut=populations/gametes/pops/h"$h"_maf"$maf"_N"$N"
gametesFile="$gametesOut"_EDM-"$modelNo"/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".txt
aesFile=$gametesFile.aes
aesOut=populations/gametes/aes/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".aes.txt

sed 's/\t/,/g' $gametesFile >$aesFile

## run aes in a box
box=h"$h"_maf"$maf"_N"$N"_"$modelNo"_"$repNo"_2

mkdir $box
cd $box

sed "s,INPUT_FILE,../$aesFile," ../libs/AntEpiSeeker1.0_linux/parameters.txt |
sed "s,OUTPUT_FILE,../$aesOut," |
sed "s,NUM_SNPS,$N," >parameters.txt

../libs/AntEpiSeeker1.0_linux/AntEpiSeeker

cd ..
rm -r $box

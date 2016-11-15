#!/bin/bash
#$ -N gametes
#$ -cwd
#$ -V
#$ -S /bin/sh
# maximum number of nodes to take
#$ -tc 100

#####################
#     GET PARAMS    #
#####################

param_file="populations/gametes/parameters_short.txt"

h=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -f1`
maf=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -f2`
N=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -f3`

#####################
#    RUN GAMETES    #
#####################

gametesOut=populations/gametes/pops/h"$h"_maf"$maf"_N"$N"

# check if one model eg model #10 exists
if [ ! -d "$gametesOut"_EDM-10 ]; then
  java -jar libs/GAMETES/GAMETES_2.1.jar -M " -h 0.$h -a 0.$maf -a 0.$maf -o $gametesOut" -q 10 -p 1000 -t 100000
  java -jar libs/GAMETES/GAMETES_2.1.jar -i "$gametesOut"_Models.txt -D " -n 0.01 -x 0.5 -a $N -s 1000 -w 1000 -r 100 -o $gametesOut"
fi

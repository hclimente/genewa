#!/bin/bash
#$ -N gametes
#$ -cwd
#$ -V
#$ -S /bin/sh
# maximum number of nodes to take
#$ -tc 100

#####################
#    ENVIRONMENT    #
#####################

mkdir -p populations/gametes/mdr
mkdir populations/gametes/plink
mkdir populations/gametes/sixpac
mkdir populations/gametes/turf
mkdir populations/gametes/pops

# make sure no interactive mode is active
unalias rm

#####################
#     GET PARAMS    #
#####################

param_file="populations/gametes/parameters.txt"

h=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -f1`
maf=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -f2`
N=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -f3`

echo Analyzing population h $h maf $maf N $N

#####################
#    RUN GAMETES    #
#####################

gametesOut=populations/gametes/pops/h"$h"_maf"$maf"_N"$N"
java -jar libs/GAMETES/GAMETES_2.1.jar -M " -h 0.$h -a 0.$maf -a 0.$maf -o $gametesOut" -q 10 -p 1000 -t 100000
java -jar libs/GAMETES/GAMETES_2.1.jar -i "$gametesOut"_Models.txt -D " -n 0.01 -x 0.5 -a $N -s 1000 -w 1000 -r 100 -o $gametesOut"

########################
#  LOOK FOR EPISTASIS  #
########################

for modelNo in $(seq -f "%02g" 1 10)
do
  for repNo in $(seq -f "%03g" 1 100)
  do
    # generate PED file
    ped="$gametesOut"_EDM-"$modelNo"/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".txt.ped
    map=populations/gametes/pops/map_$N.txt
    gametesFile="$gametesOut"_EDM-"$modelNo"/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".txt
    scripts/gametes2ped.R $gametesFile $N $map $ped

    # PLINK
    plinkOut=populations/gametes/plink/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".plink.txt
    plink --noweb --epistasis --ped $ped --map $map --out $plinkOut --epi1 1 --epi2 1 --allow-no-sex

    # SIXPAC
    # sixpacOut=populations/gametes/sixpac/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".sixpac.txt
    # java -jar libs/Sixpac.jar --raw $ped --map populations/gametes/map.txt --out $sixpacOut --mode complete

    # MDR
    mdrOut=populations/gametes/mdr/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".mdr.txt
    java -jar libs/mdr_3.0.2/mdr_3.0.2.jar -min=1 -max=2 -cv=10 -table_data=true $gametesFile | \
    ## crop the Top models table
    awk '/### Top Models ###/{f=1} /FINISHED/{f=0;print} f' | \
    ## remove first and last line and blank lines
    grep -v '###' | sed '/^$/d' | sed '$ d' | \
    ## format header
    sed 's/# /No/' | sed 's/ /_/g' >$mdrOut

    # TURF
    turfOut=populations/gametes/turf/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".turf.txt

    ## run turf in a box
    box=turf_h"$h"_maf"$maf"_N"$N"
    mkdir $box
    cd $box
    if [ "$N" == "20" ]
    then
      perl ../libs/turf/TuRF-E.pl -f ../$gametesFile -o ../$turfOut -t 20
    else
      perl ../libs/turf/TuRF-E.pl -f ../$gametesFile -o ../$turfOut
    fi
    cd ..
    rm -r $box
  done
done

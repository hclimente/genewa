#!/bin/bash
#$ -N gametesBenchmark
#$ -cwd
#$ -V
#$ -S /bin/sh
# maximum number of nodes to take
#$ -tc 100
# qsub -t 1-54000 -e populations/gametes/logs -o populations/gametes/logs scripts/gametesBenchmark.sh

#####################
#    ENVIRONMENT    #
#####################

mkdir -p populations/gametes/pops
mkdir populations/gametes/mdr
mkdir populations/gametes/plink
mkdir populations/gametes/beam
mkdir populations/gametes/turf
mkdir populations/gametes/aes
mkdir populations/gametes/gwggi

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

# generate PED file
#######################
gametesOut=populations/gametes/pops/h"$h"_maf"$maf"_N"$N"
gametesFile="$gametesOut"_EDM-"$modelNo"/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".txt
ped=$gametesFile.ped
map=populations/gametes/pops/map_$N.txt

# if [ ! -s $ped ];
# then
  scripts/gametes2ped.R $gametesFile $N $map $ped
# fi

# PLINK
#######################
plinkOut=populations/gametes/plink/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".plink.txt

if [ ! -s $plinkOut.epi.cc ];
then
  plink --noweb --epistasis --ped $ped --map $map --out $plinkOut --epi1 1 --epi2 1 --allow-no-sex
fi

# MDR
#######################
mdrOut=populations/gametes/mdr/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".mdr.txt

if [ ! -s $mdrOut ];
then
  java -jar libs/mdr_3.0.2/mdr_3.0.2.jar -min=1 -max=2 -table_data=true $gametesFile | \
  ## crop the Top models table
  awk '/### Top Models ###/{f=1} /FINISHED/{f=0;print} f' | \
  ## remove first and last line and blank lines
  grep -v '###' | sed '/^$/d' | sed '$ d' | \
  ## format header
  sed 's/# /No/' | sed 's/ /_/g' >$mdrOut
fi

# TURF
#######################
turfOut=../populations/gametes/turf/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".turf.txt

## run turf in a box
box=h"$h"_maf"$maf"_N"$N"_"$modelNo"_$repNo
mkdir $box
cd $box

if [ ! -s $turfOut ];
then
  if [ "$N" == "20" ]
  then
    perl ../libs/turf/TuRF-E.pl -f ../$gametesFile -t 20
  else
    perl ../libs/turf/TuRF-E.pl -f ../$gametesFile
  fi

  perl ../libs/turf/DataConverter/ARFF2MDR.pl SNP_filtered.arff >SNP_filtered.mdr
  java -jar ../libs/mdr_3.0.2/mdr_3.0.2.jar -min=1 -max=2 -table_data=true SNP_filtered.mdr | \
  ## crop the Top models table
  awk '/### Top Models ###/{f=1} /FINISHED/{f=0;print} f' | \
  ## remove first and last line and blank lines
  grep -v '###' | sed '/^$/d' | sed '$ d' | \
  ## format header
  sed 's/# /No/' | sed 's/ /_/g' >$turfOut
fi

cd ..
rm -r $box

# BEAM
#######################
mkdir $box

beamIn=../$gametesFile.beam
beamRawOut=../populations/gametes/beam/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".beam.raw.txt
beamOut=../populations/gametes/beam/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".beam.txt

cd $box

if [ ! -s $beamOut ];
then
  ../scripts/ped2beam.R ../$ped $beamIn

  sed "s/INSERT_THIN/$N/" ../libs/beam/parameters.txt >parameters.txt
  ../libs/beam/BEAM $beamIn $beamRawOut

  ## take the relevant table
  awk '/DETECTED_ASSOCIATIONS/{f=1} /POSTERIOR_DISTRIBUTION/{f=0;print} f' $beamRawOut | \
  ## remove delimiters and empty lines
  grep -v "\[" | sed '/^$/d' | \
  ## format marker combinations
  sed 's/( //' | sed 's/ )//' | sed 's/ /,/g' >$beamOut
fi

cd ..
rm -r $box

# AntEpiSeeker
#######################
aesFile=$gametesFile.aes
aesOut=../populations/gametes/aes/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".aes.txt

sed 's/\t/,/g' $gametesFile >$aesFile

## run aes in a box
mkdir $box
cd $box

if [ ! -s $aesOut ];
then
  sed "s,INPUT_FILE,../$aesFile," ../libs/AntEpiSeeker1.0_linux/parameters.txt |
  sed "s,OUTPUT_FILE,$aesOut," |
  sed "s,NUM_SNPS,$N," >parameters.txt

  ../libs/AntEpiSeeker1.0_linux/AntEpiSeeker
fi

cd ..
rm -r $box

# GWGGI
#######################
gwggiOut=../populations/gametes/gwggi/h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo".gwggi.txt
gwggiIn=h"$h"_maf"$maf"_N"$N"_EDM-"$modelNo"_"$repNo"

mkdir $box
cd $box

if [ ! -s $gwggiOut.tamw.rst ];
then
  cut -f1-6 ../$ped >$gwggiIn.ped
  grep -v ^N ../$gametesFile | sed 's/\t.$//' >$gwggiIn.gen
  sed 's/$/A\tT/' ../$map >$gwggiIn.map

  ../libs/gwggi.unix/gwggi --file $gwggiIn --tamw --converge --ntree 4000 --clsf-sqrt --tree-depth 4 --showntree --out $gwggiOut
fi

cd ..
rm -r $box

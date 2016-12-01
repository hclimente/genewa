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

baseDir=`pwd`

# number of samples
n=$1
param_file="$baseDir/populations/gametes/parameters_long.txt"

h=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -d' ' -f1`
maf=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -d' ' -f2`
p=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -d' ' -f3`
modelNo=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -d' ' -f4`
repNo=`head -n$SGE_TASK_ID $param_file | tail -n1 | cut -d' ' -f5`

# generate PED file
#######################
gametesOut=$baseDir/populations/gametes/pops/h"$h"_maf"$maf"_n"$n"_p"$p"
gametesFile="$gametesOut"_EDM-"$modelNo"/h"$h"_maf"$maf"_n"$n"_p"$p"_EDM-"$modelNo"_"$repNo".txt
ped=$gametesFile.ped
map=$baseDir/populations/gametes/pops/map_$p.txt

# if [ ! -s $ped ];
# then
  scripts/gametes2ped.R $gametesFile $n $p $map $ped
# fi

# PLINK
#######################
plinkOut=$baseDir/populations/gametes/plink/h"$h"_maf"$maf"_n"$n"_p"$p"_EDM-"$modelNo"_"$repNo".plink.txt

if [ ! -s $plinkOut.epi.cc ];
then
  plink --noweb --epistasis --ped $ped --map $map --out $plinkOut --epi1 1 --epi2 1 --allow-no-sex
fi

# MDR
#######################
mdrOut=$baseDir/populations/gametes/mdr/h"$h"_maf"$maf"_n"$n"_p"$p"_EDM-"$modelNo"_"$repNo".mdr.txt

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
turfOut=$baseDir/populations/gametes/turf/h"$h"_maf"$maf"_n"$n"_p"$p"_EDM-"$modelNo"_"$repNo".turf.txt

## run turf in a box
box=$baseDir/tmp/h"$h"_maf"$maf"_n"$n"_p"$p"_"$modelNo"_$repNo
mkdir $box
cd $box

if [ ! -s $turfOut ];
then
  if [ "$p" == "20" ]
  then
    perl $baseDir/libs/turf/TuRF-E.pl -f $gametesFile -t 20
  else
    perl $baseDir/libs/turf/TuRF-E.pl -f $gametesFile
  fi

  perl $baseDir/libs/turf/DataConverter/ARFF2MDR.pl SNP_filtered.arff >SNP_filtered.mdr
  java -jar $baseDir/libs/mdr_3.0.2/mdr_3.0.2.jar -min=1 -max=2 -table_data=true SNP_filtered.mdr | \
  ## crop the Top models table
  awk '/### Top Models ###/{f=1} /FINISHED/{f=0;print} f' | \
  ## remove first and last line and blank lines
  grep -v '###' | sed '/^$/d' | sed '$ d' | \
  ## format header
  sed 's/# /No/' | sed 's/ /_/g' >$turfOut
fi

cd $baseDir
rm -r $box

# BEAM
#######################
mkdir $box

beamIn=$gametesFile.beam
beamRawOut=$baseDir/populations/gametes/beam/h"$h"_maf"$maf"_n"$n"_p"$p"_EDM-"$modelNo"_"$repNo".beam.raw.txt
beamOut=$baseDir/populations/gametes/beam/h"$h"_maf"$maf"_n"$n"_p"$p"_EDM-"$modelNo"_"$repNo".beam.txt

cd $box

if [ ! -s $beamOut ];
then
  $baseDir/scripts/ped2beam.R $ped $beamIn

  sed "s/INSERT_THIN/$p/" $baseDir/libs/beam/parameters.txt >parameters.txt
  $baseDir/libs/beam/BEAM $beamIn $beamRawOut

  ## take the relevant table
  awk '/DETECTED_ASSOCIATIONS/{f=1} /POSTERIOR_DISTRIBUTION/{f=0;print} f' $beamRawOut | \
  ## remove delimiters and empty lines
  grep -v "\[" | sed '/^$/d' | \
  ## format marker combinations
  sed 's/( //' | sed 's/ )//' | sed 's/ /,/g' >$beamOut
fi

cd $baseDir
rm -r $box

# AntEpiSeeker
#######################
aesFile=$gametesFile.aes
aesOut=$baseDir/populations/gametes/aes/h"$h"_maf"$maf"_n"$n"_p"$p"_EDM-"$modelNo"_"$repNo".aes.txt

sed 's/\t/,/g' $gametesFile >$aesFile

## run aes in a box
mkdir $box
cd $box

if [ ! -s $aesOut ];
then
  sed "s,INPUT_FILE,$aesFile," $baseDir/libs/AntEpiSeeker1.0_linux/parameters.txt |
  sed "s,OUTPUT_FILE,$aesOut," |
  sed "s,NUM_SNPS,$p," >parameters.txt

  $baseDir/libs/AntEpiSeeker1.0_linux/AntEpiSeeker
fi

cd $baseDir
rm -r $box

# GWGGI
#######################
gwggiOut=$baseDir/populations/gametes/gwggi/h"$h"_maf"$maf"_n"$n"_p"$p"_EDM-"$modelNo"_"$repNo".gwggi.txt
gwggiIn=h"$h"_maf"$maf"_p"$p"_EDM-"$modelNo"_"$repNo"

mkdir $box
cd $box

if [ ! -s $gwggiOut.tamw.rst ];
then
  cut -f1-6 $ped >$gwggiIn.ped
  grep -v ^N $gametesFile | sed 's/\t.$//' >$gwggiIn.gen
  sed 's/$/A\tT/' $map >$gwggiIn.map

  $baseDir/libs/gwggi.unix/gwggi --file $gwggiIn --tamw --converge --ntree 4000 --clsf-sqrt --tree-depth 4 --showntree --out $gwggiOut
fi

cd $baseDir
rm -r $box

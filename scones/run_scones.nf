#!/usr/bin/env nextflow

params.input = "$HOME/genewa/populations/gametes/pops/h*/*.txt"
gametes_populations = Channel
                          .fromPath(params.input)
                          .map { file -> tuple(file.baseName, file)}

process gametes2ped {

  input:
    set val(name), file(gametes) from gametes_populations
  output:
    file "${name}.ped" into ped
    file "${name}.map" into map

  """
  ped=${name}.ped
  map=${name}.map

  $HOME/genewa/scripts/gametes2ped.R $gametes 20 \$map \$ped

  cut -f7- \$ped | sed 's/\\t1\\t/\\tA\\t/g' | sed 's/\\t2\\t/\\tT\\t/g' > kk
  cut -f1-6 \$ped > peo

  paste peo kk >\$ped

  """

}

process ped2hdf5 {
  input:
    file ped
    file map
  output:
    file "genotype.hdf5" into hdf5In_default, hdf5In_emmax

  """
  wDir=`pwd`

  cd ~/genewa/libs/easyGWASCore
  ~/anaconda2/bin/python python/easygwascore.py data --plink2hdf5 --plink_data data/testing/scones/genotype --plink_phenotype data/testing/scones/phenotype.txt --hout \$wDir/genotype.hdf5

  """
}

process easyGWAS_default {

  input:
    file hdf5In_default
  output:
    file 'output.hdf5' into hdf5Out_default

  """
  wDir=`pwd`

  cd ~/genewa/libs/easyGWASCore
  ~/anaconda2/bin/python python/easygwascore.py gwas --hdata \$wDir/$hdf5In_default --out \$wDir

  cd \$wDir
  mv 0.hdf5 output.hdf5

  """

}

process easyGWAS_emmax {

  input:
    file hdf5In_emmax
  output:
    file 'output.hdf5' into hdf5Out_emmax

  """
  wDir=`pwd`

  cd ~/genewa/libs/easyGWASCore
  ~/anaconda2/bin/python python/easygwascore.py gwas --hdata \$wDir/$hdf5In_emmax --out \$wDir --algorithm EMMAX

  cd \$wDir
  mv 0.hdf5 output.hdf5

  """

}

process hdf52csv {

  publishDir="$HOME/genewa/scones"

  input:
    file hdf5Out_default .mix(hdf5Out_emmax)
  output:
    file 'output.csv' into csv

  """
  wDir=`pwd`

  cd ~/genewa/libs/easyGWASCore
  ~/anaconda2/bin/python python/easygwascore.py data --hfile \$wDir/$hdf5Out --csv --cout \$wDir

  cd \$wDir
  mv TEST.csv output.csv

  """

}

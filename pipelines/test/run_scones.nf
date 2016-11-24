#!/usr/bin/env nextflow

params.input = "$HOME/genewa/populations/gametes/pops/h*/*2.txt"
gametes_populations = Channel
                          .fromPath(params.input)
                          .map { file -> tuple(file.baseName, file)}

process gametes2ped {

  input:
    set val(name), file(gametes) from gametes_populations
  output:
    file "${name}.ped" into ped
    file "${name}.map" into map
    val name

  """
  ped=${name}.ped
  map=${name}.map

  $HOME/genewa/scripts/gametes2ped.R $gametes 20 \$map \$ped

  cut -f7- \$ped | sed 's/1/A/g' | sed 's/2/T/g' > gt
  cut -f1-6 \$ped > samples

  paste samples gt | sed 's/\\t/ /g' >\$ped
  sed 's/\\t/ /g' \$map >map.tmp

  mv map.tmp \$map
  """

}

process ped2hdf5 {
  input:
    file ped
    file map
    val name
  output:
    file "genotype.hdf5" into hdf5In_default, hdf5In_emmax

  """
  wDir=`pwd`

  cd ~/genewa/libs/easyGWASCore
  ~/anaconda2/bin/python python/easygwascore.py data --plink2hdf5 --plink_data \$wDir/$name --hout \$wDir/genotype.hdf5

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

process hdf5ToCsv {

  publishDir="$HOME/genewa/scones"

  input:
    file hdf5Out from hdf5Out_default .mix(hdf5Out_emmax)
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

process runScones {

  input:
    file "/bioinfo/users/hcliment/genewa/libs/easyGWASCore/data/testing/scones/genotype.ped"
    file "/bioinfo/users/hcliment/genewa/libs/easyGWASCore/data/testing/scones/genotype.map"
    file "/bioinfo/users/hcliment/genewa/libs/easyGWASCore/data/testing/scones/phenotype.txt"
    file "/bioinfo/users/hcliment/genewa/libs/easyGWASCore/data/testing/scones/network.txt"
  output:
    file "TEST.scones.out.txt" into out
    file "TEST.scones.pmatrix.txt" into pmatrix

  """

  ~/genewa/libs/easyGWASCore/bin/linux2/tools/scones genotype phenotype.txt network.txt 0.05 `pwd` additive 0

  """

}

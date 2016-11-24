#!/usr/bin/env nextflow

ped = file("$HOME/genewa/data/genesis/Genesis.ped")
map = file("$HOME/genewa/data/genesis/Genesis.map")
tab = file("$HOME/genewa/data/genesis/BIOGRID-ORGANISM-Homo_sapiens-3.4.138.tab.txt")

process create_sparse_matrix {
  input:
    file tab
  output:
    file "sparse_matrix.txt" into net

  """
  touch sparse_matrix.txt
  """

}

process get_phenotypes {
  input:
    file ped
  output:
    file "phenotype.txt" into pheno

  """
  touch phenotype.txt
  """
}

process convert2ACGT {
  input:
    file ped
    file map
  output:
    file "Genesis.ACGT.ped" into pedACGT
    file "Genesis.ACGT.map" into mapACGT

  """
  cut -f7- $ped | tr '1234' 'ACGT' > gt
  cut -f1-6 $ped > samples

  paste -d' ' samples gt >Genesis.ACGT.ped

  ln -s $map Genesis.ACGT.map
  """
}

process run_scones {

  input:
    file pedACGT
    file mapACGT
    file pheno
    file net

  """
  $HOME/genewa/libs/easyGWASCore/bin/linux2/tools/scones ${pedACGT.baseName} $pheno $net 0.05 `pwd` additive 0
  """
}

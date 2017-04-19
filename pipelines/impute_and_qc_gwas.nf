imputationScript = file("prepare_ped.nf")
nextflow_cfg = file("nextflow.config")

autosomalChromosomes = Channel.from( 1..22 )
xChromosomeRegions = [ "X_NONPAR","X_PAR1","X_PAR2" ]

process imputeAutosomal {

  input:
    file imputationScript
    file nextflow_cfg
    val chr from autosomalChromosomes
  output:
    file "chr${chr}.processed.ped" into ped_auto
    file "chr${chr}.processed.map" into map_auto
    file "warning_info.txt" into warning_info_auto

  """
  nextflow run $imputationScript --chr $chr -profile curie -resume
  """

}

process imputeX {

  input:
    file imputationScript
    file nextflow_cfg
  output:
    file "chrX_NONPAR.processed.ped" into ped_X
    file "chrX_NONPAR.processed.map" into map_X
    file "warning_info.txt" into warning_info_X

  """
  nextflow run $imputationScript --chr 23 --chr_region X_NONPAR -profile curie -resume
  """

}

process joinChromosomes {

  publishDir "~/genewa/data/genesis", mode: 'copy', overwrite: true

  input:
    file ped_chr from ped_auto .mix(ped_X) .toList()
    file map_chr from map_auto .mix(map_X) .toList()
  output:
    file "genesis.processed.ped" into ped
    file "genesis.processed.map" into map

  """
  # start by the X, as is the one we cannot use in the loop
  cp chrX_NONPAR.processed.ped genesis.processed.ped
  cp chrX_NONPAR.processed.map genesis.processed.map

  for i in `seq 1 22`
  do
    echo CHROMOSOME \$i
    paste -d' ' genesis.processed.ped <(cut -f7- chr\$i.processed.ped -d' ') > tmp.ped
    mv tmp.ped genesis.processed.ped

    cat chr\$i.processed.map >>genesis.processed.map
  done
  """

}

process studyPopulationStructure {

  publishDir "~/genewa/data/genesis", mode: 'copy', overwrite: true

  input:
    file ped
    file map
  output:
    file "genesis.processed.plot.pdf" into pca
    file "genesis.processed.lambda" into lambda

  """
  smartpca.perl -i $ped -a $map -b $ped -o genesis.processed.pca -p genesis.processed.plot -e genesis.processed.eval -l genesis.processed.pca.log
  smarteigenstrat.perl -i $ped -a $map -b $ped -p genesis.processed.pca -l genesis.processed.eigenstrat.log -o genesis.processed.chisq
  gc.perl genesis.processed.chisq genesis.processed.lambda
  """

}

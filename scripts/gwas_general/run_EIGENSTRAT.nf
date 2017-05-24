ped = file("$params.ped")
map = file("$params.map")
params.out = "."

process pca {

  container = 'mercury/eigensoft'

  input:
    file ped
    file map

  output:
    file "${ped.baseName}.pca" into pca
    file "${ped.baseName}.pca.evec" into evec

  """
  smartpca.perl -i $ped -a $map -b $ped -o ${ped.baseName}.pca -p ${ped.baseName}.plot -e ${ped.baseName}.eval -l ${ped.baseName}.pca.log
  """

}

process plotPCA {

  publishDir "$params.out", overwrite: true, mode: "copy"

  input:
    file evec
  output:
    file "pca.png"

  """
  #!/usr/bin/env Rscript

  library(ggplot2)

  ev <- read.delim("$evec", sep = "", header = F, skip = 1)

  ncols <- ncol(ev)
  colnames(ev) <- c("sample", paste0("pc", seq(1, ncols - 2)), "phenotype")

  ggplot(ev, aes(x = pc1, y = pc2, color = phenotype)) +
    geom_point() +
    theme_minimal()

  ggsave("pca.png")
  """


}

process eigenstrat {

  container = 'mercury/eigensoft'

  input:
    file ped
    file map
    file pca

  output:
    file "${ped.baseName}.chisq" into chisq

  """
  smarteigenstrat.perl -i $ped -a $map -b $ped -p $pca -l ${ped.baseName}.eigenstrat.log -o ${ped.baseName}.chisq
  """

}

process getLambda {

  publishDir "$params.out", overwrite: true, mode: "copy"
  container = 'mercury/eigensoft'

  input:
    file chisq
  output:
    file "${ped.baseName}.lambda" into lambda

  """
  gc.perl $chisq ${ped.baseName}.lambda
  """

}

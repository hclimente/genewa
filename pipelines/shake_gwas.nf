#!/usr/bin/env nextflow
// files and working directories
params.wd = "."
params.genewawd = "/Users/hclimente/projects/genewa"

wd = params.wd
genewawd = params.genewawd
srcReadPed = file("$genewawd/martiniflow/io/read_ped.nf")
srcGetNetwork = file("$genewawd/martiniflow/io/get_network.nf")

// GWAS files
ped = file("$genewawd/${params.geno}.ped")
map = file("$genewawd/${params.geno}.map")
snp2gene = file("$genewawd/$params.snp2gene")
tab = file("$genewawd/$params.tab")

// evo params
nets = ["gs","gm","gi"]
associationScore = "chi2"
modelScore = "consistency"

process readData {

  input:
    file srcReadPed
    file ped
    file map

  output:
    file "gwas.RData" into rgwas

  """
  nextflow run $srcReadPed --ped $ped --map $map -profile cluster
  """

}

rgwas.into { rgwas_getNetwork; rgwas_evo }

process getNetwork {

  input:
    file srcGetNetwork
    val net from nets
    file rgwas_getNetwork
    file snp2gene
    file tab

  output:
    file "net.RData" into rnets

  """
  nextflow run $srcGetNetwork --gwas $rgwas_getNetwork --net $net --snp2gene $snp2gene --tab $tab -profile cluster
  """

}

process run_evo {

  publishDir "$params.out", overwrite: true, mode: "copy"

  input:
    file rgwas_evo.first()
    file rnet from rnets

  output:
    file "*.RData" into analyses

    """
    #!/usr/bin/env Rscript
    library(martini)
    load("$rgwas")
    load("$rnet")

    test <- "evo.${net.baseName}"

    start.time <- Sys.time()
    cones <- search_cones(gwas, net, associationScore = "$associationScore", modelScore = "$modelScore")
    end.time <- Sys.time()
    time.taken <- end.time - start.time

    save(test, info, cones, time.taken, file = paste("cones", test, info\$id, "RData", sep = "."))
    """

}

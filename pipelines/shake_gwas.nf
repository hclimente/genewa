#!/usr/bin/env nextflow
// files and working directories
params.wd = "."
params.genewawd = "/Users/hclimente/projects/genewa"

wd = params.wd
genewawd = params.genewawd
readDataRScript = file("$genewawd/scripts/scones/r_read_gwas.nf")

// GWAS files
ped = file("$genewawd/${params.geno}.ped")
map = file("$genewawd/${params.geno}.map")
snp2gene = file("$genewawd/$params.snp2gene")
tab = file("$genewawd/$params.tab")

// evo params
associationScore = "chi2"
modelScore = "consistency"

process readData {

  input:
    each net from ["gs","gm","gi"]
    file readDataRScript
    file ped
    file map
    file snp2gene
    file tab

  output:
    file "gwas.*.RData" into rgwas
    file "net.*.RData" into rnet

  """
  nextflow run $readDataRScript --ped $ped --map $map --net $net --snp2gene $snp2gene --tab $tab -profile bigmem
  """

}

process run_evo {

  input:
    file rgwas
    file rnet

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


    save(test, cones, time.taken, file = paste(test, id, "RData", sep = "."))
    """

}

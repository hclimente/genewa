#!/usr/bin/env nextflow
// files and working directories
params.wd = "."
params.out = "."
params.genewawd = "/Users/hclimente/projects/genewa"

wd = params.wd
genewawd = params.genewawd
srcReadPed = file("$genewawd/martiniflow/io/read_ped.nf")
srcGetNetwork = file("$genewawd/martiniflow/io/get_network.nf")
srcRunEvo = file("$genewawd/martiniflow/analyze/quick_evo.nf")

// GWAS files
ped = file("$genewawd/${params.geno}.ped")
map = file("$genewawd/${params.geno}.map")
snp2gene = file("$genewawd/$params.snp2gene")
tab = file("$genewawd/$params.tab")
params.rld = "None"

// evo params
params.encoding = "additive"
params.nets = "gs,gm,gi"
params.prune = "FALSE"
params.associationScore = "chi2"
params.modelScore = "consistency"

nets = params.nets.split(",")
prune = params.prune
associationScore = params.associationScore
modelScore = params.modelScore
encoding = params.encoding

process readData {

    input:
        file srcReadPed
        file ped
        file map

    output:
        file "gwas*.RData" into rgwas

    """
    nextflow run $srcReadPed --ped $ped --map $map -profile bigmem
    """

}

rgwas.into { rgwas_getNetwork; rgwas_evo }

if (params.rld != "None") {

  rld = file("$genewawd/$params.rld")

  process getLDNetwork {

    input:
      file srcGetNetwork
      val net from nets
      file rgwas_getNetwork
      file snp2gene
      file tab
      file rld

    output:
      file "net*.RData" into rnets

    """
    nextflow run $srcGetNetwork --gwas $rgwas_getNetwork --net $net --snp2gene $snp2gene --tab $tab --prune $prune --rld $rld -profile bigmem
    """

  }

} else {

    process getNetwork {

        input:
            file srcGetNetwork
            val net from nets
            file rgwas_getNetwork
            file snp2gene
            file tab

        output:
            file "net*.RData" into rnets

        """
        nextflow run $srcGetNetwork --gwas $rgwas_getNetwork --net $net --snp2gene $snp2gene --tab $tab --prune $prune -profile bigmem
        """

    }
}

process run_evo {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file rgwas_evo
        file rnet from rnets
        file srcRunEvo

    output:
        file "cones.*.RData" into rcones
        file "cones.*.tsv" into cones

    """
    nextflow run $srcRunEvo --rgwas $rgwas_evo --rnet $rnet --associationScore $associationScore --modelScore $modelScore --encoding $encoding -profile bigmem
    """

}

if (associationScore == "chi2") {

    process ld_clump {

        publishDir "$params.out", overwrite: true, mode: "copy"

        input:
            file ped
            file map
            file cones

        output:
            file '${cones.baseName}.clumped_markers.fwf' into markers

        """
        #!/usr/bin/env Rscript

        library(tidyverse)

        read_tsv("$cones", col_types="ciiiccdld") %>%
            filter(selected) %>%
            mutate(p = pchisq(c, 1, lower.tail = FALSE)) %>%
            rename(SNP = snp, P = p) %>%
            select(SNP, P) %>%
            write_tsv("report.tsv")

        system2("plink", c("--file", "${ped.baseName}", "--filter-controls", "--clump", "report.tsv"))
        file.rename("plink.clumped", "${cones.baseName}.clumped_markers.fwf")
        """

    }

}

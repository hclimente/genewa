params.out = '.'
tpeds = Channel.fromPath('./*.tped')
tfam = file('./*.tfam')

process tped2bim {

    input:
        each tped from tpeds
        file tfam

    output:
        set file("${tped.baseName}.bed"), file("${tped.baseName}.bim"), file("${tped.baseName}.fam") into chrs

    """
    plink --tped $tped --tfam $tfam --make-bed --out ${tped.baseName}
    """

}

process bims2ped {

    publishDir "$params.out", mode: 'move', overwrite: true

    input:
        file '*' from chrs. collect()

    output:
        file 'out.ped' into ped
        file 'out.map' into map

    """
    ls | grep bed\$ | sed 's/.bed//' >chrs
    plink --merge-list chrs --recode --out out
    """

}

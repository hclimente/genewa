ped = file("${params.gen}.ped")
map = file("${params.gen}.map")
params.out = "."

process ped2bed {

  publishDir "$params.out", overwrite: true, mode: "copy"

  input:
    file ped
    file map

  output:
    file "${ped.baseName}.bed" into bed
    file "${ped.baseName}.bim" into bim
    file "${ped.baseName}.fam" into fam

  """
  plink --noweb --file ${ped.baseName} --make-bed --out ${ped.baseName}
  """

}

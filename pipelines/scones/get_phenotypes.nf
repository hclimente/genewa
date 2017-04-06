params.wd = "."
params.ped = "genotypes.ped"

ped = file("$params.ped")

process get_phenotypes {

  publishDir "$params.wd", overwrite: true, mode: "copy"

  input:
    file ped
  output:
    file "phenotype.txt" into pheno

  """
  echo FID IID TOY > phenotype.txt
  cut -d' ' -f1,2,6 $ped >> phenotype.txt
  """
}

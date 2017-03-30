ped = file("$params.ped")

process get_phenotypes {

  publishDir ".", overwrite: true

  input:
    file ped
  output:
    file "phenotype.txt" into pheno

  """
  echo FID IID TOY > phenotype.txt
  cut -d' ' -f1,2,6 $ped >> phenotype.txt
  """
}

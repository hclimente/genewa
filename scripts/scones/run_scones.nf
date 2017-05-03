// scones parameters
association_score = params.association_score
model_selection   = params.model_selection
depth             = params.depth
maf               = params.maf
lambda            = params.lambda
eta               = params.eta
outdir            = params.outdir
encoding          = params.encoding
pc                = params.pc
seed              = params.seed

// get files
ped = file("${params.gen}.ped")
map = file("${params.gen}.map")
phen = file("$params.phen")
net = file("$params.net")

process run_scones {

  publishDir "$params.out", overwrite: true, mode: "copy"

  input:
    file ped
    file map
    file phen
    file net
  output:
    file "selected_snps.txt" into selectedSnps
    file "pmatrix.txt" into pmatrix

  """
  scones2 \
    --ped ${ped.baseName} \
    --pheno $phen \
    --net $net \
    --association_score $association_score \
    --model_selection $model_selection \
    --depth $depth \
    --maf $maf \
    --lambda $lambda \
    --eta $eta \
    --outdir $outdir \
    --encoding $encoding \
    --pc $pc \
    --seed $seed

  mv *.scones.pmatrix.txt pmatrix.txt
  mv *.scones.out.txt selected_snps.txt
  """
}

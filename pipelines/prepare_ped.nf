#!/usr/bin/env nextflow

ped = file("$HOME/genewa/data/genesis/Genesis.ped")
map = file("$HOME/genewa/data/genesis/Genesis.map")
all_variants = Channel.fromPath("$HOME/genewa/data/1000GP_Phase3/1000GP_Phase3_chr*.legend.gz")
icogs_info = file("$HOME/genewa/data/genesis/icogs_snp_list.csv")

plink = "$HOME/genewa/libs/plink-1.07-x86_64/plink"
gtool = "$HOME/genewa/libs/gtool/gtool"
impute = "$HOME/genewa/libs/impute_v2.3.2_x86_64_static/impute2"

/* TODO
- Check how missing values are handled. There are supposed to be missing values in the output still...
- Check how many SNPs are lost (if they are not in the reference panel)
- gen2ped  --phenotype phenotype_1 --threshold
*/

process ped2gen {
  input:
    file ped
    file map
  output:
    file "out.gen" into gen
    file "out.sample" into sample

  """
  $gtool -P --ped $ped --map $map --og out.gen --os out.sample --binary_phenotype
  """

}

process get_snp_strand {

  input:
    file icogs_info

  output:
    file "icogs_variants" into icogs_variants

  """
  cut -d',' -f4,6,8 $icogs_info | sed 's/,/ /g' >icogs_variants
  """

}

process chunk_chromosomes {

  tag { chr_variants }

  input:
    file chr_variants from all_variants
  output:
    file "chunks" into chunks

  """
  zcat $chr_variants | cut -d' ' -f2 | sort -n >sortedPositions

  begin=`head -n2 sortedPositions | tail -n1`
  end=`tail -n1 sortedPositions`
  chr=`echo $chr_variants | cut -d'.' -f1 | sed 's/.\\+chr//'`

  for pos in `seq \$begin 6000000 \$end`
  do
    echo \$chr \$pos >>chunks
  done

  """
}

process prephase {

  clusterOptions = '-l mem=1G'
  queueSize = 1

  input:
    file gen from gen.first()
    file "chunk_parameters" from chunks.splitText()

  output:
    file "chunk_parameters" into chunk_parameters
    file "prephased.gen" into gen_prephased


  """
  chr=`cut -d' ' -f1 chunk_parameters`
  pos=`cut -d' ' -f2 chunk_parameters`

  map=$HOME/genewa/data/1000GP_Phase3/genetic_map_chr"\$chr"_combined_b37.txt

  chrXflags=''
  if [[ \$chr == *X* ]]
  then
    chrXflags="-chrX"
    if [[ \$chr != *NONPAR* ]]
    then
      chrXflags="\$chrXflags -Xpar"
    fi
  fi

  $impute \
    -prephase_g \
    # Required
    -m \$map -g $gen -int \$pos \$((\$pos + 6000000)) \
    # Basic options
    -Ne 20000 -verbose \
    # Chromosome X options
    \$chrXflags \
    # MCMC options
    ## number of samples for the reference with european origin
    -k_hap 503 \
    # Output options
    -o prephased.gen

  """
}

process impute {

  tag { chunk_parameters.text }

  clusterOptions = '-l mem=30G'
  queueSize = 10

  input:
    file chunk_parameters
    file icogs_variants
    file gen_prephased

  output:
    file "imputed" into gens_imputed

  """
  chr=`cut -d' ' -f1 $chunk_parameters`
  pos=`cut -d' ' -f2 $chunk_parameters`

  haps=$HOME/genewa/data/1000GP_Phase3/1000GP_Phase3_chr"\$chr".hap.gz
  legend=$HOME/genewa/data/1000GP_Phase3/1000GP_Phase3_chr"\$chr".legendend.gz
  map=$HOME/genewa/data/1000GP_Phase3/genetic_map_chr"\$chr"_combined_b37.txt

  grep "^\$chr " $icogs_variants | sed "s/^\$chr //" >strand_info

  chrXflags=''
  if [[ \$chr == *X* ]]
  then
    chrXflags="-chrX"
    if [[ \$chr != *NONPAR* ]]
    then
      chrXflags="\$chrXflags -Xpar"
    fi
  fi

  # impute only on the screened SNPs
  $impute \
    # Required
    -m \$map -known_haps_g $prephased -int \$pos \$((\$pos + 6000000)) \
    # Input options
    -h \$haps -l \$legend \
    # Basic options
    -Ne 20000 -verbose \
    # Strand options
    -strand_g strand_info \
    # Chromosome X options
    \$chrXflags \
    # MCMC options
    ## number of samples for the reference with european origin
    -k_hap 503 \
    # Output options
    ## only genotyped snps are included in the panel
    -os 2 -o imputed

  """

}

process gen2ped {
  input:
    file gens_imputed
  output:
    file "genesis.imputed.ped" into ped_imputed
    file "genesis.imputed.map" into map_imputed

  """
  cat out.chr*.one.phased.impute2.gen >imputed.gen
  cat out.chr*.one.phased.impute2.sample >imputed.sample

  $gtool -G --g imputed.gen --s imputed.sample --ped genesis.imputed.ped --map genesis.imputed.map --phenotype phenotype_1 --threshold 0.9
  """
}

process numeric2acgt {
  publishDir "$HOME/genewa/data/genesis"

  input:
    file ped_imputed
    file map_imputed
  output:
    file "genesis.imputed.filtered.ped" into ped_filtered
    file "genesis.imputed.filtered.map" into map_filtered

  """
  # remove indels for the moment
  grep '\\[D/I\\]\\|\\[I/D\\]' ~/genewa/data/genesis/icogs_snp_list.csv | cut -d',' -f2 >indels.txt
  $plink --file ${ped_imputed.baseName} --recode --alleleACGT --exclude indels.txt --out genesis.imputed.filtered --noweb
  """
}

#!/usr/bin/env nextflow

// parameters
params.chr = 21
params.chr_region = params.chr
chunk_size = 5000000
chr = params.chr
chr_region = params.chr_region

// GENESIS data
ped = file("$HOME/genewa/data/genesis/Genesis.ped")
map = file("$HOME/genewa/data/genesis/Genesis.map")
icogs_info = file("$HOME/genewa/data/genesis/icogs_snp_list.csv")

// 1kG data
legend = file("$HOME/genewa/data/1000GP_Phase3/1000GP_Phase3_chr${chr_region}.legend.gz")
haps = file("$HOME/genewa/data/1000GP_Phase3/1000GP_Phase3_chr${chr_region}.hap.gz")
genetic_map = file("$HOME/genewa/data/1000GP_Phase3/genetic_map_chr${chr_region}_combined_b37.txt")

process ped2gen {

  clusterOptions = '-l mem=20G'

  input:
    file ped
    file map
  output:
    file "chr${chr_region}.gen" into gen
    file "chr${chr_region}.sample" into samples_Xi, samples

  """
  awk '\$1==$chr' $map | cut -f2 >allVariants.txt

  # remove indels for the moment
  grep '\\[D/I\\]\\|\\[I/D\\]' $icogs_info | cut -d',' -f2 >indels.txt
  diff --new-line-format="" --unchanged-line-format="" <(sort allVariants.txt) <(sort indels.txt) | cat > snps.txt

  gtool -P --ped $ped --map $map --og tmp.gen --os chr"$chr_region".sample --binary_phenotype --family
  gtool -S --g tmp.gen --s tmp.sample --og chr"$chr_region".gen --inclusion snps.txt
  """

}

process get_snp_strand {

  input:
    file icogs_info

  output:
    file "strand_info" into strand_info

  """
  chr_clean=`echo $chr_region | sed 's/_.\\+//'`
  cut -d',' -f4,6,8 $icogs_info | sed 's/,/ /g' | grep "^\$chr_clean " | sed "s/^\$chr_clean //" >strand_info
  """

}

process chunk_chromosome {

  input:
    file legend
    file map
  output:
    file "chunks" into chunks

  """
  zcat $legend | cut -d' ' -f2 | sort -n >sortedPositions

  first_snp_pos=`head -n2 sortedPositions | tail -n1`
  last_snp_pos=`tail -n1 sortedPositions`
  chr=`echo $legend | cut -d'.' -f1 | sed 's/.\\+chr//'`
  chr_clean=`echo \$chr | sed 's/_.\\+//'`

  for pos in `seq \$first_snp_pos $chunk_size \$last_snp_pos`
  do
    begin=\$pos
    end=\$((\$pos + $chunk_size - 1))

    any_snp=`awk -v b="\$begin" -v e="\$end" '\$1==$chr && \$4 >= b && \$4 <= e' $map`

    if [[ ! -z \$any_snp ]]
    then
      if [ "\$end" -gt "\$last_snp_pos" ]
      then
        end=\$last_snp_pos
      fi
      echo \$chr \$begin \$end >>chunks
    fi
  done
  """
}

process impute {

  clusterOptions = '-l mem=30G'
  queueSize = 15

  input:
    file sample from samples_Xi.first()
    file strand_info from strand_info.first()
    file gen from gen.first()
    file "chunk_info" from chunks.splitText()

  output:
    file "chr${chr_region}.*.imputed.gen" into imputed_gens
    file "chr${chr_region}.*.imputed.gen_warnings" into warnings

  """
  chr=$chr_region
  begin=`cut -d' ' -f2 chunk_info`
  end=`cut -d' ' -f3 chunk_info`

  chrXflags=''
  if [[ \$chr == *X* ]]
  then
    chrXflags="-sample_g $sample -chrX"
    if [[ \$chr != *NONPAR* ]]
    then
      chrXflags="\$chrXflags -Xpar"
    fi
  fi

  # impute only on the screened SNPs
  ## only genotyped snps are included in the panel
  ## use 503 samples for the reference with european origin
  impute2 \
    -g $gen \
    -m $genetic_map \
    -int \$begin \$end \
    -h $haps \
    -l $legend \
    -Ne 20000 \
    -verbose \
    -strand_g $strand_info \
    \$chrXflags \
    -k_hap 503 \
    -os 2 -o chr$chr_region.\$begin.\$end.imputed.gen

  """

}

process gen2ped {
  publishDir ".", mode: 'move', overwrite: true

  input:
    file imputed_gen from imputed_gens.toList()
    file("chr${chr_region}.imputed.sample") from samples
  output:
    file "chr${chr_region}.processed.ped" into ped_imputed
    file "chr${chr_region}.processed.map" into map_imputed

  """
  cat ${imputed_gen} >chr"$chr_region".imputed.gen

  gtool -G --g chr"$chr_region".imputed.gen --s chr"$chr_region".imputed.sample --ped chr"$chr_region".processed.ped --map chr"$chr_region".processed.map --phenotype phenotype --threshold 0.95 --snp

  # recover family information
  sed 's/\\t/ /g' chr"$chr_region".processed.ped | sed 's/_[^ ]\\+//' | sed 's/[^ ]\\+_//' >tmp && mv tmp chr"$chr_region".processed.ped
  """
}

process extractWarningPositions {

  input:
    file all_warns from warnings.toList()

  output:
    file "warning_positions.txt" into warning_positions

  """
  sed 's/[^0-9]\\+//' *.imputed.gen_warnings | sed 's/in Panel 2.\\+//' >warning_positions.txt
  """

}

process getWarningInformation {
  publishDir ".", mode: 'move', overwrite: true

  input:
    file warning_positions

  output:
    file "warning_info.txt" into warning_information

  """
  #!/usr/bin/env Rscript

  library(readr)
  library(magrittr)
  library(dplyr)

  warns <- read_tsv("$warning_positions", col_names=F)

  icogs <- read_csv("~/genewa/data/genesis/icogs_snp_list.csv") %>%
    filter(Chromosome == $chr & Build37_Position %in% warns\$X1) %>%
    select(Illumina_SNP_Name,Alternative_SNP_Name,Chromosome,Build37_Position,Strand,ForwardAlleles,DesignAlleles) %>%
    write_tsv("warning_info.txt")


  """

}

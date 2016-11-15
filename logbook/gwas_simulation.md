# Choice of simulation tool

NIH has a [nice comparative](https://popmodels.cancercontrol.cancer.gov/gsr/search/) of genome simulation tools. We filtered by the following criteria:

* Select those that have a GSR Certification, which has is interesting from both documentation and support perspectives.
* Tested in Linux
* Type of Simulated Data is Genotype at Genetic Markers.
* Phenotype is determined by Gene-Gene Interaction.

Two well documented tools emerge [simuPOP](http://simupop.sourceforge.net) and [SLiM](https://messerlab.org/slim/). While the former one has a Python interface, the latter has a ad hoc scripting language. Hence, the first one looks more interesting. The combo [simuPOP](http://simupop.sourceforge.net) + [simuGWAS](http://simupop.sourceforge.net/Cookbook/SimuGWAS) + [GENS2](https://sourceforge.net/projects/gensim/) does not seem to be a standard in the literature (6 papers use it). Indeed, [HAPGEN 2](https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html) seem to be the preferred option for people developing GWAS methods (66 articles in the last two years use HAPGEN2 for simulated GWAS studies). Additionally, [one published paper](http://link.springer.com/article/10.1007/s00702-014-1341-9) uses it with a purpose similar to ours, which we can use as guideline.

# HAPGEN 2

HAPGEN requires 3 files with genomic information, all of which can [easily be obtained](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#reference) from the authors. Those are:

* **Known haplotypes**, with one row per SNP and one column per haplotype. Every haplotype file needs a corresponding legend file (see below), and all alleles must be coded as 0 or 1 -- no other values are allowed.
* A legend file for the **SNP markers**. This file should have 4 columns with one line for each SNP:
  1. ID for each SNP i.e. rs id of the marker.
  2. The base pair position of each SNP.
  3. Base represented by 0.
  4. Base represented by 1.
* A file containing the fine-scale **recombination rate** across the region. This file should have 3 columns with one line for each SNP:
  1. physical location,
  2. rate in cM/Mb to the right of the marker
  3. cumulative rate in cM to the left of the marker.

Then, it requires some user-specified parameters, such as the number of cases and controls to be simulated, relative risks associated with disease SNPs, a subset of SNPs to output, etc.

We downloaded the 1000G phase 3 data, and simulated 1000 cases and 1000 controls for the chromosome 20. We arbitrarily picked the biallellic SNP rs13039134 (pos 92366), and tagged the minor allele G as pathogenic (heterozygous disease risk = 1.5; homozygous disease risk = 2.25). We ran the simulation using the following command:

```bash

echo "libs/hapgen2/hapgen2 -m data/genetic_map_chr20_combined_b37.txt -h data/1000GP_Phase3_chr20.hap -l data/1000GP_Phase3_chr20.legend -o populations/chr20 -n 1000 1000 -dl 92366 1 1.5 2.25" | qsub -cwd -S /bin/bash -V -o o.chr20.100.txt -e e.chr20.100.txt -N hapgen.chr20

```
The output consist of the following files:

* populations/chr20.[cases|controls].haps: simulated haplotype data.
* populations/chr20.legend: information about the SNPs in the .haps files.
* populations/chr20.[cases|controls].gen: simualted genotype data
* populations/chr20.[cases|controls].sample: a sample file for the simulated genotype data.
* populations/chr20.[cases|controls].tags (if applicable): the genotype data limited to the subset of SNPs specified by the file supplied to the -t flag.

# SimulatePhenotypes

*to do*

# GAMETES

We discovered [GAMETES](https://sourceforge.net/projects/gametes/?source=navbar) ([paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3605108/)) when reviewing the bibliography of [Moore](https://scholar.google.fr/citations?user=mE1Te78AAAAJ&hl=en&oi=ao). It generates pure, random and strict epistatic models. GAMETES generates a number of random model arquitectures under the specified constraints (MAF, Heritability...); the size of this sampe is specified in the parameter *Quantile population size*. Then, those models are scored according to the metric specified under *Quantile* (EDM or Odds ratio). The higher those values are, the easier the underlying epistasis is to detect. The quantiles specified in *Quantile count* are picked from the distribution of scores. Then it generates sample datasets for each of the generated models. Those datasets consist of a number of cases and controls, and a specified number of cohorts.

We generate a set of models with different heritabilities and MAFs. We allow prevalence be chosen by the program because, as stated in the paper, it doesn't make much difference. Then, for each of them, we generate 100 replicates of 1000 cases and 1000 controls. We will do it for a variable number of SNPs: 20, 100, 1000.

```bash
rm populations/gametes/parameters_short.txt

for h in "005" "01" "025" "05" "1" "2"
do
	for maf in "05" "2" "4"
	do
    for N in "1000" "100" "20"
    do
      echo $h $maf $N >>populations/gametes/parameters_short.txt
    done
  done
done

qsub -t 1-54 -e populations/gametes/logs -o populations/gametes/logs scripts/gametes.sh
```

Example of the kind of data generated by GAMETES:

```r
library(magrittr)
library(dplyr)

# pool the 100 replicates
for (repNo in formatC(1:100, width = 3, flag = "0")){
  d <- read.delim(paste0("populations/gametes/pops/h2_maf2_N100_EDM-01/h2_maf2_N100_EDM-01_",repNo,".txt"))

  c0.d <- d %>% filter(Class == 0) %>% select(M0P0,M0P1) %>% table
  c1.d <- d %>% filter(Class == 1) %>% select(M0P0,M0P1) %>% table

  if (exists("c0")){
    c0 <- c0 + c0.d
    c1 <- c1 + c1.d
  } else {
    c0 <- c0.d
    c1 <- c1.d
  }

}

#ratio between both contingency tables. They both include 100*1000 individuals
c1/c0

# M0P1
# M0P0          0          1          2
# 0  0.8341160  1.5908712  0.4084012
# 1  1.6137866  0.2437159  6.7852619
# 2  0.3246154 13.4517134  0.9432624

```
The contingency matrixes have this "diamond" unbalance, meaning that only a combination of heterozygosity in one alele and homozygosity in the other one will increase the chances of developing the disease. This is because GAMETES generates pure strict epistatic models, where no marginal effects are to be observed.

# Open questions

* HAPGEN2: do we need to impute the haplotypes to get SNPs on all the genes of the haplotype, instead of the haplotypes? Else it will be difficult to build any kind of biological network.

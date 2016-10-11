# gwas_simulation

# Choice of simulation tool

NIH has a [nice comparative](https://popmodels.cancercontrol.cancer.gov/gsr/search/) of genome simulation tools. We filtered by the following criteria:

* Select those that have a GSR Certification, which has is interesting from both documentation and support perspectives.
* Tested in Linux
* Type of Simulated Data is Genotype at Genetic Markers.
* Phenotype is determined by Gene-Gene Interaction.

Two well documented tools emerge [simuPOP](http://simupop.sourceforge.net) and [SLiM](https://messerlab.org/slim/). While the former one has a Python interface, the latter has a ad hoc scripting language. Hence, the first one looks more interesting. The combo [simuPOP](http://simupop.sourceforge.net) + [simuGWAS](http://simupop.sourceforge.net/Cookbook/SimuGWAS) + [GENS2](https://sourceforge.net/projects/gensim/) does not seem to be a standard in the literature (6 papers use it). Indeed, [HAPGEN 2](https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html) seem to be the preferred option for people developing GWAS methods (66 articles in the last two years use HAPGEN2 for simulated GWAS studies). Additionally, [one published paper](http://link.springer.com/article/10.1007/s00702-014-1341-9) uses it with a purpose similar to ours, which we can use as guideline.

## HAPGEN 2 + SimulatePhenotypes

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
  
## 


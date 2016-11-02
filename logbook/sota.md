We want to compile a list of existing tools to detect epistasis. Firstly, we want to test the performance of existing tools. Not only under the lenses of accuracy and type I errors, but also of **stability** ie how robust are the results. In order to do this, we will sample our simulated data and compare the results (eg correlation between the p-values among the different runs). Secondly, we want to test how different simulation settings affect the outcome. For example, what are the lower boundaries of incidence to detect epistasis.

# List of tools for disease (case/control)

Curated collection of tools compiled from [several](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4564769/) [papers](http://www.nature.com/nrg/journal/v15/n11/full/nrg3747.html) and [OMICtools](https://omictools.com/epistasis-detection-category). Only tools still accessible are kept.

## Exhaustive search

* Regression-based methods: Methods based on logistic regression. Time consuming and low powered.

  * **[PLINK epistasis-module](http://pngu.mgh.harvard.edu/~purcell/plink/epi.shtml)**
  * [BOOST](http://bioinformatics.ust.hk/BOOST.html). *Download site is recognised as malware by our proxy. His email server rejects our emails.*

* LD- and haplotype-based methods: They check if the joint co-occurrence is different between cases and controls. Faster and higher powered. Works well for *unlinked* loci (eg separate chromosomes) in *rare* diseases.

  * [SIXPAC](http://www.cs.columbia.edu/~snehitp/sixpac/). *Unclear problem with input format, and lack of documentation.*

## Non-exhaustive search

* Data-filtering methods: hypothesis-driven approaches aim to select a subset of SNPs for interaction tests on the basis of existing biological knowledge (for example, databases of pathways and protein–protein interactions), statistical features (for example, marginal effects and SNP genotype frequencies) or fast algorithms.

  * **[ReliefF series](https://code.google.com/archive/p/ensemble-of-filters/)**
  * [SNPHarvester](http://bioinformatics.ust.hk/SNPHarvester.html)
  * [EDCF](http://www.cs.ucr.edu/~minzhux/EDCF.zip)

* Artificial intelligence algorithms: use classifiers for data reduction and/or feature selection to reduce both the computational burden and the statistical burden of an exhaustive search.

  * **[MDR](https://sourceforge.net/projects/mdr/)**
  * Random forest
    * [SNPInterForest](https://gwas.biosciencedbc.jp/SNPInterForest/index.html)
    * [GWGGI](https://github.com/changshuaiwei/gwggi)
  * Bayesian network-based
    * **[BEAM](http://sites.stat.psu.edu/~yuzhang/)**: reference among Bayesian-based classifiers. Needed to compile the Linux version.
    * EpiBN: apparently outperforms BEAM. However, we couldn't find any download links on the paper or on the internet.
    * DASSO-MB: also outperforms BEAM and is nowhere to be found.
  * Ant colony optimization
    * [MACOED](http://www.csbio.sjtu.edu.cn/bioinf/MACOED/): available for Windows.
    * [AntEpiSeeker](https://github.com/wyp1125/AntEpiSeeker)
  * [SNPruler](http://bioinformatics.ust.hk/Software.html)
  * [ranger](https://github.com/imbs-hl/ranger)

# Evaluation of tools

## PLINK

![PLINK gametes accuracy](results/sota_benchmark/plink.gametes.accuracy.png)

![PLINK gametes sensitivity](results/sota_benchmark/plink.gametes.sensitivity.png)

![PLINK gametes specificity](results/sota_benchmark/plink.gametes.specificity.png)

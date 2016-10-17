We want to compile a list of existing tools to detect epistasis. Firstly, we want to test the performance of existing tools. Not only under the lenses of accuracy and type I errors, but also of **stability** ie how robust are the results. In order to do this, we will sample our simulated data and compare the results (eg correlation between the p-values among the different runs). Secondly, we want to test how different simulation settings affect the outcome. For example, what are the lower boundaries of incidence to detect epistasis.

# List of tools for disease (case/control)

Curated collection of tools compiled from the [literature](http://www.nature.com/nrg/journal/v15/n11/full/nrg3747.html) and [OMICtools](https://omictools.com/epistasis-detection-category). Only tools still accessible are kept.

## Regression-based methods

Methods based on logistic regression. Time consuming and low powered.

* [PLINK epistasis-module](http://pngu.mgh.harvard.edu/~purcell/plink/epi.shtml)
* [FastEpistasis](http://www.vital-it.ch/software/FastEpistasis)
* [BOOST](http://bioinformatics.ust.hk/BOOST.html)
* [SNP-SNP Interactions](https://sites.ualberta.ca/~yyasui/software.html#snpGWAS)
* [BiForce](http://bioinfo.utu.fi/biforcetoolbox)
* [PIAM](http://www.ihs.ac.cn/xykong/PIAM.zip)

## LD- and haplotype-based methods

They check if the joint co-occurrence is different between cases and controls. Faster and higher powered. Works well for *unlinked* loci (eg separate chromosomes) in *rare* diseases.

* [SIXPAC](http://www.cs.columbia.edu/~snehitp/sixpac/)
* [iLOCi](http://www4a.biotec.or.th/GI/tools/iloci)

## Bayesian methods

Classify pairs of unlinked loci in one of three non-overlapping groups (unassociated, associated by marginal effects and associated by joint effects) according to their posterior probabilities.

* [BEAM](http://sites.stat.psu.edu/~yuzhang/)
* [GenomeMatrix](https://sph.uth.edu/hgc/faculty/xiong/software-B.html)
* [IndOR](http://emily.perso.math.cnrs.fr/IndOR/IndOR/IndOR.html)
* [HAPAL](https://www.stt.msu.edu/~cui/software.html)

## Data-filtering methods

Hypothesis-driven approaches aim to select a subset of SNPs for interaction tests on the basis of existing biological knowledge (for example, databases of pathways and protein–protein interactions), statistical features (for example, marginal effects and SNP genotype frequencies or fast algorithms.

* [SNPHarvester](http://bioinformatics.ust.hk/SNPHarvester.html)
* [EDCF](http://www.cs.ucr.edu/~minzhux/EDCF.zip)
* [Relief series](https://code.google.com/archive/p/ensemble-of-filters/)

## Artificial intelligence algorithms

Use classifiers for data reduction and/or feature selection to reduce both the computational burden and the statistical burden of an exhaustive search.

* [SNPruler](http://bioinformatics.ust.hk/Software.html)
* [SNPInterForest](https://gwas.biosciencedbc.jp/SNPInterForest/index.html)
* [ranger](https://github.com/imbs-hl/ranger)

# Evaluation of tools

*to do*

In this repository, we evaluate six network-based GWAS tools on the [GENESIS dataset](http://bmccancer.biomedcentral.com/articles/10.1186/s12885-015-2028-9). It containts the accompanying scripts, results, and laboratory notebooks of the following article:

> **Biological networks and GWAS: comparing and combining network methods to understand the genetics of familial breast cancer susceptibility in the GENESIS study.** Héctor Climente-González, Christine Lonjou, Fabienne Lesueur, GENESIS Study collaborators, Dominique Stoppa-Lyonnet, Nadine Andrieu, Chloé-Agathe Azencott,  GENESIS study group. *bioRxiv* 2020.05.04.076661; doi: https://doi.org/10.1101/2020.05.04.076661

This repository is organized into four main subfolders:

- [scripts/](scripts): contains the mai pieces code using to run the experiments.
- [results/](results): contains the results of the experiments. Generally, every subfolder contains a script `run.sh` which, upon running, produced all the files within that directory.
- [doc/](doc): contains the Jupyter notebooks in which the results are analyzed (see [below](#analyses)).

## Analyses

* [Conventional GWAS](doc/gwas.ipynb), including [SNP-](results/conventional_gwas/univariate_models.no_covars.tsv) and [gene-level](results/preprocessing/scored_genes.vegas.txt) summary statistics.
* [Recovering known biomarkers](doc/genesis_bcac_comparison.ipynb)
* [Method benchmark](doc/benchmark.ipynb)
* [SConES results](doc/scones.ipynb)
* [Consensus network](doc/consensus.ipynb)
* [Stable consensus](doc/stability.ipynb)
* [Genes excluded from the analysis](doc/bad_genes.ipynb)

## Manuscript

* [Manuscript](manuscript/plos/review_1.pdf)
* [Figures](manuscript/plos/figures.ipynb)
* [Supplementary figures](manuscript/plos/supp_figures.ipynb)

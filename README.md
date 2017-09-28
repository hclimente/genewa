# genewa project

## Logbook

1. GENESIS dataset

    * [Description](logbook/genesis_data.md)
    * [Data overview + SConES](logbook/genesis_exploration_1.ipynb)

2. GIN
    * Playing with the simplest example: [1](logbook/simplest_example_1.ipynb) and [2](logbook/simplest_example_2.ipynb).
    * Playing with little green men: [1](logbook/little_green_men_example_1.ipynb).
    * Preliminary analysis on GENESIS for [JED1A](logbook/genesis_gin_skat_aicc_1.ipynb).

3. EVO: rewriting easyGWAS for extension and mantainability.

    * [Benchmark](results/heritability/benchmark.ipynb) at h<sup>2</sup> = 1 to pick the best experimental conditions, which turned out to be chi2 and consistency.
    * [Evo on GENESIS](results/evo/evo_analysis.ipynb)

4. Epistasis
    * [Data simulation](logbook/gwas_simulation.md)
    * [State of the art](logbook/sota.md)


## Requirements

Required files:

* libs
  * [hapgen2](https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html#Download_and_Compilation)
  * [SimulatePhenotypes](https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/download/R_scripts/SimulatePhenotypes_1.0.tar.gz)
* data
  * 1000GP_Phase3_chr20.hap
  * 1000GP_Phase3_chr20.legend
  * genetic_map_chr20_combined_b37.txt
* populations

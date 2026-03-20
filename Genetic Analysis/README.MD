# Genetic Analysis

## Requirements
* Download [STRUCTURE](https://web.stanford.edu/group/pritchardlab/structure.html) to `./STRUCTURE/structure`
* Download [ADMIXTURE](https://dalexander.github.io/admixture/) to `./ADMIXTURE/admixture`
* Download [PLINK](https://www.cog-genomics.org/plink/) to `./PLINK/plink`

## Data
* Data: `Report_DAc23-8280_SNP_2.csv`
* Metadata: `ind_metrics_site.csv`

## `SNP data analysis.rmd`
Code for SNP quality control, filtering, visualisation and analysis to infer population structure. See comments within code for each analysis step. 

## `ADMIXTURE/ADMIXTURE.sh`
Referenced in `SNP data analysis.rmd`
* Runs ADMIXTURE in a loop with multiple values of K and multiple runs per L
* Checks cross-validation error to find the best value of K
* Visualises using PONG
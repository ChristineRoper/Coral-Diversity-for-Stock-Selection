## 1. Fluorcam data folder

Extracts data from Fluorcam files

## 2. Fit parameters to thermal curves

Run ed50.r to create Fq´/Fm´ log logistic curves to extract ED50 values

   - Plots and csv files will be output to the `outputs` folder.

## 3. Statistical Analysis

   - Mean ED50 compared across sub-populations using a Kruskal-Wallis test by ranks and a Dunn’s Multiple Comparison post-hoc test (with Benjamini-Hochberg correction). Differences in sub-population variances were examined using a Levene’s test, followed by Tukey’s HSD post-hoc test.

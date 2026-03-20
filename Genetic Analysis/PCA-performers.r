# Install
## One-time install
# install.packages("dartR")
# install.packages("BiocManager")
# install.packages("devtools")
# BiocManager::install(c("SNPRelate", "qvalue"))
# gl.install.vanilla.dartR()

library(readr)
library(dplyr)
library(dartR)

# Import and check genlight object

gl <- gl.read.dart(
  filename = "Genetic Analysis/Report_DAc23-8280_SNP_2.csv",
  ind.metafile = "Genetic Analysis/ind_metrics_site.csv"
)

# filter out Low Isles slope pale as we are not using the matching phenotypic data, filter mangrove coral as 1 individual is not valid for a population analysis
gl <- gl.drop.pop(gl, pop.list = "Low Isles slope pale")
gl <- gl.drop.pop(gl, pop.list = "Low Isles mangrove")
# filer out individuals that have been identified as outliers in the PCA RB F1, MO S30, MO S23, MO S24, RB F14, MO S26 - these individuals have also been removed from phenotype analysis
gl <- gl.drop.ind(gl, ind.list = c("89", "146", "139", "140", "102", "142"))

# to check genlight object is compliant
gl <- gl.compliance.check(gl)

# save(gl, file="gl_dartR.rdata")


# Filtering data
## 1. FILTER SECONDARIES
# gl.report.secondaries(gl) # report
gl_f <- gl.filter.secondaries(gl, verbose = 3) # filtering
# this is standard, no threshold setc

## 2. FILTER REPRODUCIBILITY
# gl.report.reproducibility(gl_f) # report
gl_f <- gl.filter.reproducibility(gl_f, threshold = 0.99, verbose = 3) # filtering threshold 0.99

## 3. FILTER CALL RATE BY LOCI
# gl.report.callrate(gl_f, method = "loc") # report
gl_f <- gl.filter.callrate(gl_f, method = "loc", threshold = 0.80, verbose = 3) # filtering
# Note: this has been changed to 0.8 due to significant loss of data at higher thresholds
# This retains a good number of loci, 0.93 also seems appropriate, 0.95 seems to strict (only 50 loci remain)

## 4. FILTER CALL RATE BY INDIVIDUALS
# gl.report.callrate(gl_f, method = "ind") # report
gl_f <- gl.filter.callrate(gl_f, method = "ind", threshold = 0.80, verbose = 3) # filtering

## 5. FILTER READ DEPTH TO REMOVE LOW COVERAGE LOCI
# gl.report.rdepth(gl_f) # report
gl_f <- gl.filter.rdepth(gl_f, lower = 5, upper = 85, verbose = 3) # filtering
# selecting a threshold of 10% either side (data between the 10th and 90th percentile was kept)

## 6. FILTER MONOMORPHIC LOCI
# gl.report.monomorphs(gl_f) # report
gl_f <- gl.filter.monomorphs(gl_f, verbose = 3) # filtering
# Filtered prior to MAF as it was generating a warning about monomorphoc loci

## 7. FILTER Minor Allele Frequencies (MAF)
# gl.report.maf(gl_f) # report
gl_f <- gl.filter.maf(gl_f, threshold = 0.01, verbose = 3) # filtering

# filtering complete.

gl_f <- gl.recalc.metrics(gl_f) # re-calculate locus metadata now that filtering is complete


# Test for clones - none under the threshold 0.45
# 80% & 70% loci coverage
gl_grm <- gl.grm(gl_f, verbose = 3)

# This function calculates the mean probability of identity by state (IBS) across loci that would result from all the possible crosses of the individuals analyzed.
# IBD is calculated by an additive relationship matrix approach developed by Endelman and Jannink (2012)

# Represents a genomic relationship matrix (GRM) as a network
gl_grm_net <- gl.grm.network(gl_grm, gl_f, method = "mds", link.size = 2, relatedness_factor = 0.45) # decided on value between sibling/half sibling


# Impute missing data using the nearest neighbour function - needed for some downstream analyses (e.g. Fst)
gl_f2 <- gl.impute(gl_f, method = "neighbour", verbose = 5)
# re-calculating the locus metrics after manipulation
gl_f2 <- gl.recalc.metrics(gl_f2)

# 2244 / (nInd(gl_f) * nLoc(gl_f)) # Percentage imputed


# Check for Symbiodiniaceae contamination - no symbiont sequences aligned with filtered data indicating no contamination
# Download Symbiodiniaceae genomes here: sampgr.org.cn/index.php/download
# Download & install BLAST (needed to run analysis) https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
# select this file: ncbi-blast-2.14.1+.dmg
# Run the BLAST function for each genome
gl.blast(gl_f, "Symbiontgenomes/symC_scaffold_40.fasta", verbose = 5)
gl.blast(gl_f, "Symbiontgenomes/Symbiodinium_microadriacticum_genome.scaffold.fasta", verbose = 5)
gl.blast(gl_f, "Symbiontgenomes/symA3_scaffold_37.fasta", verbose = 5)
gl.blast(gl_f, "Symbiontgenomes/Fugacium_kawagutii_V3_genome_scaffold.fasta", verbose = 5)
gl.blast(gl_f, "Symbiontgenomes/Fugacium_kawagutii_V2_genome_Scaffolds.fasta", verbose = 5)
gl.blast(gl_f, "Symbiontgenomes/Fugacium_kawagutii_V1_genome_scaffold.fasta", verbose = 5)
gl.blast(gl_f, "Symbiontgenomes/Cladocopium_goreaui_Genome.Scaffolds.fasta", verbose = 5)
gl.blast(gl_f, "Symbiontgenomes/Breviolum_minutum.v1.0.genome.fa", verbose = 5)




# Load phenotypic performance data
sample_to_numeric_ids <- read_csv("Genetic Analysis/ind_metrics_site.csv") %>%
  select(id, "Sample ID") %>%
  rename(SampleID = "Sample ID")


ed50s <- read_csv("Phenotype Analysis/3. Thermal Curves/outputs/FqFm_T3/ed50.csv") %>%
  left_join(sample_to_numeric_ids, by = "SampleID")

top20 <- quantile(ed50s$ed50, 0.8)
bottom20 <- quantile(ed50s$ed50, 0.2)

ed50s = ed50s %>%
  mutate(id = as.numeric(id), pop = case_when(
    ed50 > top20 ~ "top",
    ed50 < bottom20 ~ "bottom",
    TRUE ~ "middle"
  )) %>%
  select(id, pop) %>%
  filter(!is.na(id)) 

pop(gl_f) <- ed50s$pop[match(indNames(gl_f), ed50s$id)]

# PCA
library(adegenet)
source("../default_theme.r")

pcoa <- gl.pcoa(gl_f, verbose = 3)

# gl.pcoa.plot won't return the plot for us to work with
# We can tell it to save it in a temporary directory though
gl.pcoa.plot(pcoa, gl_f, pop.labels = "pop", ellipse = FALSE, save2tmp = TRUE)
# Now we need to get the list of repots it has saved and find the last plot
gl.list.reports()
# put the number of that plot here so we can finally get the plot and work with it
plot <- gl.print.reports(1)

# Remove the horizontal and vertical lines
plot$layers[[4]] <- NULL
plot$layers[[3]] <- NULL

# output the plot
plot +
  # stat_ellipse(type = "norm", level = 0.95) +
  xlim(-5, 10) + ylim(-10, 10) +
  default_theme

ggsave("Genetic Analysis/plots/pca-performers.png", width = 6, height = 5)

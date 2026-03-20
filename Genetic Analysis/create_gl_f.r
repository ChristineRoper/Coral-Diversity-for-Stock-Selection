# Load Report_DAc23-8280_SNP_2.csv and ind_metrics_site.csv
# Filter out low isles slope pale, low isles mangrove, and individuals with missing data
# Filter out secondaries, reproducibility, call rate, rdepth, monomorphs, and MAF
# Save filtered data to gl_f.RData for use in other scripts

library(dartR)

# Load the genlight object (same filtering as main analysis)
gl <- gl.read.dart(
    filename = "Genetic Analysis/Report_DAc23-8280_SNP_2.csv",
    ind.metafile = "Genetic Analysis/ind_metrics_site.csv"
)

# Apply same filtering as main analysis
gl <- gl.drop.pop(gl, pop.list = "Low Isles slope pale")
gl <- gl.drop.pop(gl, pop.list = "Low Isles mangrove")
gl <- gl.drop.ind(gl, ind.list = c("89", "146", "139", "140", "102", "142"))
gl <- gl.compliance.check(gl)

save(gl, file = "gl.RData")

# Filter SNPs
gl_f <- gl.filter.secondaries(gl, verbose = 3)
gl_f <- gl.filter.reproducibility(gl_f, threshold = 0.99, verbose = 3)
gl_f <- gl.filter.callrate(gl_f, method = "loc", threshold = 0.80, verbose = 3)
gl_f <- gl.filter.callrate(gl_f, method = "ind", threshold = 0.80, verbose = 3)
gl_f <- gl.filter.rdepth(gl_f, lower = 5, upper = 85, verbose = 3)
gl_f <- gl.filter.monomorphs(gl_f, verbose = 3)
gl_f <- gl.filter.maf(gl_f, threshold = 0.01, verbose = 3)
gl_f <- gl.recalc.metrics(gl_f)

save(gl_f, file = "gl_f.RData")

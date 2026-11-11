#!/usr/bin/env Rscript
# Draw the SIMEX subsamples from the complete cohort.
#
# For each sampling fraction lambda in the grid, N_REPS independent subsamples
# are drawn without replacement. The number of isolates retained is
# lambda / OBS_LAMBDA of the sequenced cohort, so that lambda is expressed on
# the scale of the notified case burden rather than of the sequenced set.
#
# Sites invariant within a subsample are removed from that subsample before it
# is written, so that each replicate is reconstructed from a genuinely
# variable-site alignment under the ascertainment-bias correction.
#
# Writes output/fasta/<run>.fa, output/fasta/<run>_taxa.txt and
# output/manifest.csv. The complete cohort is included as the run "observed".
suppressPackageStartupMessages({
  library(ape); library(readr); library(dplyr)
})
script_dir <- function() {
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "."
}
source(file.path(script_dir(), "common.R"))

set.seed(SEED_SUBSAMPLE)
dir.create(file.path(OUTPUT, "fasta"), showWarnings = FALSE, recursive = TRUE)
for (d in c("trees", "trees_rooted", "snppar_runs", "logs", "results")) {
  dir.create(file.path(OUTPUT, d), showWarnings = FALSE, recursive = TRUE)
}

alignment <- read_alignment()
cat("alignment:", nrow(alignment), "taxa x", ncol(alignment), "sites\n")
cat("observed sampling fraction:", OBS_LAMBDA, "\n")
cat("sampling fractions:", paste(round(LAMBDAS, 3), collapse = ", "), "\n")
cat("replicates per fraction:", N_REPS, "\n\n")

manifest <- list()
write_run <- function(run_id, lambda, keep_frac, replicate, taxa) {
  m <- drop_invariant(alignment[taxa, , drop = FALSE])
  write_fasta(m, file.path(OUTPUT, "fasta", paste0(run_id, ".fa")))
  writeLines(taxa, file.path(OUTPUT, "fasta", paste0(run_id, "_taxa.txt")))
  manifest[[length(manifest) + 1]] <<- tibble(
    run_id = run_id, true_lambda = lambda, keep_frac = keep_frac,
    rep_id = replicate, n_taxa = length(taxa), n_var_sites = ncol(m))
}

for (lambda in LAMBDAS) {
  keep_frac <- lambda / OBS_LAMBDA
  n_keep <- round(nrow(alignment) * keep_frac)
  for (replicate in seq_len(N_REPS)) {
    write_run(sprintf("lambda%.3f_rep%02d", lambda, replicate),
              lambda, keep_frac, replicate,
              sample(rownames(alignment), n_keep))
  }
}
write_run("observed", OBS_LAMBDA, 1, 0, rownames(alignment))

manifest <- bind_rows(manifest)
write_csv(manifest, file.path(OUTPUT, "manifest.csv"))

cat("runs written:", nrow(manifest), "\n\n")
print(as.data.frame(
  manifest %>%
    group_by(true_lambda) %>%
    summarise(mean_taxa = mean(n_taxa),
              mean_sites = mean(n_var_sites), .groups = "drop")))

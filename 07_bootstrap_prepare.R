#!/usr/bin/env Rscript
# Draw the bootstrap replicates.
#
# Each of N_BOOTSTRAP iterations resamples the cohort with replacement. A
# phylogeny cannot be built from duplicated sequences, so the alignment written
# for each replicate holds the unique isolates drawn, and the multiplicity of
# each isolate is stored as a resampling weight to be applied at the
# classification step. Drawing n isolates with replacement from n leaves about
# 63% of them distinct, so each replicate tree has roughly two thirds the taxa
# of the full cohort.
#
# Usage:  SIMEX_ROOT=/path/to/project Rscript 07_bootstrap_prepare.R
# Writes: output/bootstrap/fasta/, output/bootstrap/manifest.csv
suppressPackageStartupMessages({
  library(ape); library(readr); library(dplyr)
})
script_dir <- function() {
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "."
}
source(file.path(script_dir(), "common.R"))

SUB <- file.path(OUTPUT, "bootstrap")
for (d in c("fasta", "trees", "trees_rooted", "snppar_runs", "logs", "results")) {
  dir.create(file.path(SUB, d), showWarnings = FALSE, recursive = TRUE)
}
set.seed(SEED_BOOTSTRAP)

alignment <- read_alignment()
taxa <- rownames(alignment)
cat("alignment:", nrow(alignment), "taxa x", ncol(alignment), "sites\n")
cat("bootstrap iterations:", N_BOOTSTRAP, "\n\n")

manifest <- vector("list", N_BOOTSTRAP)
for (iteration in seq_len(N_BOOTSTRAP)) {
  run_id <- sprintf("boot%03d", iteration)
  draw <- sample(taxa, length(taxa), replace = TRUE)
  weights <- table(draw)
  unique_taxa <- names(weights)

  m <- drop_invariant(alignment[unique_taxa, , drop = FALSE])
  write_fasta(m, file.path(SUB, "fasta", paste0(run_id, ".fa")))
  writeLines(unique_taxa, file.path(SUB, "fasta", paste0(run_id, "_taxa.txt")))
  write_csv(tibble(taxon = unique_taxa, weight = as.integer(weights)),
            file.path(SUB, "fasta", paste0(run_id, "_weights.csv")))

  manifest[[iteration]] <- tibble(run_id = run_id,
                                  n_unique = length(unique_taxa),
                                  n_drawn = length(draw),
                                  n_var_sites = ncol(m))
}

manifest <- bind_rows(manifest)
write_csv(manifest, file.path(SUB, "manifest.csv"))

cat("replicates:", nrow(manifest), "\n")
cat("unique taxa per replicate: min", min(manifest$n_unique),
    "median", median(manifest$n_unique), "max", max(manifest$n_unique), "\n")
cat("variable sites: min", min(manifest$n_var_sites),
    "median", median(manifest$n_var_sites),
    "max", max(manifest$n_var_sites), "\n")

#!/usr/bin/env Rscript
# Weighted classification over the bootstrap replicates, with percentile
# confidence intervals.
#
# The classification rules are those of 06_classify_and_extrapolate.R, applied
# to the primary scenario, with each isolate contributing its resampling
# weight. Confidence intervals are the 2.5th and 97.5th percentiles of the
# distribution across replicates, so they carry uncertainty from both sampling
# variability and phylogenetic reconstruction, the tree being rebuilt in every
# iteration.
#
# Usage:  SIMEX_ROOT=/path/to/project Rscript 08_bootstrap_analyse.R
# Writes: output/results/bootstrap.csv
suppressPackageStartupMessages({
  library(ape); library(dplyr); library(readr); library(stringr)
  library(purrr); library(tidyr)
})
script_dir <- function() {
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "."
}
source(file.path(script_dir(), "common.R"))

SUB <- file.path(OUTPUT, "bootstrap")
meta <- load_metadata()
resistant <- meta$resistant
mappable  <- meta$mappable
previous  <- meta$previous

classify_replicate <- function(run) {
  events_path <- file.path(SUB, "snppar_runs", run, "all_mutation_events.tsv.gz")
  if (!file.exists(events_path)) return(NULL)

  events <- read_tsv(gzfile(events_path), show_col_types = FALSE, progress = FALSE)
  tree <- read.tree(file.path(SUB, "snppar_runs", run, "node_labelled_newick.tre"))
  tree$tip.label <- normalise_id(tree$tip.label)

  weights_table <- read_csv(
    file.path(SUB, "fasta", paste0(run, "_weights.csv")), show_col_types = FALSE)
  weights <- setNames(weights_table$weight, normalise_id(weights_table$taxon))

  study <- names(weights)
  resistant_here <- intersect(study, resistant)
  mappable_here  <- intersect(study, mappable)
  treated <- previous[study]
  names(treated) <- study

  events <- events %>%
    mutate(Codon = suppressWarnings(as.integer(Codon)),
           Derived_Node = normalise_id(Derived_Node)) %>%
    filter(Change == "NS", is_rrdr(Codon), Gene == RESISTANCE_GENE)

  tip_events <- unique(events$Derived_Node[
    !is_internal_node(events$Derived_Node) &
      events$Derived_Node %in% mappable_here])
  internal_nodes <- unique(events$Derived_Node[
    is_internal_node(events$Derived_Node)])

  descendants <- descendants_by_node(tree)
  n_internal_treated <- 0L
  for (node in internal_nodes) {
    below <- intersect(descendants[[node]], mappable_here)
    if (!length(below)) next
    if (any(treated[below] == 1L, na.rm = TRUE)) {
      n_internal_treated <- n_internal_treated + 1L
    }
  }

  trials <- sum(weights[resistant_here])
  acquired_N <- sum(weights[tip_events][treated[tip_events] == 1L], na.rm = TRUE)
  acquired_N1 <- min(acquired_N + n_internal_treated, trials)

  tibble(run_id = run, n_unique = length(study), trials = trials,
         acquired_N = acquired_N, acquired_N1 = acquired_N1,
         transmitted_N  = 100 * (1 - acquired_N / trials),
         transmitted_N1 = 100 * (1 - acquired_N1 / trials))
}

manifest <- read_csv(file.path(SUB, "manifest.csv"), show_col_types = FALSE)
results <- map_dfr(manifest$run_id, classify_replicate)
write_csv(results, file.path(OUTPUT, "results", "bootstrap.csv"))

percentiles <- function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE)

cat("\n=== BOOTSTRAP (", nrow(results), "replicates ) ===\n\n")
cat(sprintf("unique taxa per replicate : median %.0f\n", median(results$n_unique)))
cat(sprintf("weighted denominator      : median %.0f\n\n", median(results$trials)))

for (column in c("transmitted_N", "transmitted_N1")) {
  v <- percentiles(results[[column]])
  label <- if (column == "transmitted_N") "N" else "N-1"
  cat(sprintf("%-4s transmitted: median %.1f%%  95%% CI %.1f-%.1f%%\n",
              label, v[2], v[1], v[3]))
  cat(sprintf("%-4s acquired   : median %.1f%%  95%% CI %.1f-%.1f%%\n",
              "", 100 - v[2], 100 - v[3], 100 - v[1]))
}
cat(sprintf("\nplausible range of acquired resistance (N to N-1): %.1f-%.1f%%\n",
            100 - percentiles(results$transmitted_N)[2],
            100 - percentiles(results$transmitted_N1)[2]))

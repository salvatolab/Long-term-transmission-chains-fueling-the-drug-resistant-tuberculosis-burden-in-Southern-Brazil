#!/usr/bin/env Rscript
# Classify rifampicin resistance as transmitted or de novo acquired, and fit
# the SIMEX extrapolation.
#
# Classification rules
# --------------------
# The denominator is the set of rifampicin-resistant cases, defined from the
# isolate's own genotype (any resistance-associated rpoB mutation reported by
# TB-Profiler), independently of the tree. Susceptible isolates nested inside a
# resistant clade are therefore never counted.
#
# A resistance mutation placed on a terminal branch in a previously treated
# patient is a de novo acquisition.
#
# Where a resistance mutation is placed on an internal node, the descendant in
# which resistance first emerged cannot be identified, so two bounding rules
# are applied:
#   N    all descendants are counted as transmitted. Lower bound for de novo
#        acquisition.
#   N-1  exactly one de novo acquisition is assigned per such node when at
#        least one descendant patient had previous TB treatment, irrespective
#        of how many descendants were previously treated; the remaining
#        descendants are counted as transmitted. Upper bound.
# The total number of acquisitions is capped at the number of resistant cases.
#
# Resistant isolates with no mappable event (indel-only, or minority-frequency
# alleles) remain in the denominator and are counted as transmitted by default.
# The count is reported so that this assumption is visible.
#
# Scenarios
# ---------
# The primary scenario uses the genotype-defined resistant set as denominator
# and restricts mutations to the resistance gene. Three sensitivity scenarios
# vary the denominator and the gene filter.
#
# Usage:  SIMEX_ROOT=/path/to/project Rscript 06_classify_and_extrapolate.R
# Writes: output/results/simulations.csv, output/results/extrapolation.csv
suppressPackageStartupMessages({
  library(ape); library(dplyr); library(readr); library(stringr)
  library(purrr); library(mgcv); library(tidyr)
})
script_dir <- function() {
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "."
}
source(file.path(script_dir(), "common.R"))

meta <- load_metadata()
resistant <- meta$resistant
mappable  <- meta$mappable
previous  <- meta$previous

observed_taxa <- normalise_id(
  readLines(file.path(OUTPUT, "fasta", "observed_taxa.txt")))
missing <- setdiff(observed_taxa, names(previous))
if (length(missing)) {
  stop(sprintf("%d taxa have no metadata row: %s",
               length(missing), paste(head(missing, 5), collapse = ", ")))
}

cat("rifampicin-resistant (any rpoB mutation):", length(resistant), "\n")
cat("  with a mappable RRDR SNV              :", length(mappable), "\n")
cat("  without a mappable event              :",
    length(setdiff(resistant, mappable)), "\n\n")

# ------------------------------------------------------------- classification
classify_run <- function(run, gene_filter, denominator) {
  events_path <- file.path(OUTPUT, "snppar_runs", run, "all_mutation_events.tsv.gz")
  if (!file.exists(events_path)) return(NULL)

  events <- read_tsv(gzfile(events_path), show_col_types = FALSE, progress = FALSE)
  tree <- read.tree(file.path(OUTPUT, "snppar_runs", run, "node_labelled_newick.tre"))
  tree$tip.label <- normalise_id(tree$tip.label)

  study <- normalise_id(
    readLines(file.path(OUTPUT, "fasta", paste0(run, "_taxa.txt"))))
  resistant_here <- intersect(study, resistant)
  mappable_here  <- intersect(study, mappable)
  treated <- previous[study]
  names(treated) <- study

  events <- events %>%
    mutate(Codon = suppressWarnings(as.integer(Codon)),
           Derived_Node = normalise_id(Derived_Node)) %>%
    filter(Change == "NS", is_rrdr(Codon))
  if (gene_filter) events <- filter(events, Gene == RESISTANCE_GENE)

  tip_events <- unique(events$Derived_Node[
    !is_internal_node(events$Derived_Node) &
      events$Derived_Node %in% mappable_here])
  internal_nodes <- unique(events$Derived_Node[
    is_internal_node(events$Derived_Node)])

  descendants <- descendants_by_node(tree)
  covered <- tip_events
  n_internal_treated <- 0L
  for (node in internal_nodes) {
    below <- intersect(descendants[[node]], mappable_here)
    if (!length(below)) next
    covered <- c(covered, below)
    if (any(treated[below] == 1L, na.rm = TRUE)) {
      n_internal_treated <- n_internal_treated + 1L
    }
  }
  covered <- unique(covered)

  trials <- switch(denominator,
                   resistant = length(resistant_here),
                   mappable  = length(mappable_here),
                   length(study))
  acquired_N <- sum(treated[tip_events] == 1L, na.rm = TRUE)

  tibble(run_id = run,
         n_study = length(study),
         n_resistant = length(resistant_here),
         n_mappable = length(mappable_here),
         n_covered = length(covered),
         n_tip_events = length(tip_events),
         n_internal_events = length(internal_nodes),
         trials = trials,
         acquired_N = acquired_N,
         acquired_N1 = min(acquired_N + n_internal_treated, trials))
}

# -------------------------------------------------------------- GAM and SIMEX
manifest <- read_csv(file.path(OUTPUT, "manifest.csv"), show_col_types = FALSE)
lambda_of <- setNames(manifest$true_lambda, manifest$run_id)

run_scenario <- function(gene_filter, denominator, label) {
  d <- map_dfr(manifest$run_id, classify_run,
               gene_filter = gene_filter, denominator = denominator)
  if (!nrow(d)) return(NULL)
  d <- mutate(d, lambda = lambda_of[run_id], scenario = label)
  obs <- filter(d, run_id == "observed")

  fits <- map_dfr(c("acquired_N", "acquired_N1"), function(column) {
    dd <- mutate(d, successes = trials - .data[[column]])
    # Binomial GAM with a thin-plate spline on lambda, fitted by REML, then
    # predicted at lambda = 1.0 (complete sampling).
    gam_fit <- gam(cbind(successes, trials - successes) ~
                     s(lambda, bs = "tp", k = 10),
                   family = binomial, data = dd, method = "REML")
    prediction <- predict(gam_fit, tibble(lambda = 1), type = "link",
                          se.fit = TRUE)
    fit <- as.numeric(prediction$fit)
    se  <- as.numeric(prediction$se.fit)

    # Diagnostics: significance of the linear trend, and the dispersion of a
    # quasibinomial fit. A non-significant trend, or a dispersion far from 1,
    # means the extrapolation rests on little information.
    linear_fit <- glm(cbind(successes, trials - successes) ~ lambda,
                      family = binomial, data = dd)
    quasi_fit  <- glm(cbind(successes, trials - successes) ~ lambda,
                      family = quasibinomial, data = dd)

    tibble(scenario = label,
           method = if (column == "acquired_N") "N" else "N-1",
           n_denominator = obs$trials[1],
           observed = 100 * (obs$trials[1] - obs[[column]][1]) / obs$trials[1],
           extrapolated = 100 * plogis(fit),
           ci_low  = 100 * plogis(fit - 1.96 * se),
           ci_high = 100 * plogis(fit + 1.96 * se),
           slope = coef(linear_fit)["lambda"],
           slope_p = summary(linear_fit)$coefficients["lambda", 4],
           edf = sum(summary(gam_fit)$edf),
           deviance_explained = 100 * summary(gam_fit)$dev.expl,
           dispersion = summary(quasi_fit)$dispersion)
  })
  list(data = d, fits = fits)
}

scenarios <- list(
  run_scenario(TRUE,  "resistant", "primary (denominator = resistant cases)"),
  run_scenario(TRUE,  "mappable",  "denominator = mappable subset"),
  run_scenario(TRUE,  "study",     "denominator = whole cohort"),
  run_scenario(FALSE, "study",     "no gene filter, denominator = whole cohort")
)
scenarios <- scenarios[!map_lgl(scenarios, is.null)]

fits <- map_dfr(scenarios, "fits")
data <- map_dfr(scenarios, "data")
write_csv(data, file.path(OUTPUT, "results", "simulations.csv"))
write_csv(fits, file.path(OUTPUT, "results", "extrapolation.csv"))

cat("=== % TRANSMITTED ===\n\n")
print(as.data.frame(fits %>% mutate(across(where(is.numeric), ~round(., 3)))))

primary <- filter(data, scenario == "primary (denominator = resistant cases)")
cat("\n=== means by lambda (primary scenario) ===\n\n")
print(as.data.frame(
  primary %>%
    group_by(lambda) %>%
    summarise(taxa = mean(n_study), resistant = mean(n_resistant),
              mappable = mean(n_mappable), covered = mean(n_covered),
              tip = mean(n_tip_events), internal = mean(n_internal_events),
              N  = round(100 * mean(1 - acquired_N / trials), 2),
              N1 = round(100 * mean(1 - acquired_N1 / trials), 2),
              .groups = "drop")))

obs <- filter(primary, run_id == "observed")
cat("\n=== agreement between tree and genotype (observed tree) ===\n")
cat(sprintf("rifampicin-resistant cases   : %d\n", obs$n_resistant))
cat(sprintf("  with an RRDR SNV           : %d\n", obs$n_mappable))
cat(sprintf("  with an event on the tree  : %d (%.1f%%)\n",
            obs$n_covered, 100 * obs$n_covered / obs$n_resistant))
cat(sprintf("  without a mappable event   : %d (counted as transmitted)\n",
            obs$n_resistant - obs$n_covered))

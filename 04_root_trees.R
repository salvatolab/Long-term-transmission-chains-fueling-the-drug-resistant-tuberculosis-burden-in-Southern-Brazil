#!/usr/bin/env Rscript
# Root every tree on a single common reference group.
#
# The reference group is defined once, from the smaller side of the root
# bipartition of the complete-cohort tree, and then applied to every replicate.
# Rooting each replicate independently by midpoint would place topologically
# different roots in different subsamples and would change the inferred
# direction of mutation, which is exactly the quantity being classified.
#
# A replicate that retains no member of the reference group, or that consists
# entirely of its members, falls back to midpoint rooting; the count of such
# replicates is reported and should be zero or near zero.
#
# Usage:  SIMEX_ROOT=/path/to/project Rscript 04_root_trees.R [SUBDIR]
suppressPackageStartupMessages({
  library(ape); library(readr); library(phytools)
})
script_dir <- function() {
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "."
}
source(file.path(script_dir(), "common.R"))

args <- commandArgs(trailingOnly = TRUE)
SUB <- if (length(args) >= 1) args[1] else OUTPUT
dir.create(file.path(SUB, "trees_rooted"), showWarnings = FALSE, recursive = TRUE)

reference_group_path <- file.path(OUTPUT, "root_group.txt")

# Define the reference group from the complete-cohort tree, once.
if (!file.exists(reference_group_path)) {
  observed <- read.tree(file.path(SUB, "trees", "observed.treefile"))
  observed <- midpoint.root(observed)
  write.tree(observed, file.path(SUB, "trees_rooted", "observed.treefile"))

  root_children <- observed$edge[observed$edge[, 1] == (Ntip(observed) + 1), 2]
  sides <- lapply(root_children, function(child) {
    if (child <= Ntip(observed)) observed$tip.label[child]
    else extract.clade(observed, child)$tip.label
  })
  reference_group <- sides[[which.min(lengths(sides))]]
  writeLines(reference_group, reference_group_path)
  cat("reference group defined:", length(reference_group), "taxa\n")
} else {
  reference_group <- readLines(reference_group_path)
  cat("reference group loaded:", length(reference_group), "taxa\n")
}

manifest <- read_csv(file.path(SUB, "manifest.csv"), show_col_types = FALSE)
manifest <- manifest[manifest$run_id != "observed" |
                       !file.exists(file.path(SUB, "trees_rooted",
                                              "observed.treefile")), ]

n_rooted <- 0L
n_fallback <- 0L
group_sizes <- integer(0)

for (run in manifest$run_id) {
  tree_path <- file.path(SUB, "trees", paste0(run, ".treefile"))
  if (!file.exists(tree_path)) next
  tree <- read.tree(tree_path)

  outgroup <- intersect(reference_group, tree$tip.label)
  group_sizes <- c(group_sizes, length(outgroup))

  rooted <- NULL
  if (length(outgroup) > 0 && length(outgroup) < Ntip(tree)) {
    rooted <- tryCatch(root(tree, outgroup = outgroup, resolve.root = TRUE),
                       error = function(e) NULL)
  }
  if (is.null(rooted)) {
    rooted <- midpoint.root(tree)
    n_fallback <- n_fallback + 1L
  }
  write.tree(rooted, file.path(SUB, "trees_rooted", paste0(run, ".treefile")))
  n_rooted <- n_rooted + 1L
}

cat("rooted:", n_rooted, "| midpoint fallback:", n_fallback, "\n")
if (length(group_sizes)) {
  cat("reference-group taxa per replicate: min", min(group_sizes),
      "median", median(group_sizes), "max", max(group_sizes), "\n")
}

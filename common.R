#!/usr/bin/env Rscript
# Shared configuration and helper functions for the TB-SIMEX reanalysis.
# Sourced by every R script in this directory.

# ---------------------------------------------------------------- directories
# SIMEX_ROOT is the working directory holding input/ and output/.
ROOT   <- Sys.getenv("SIMEX_ROOT", unset = ".")
INPUT  <- file.path(ROOT, "input")
OUTPUT <- file.path(ROOT, "output")

# -------------------------------------------------------------- SIMEX design
# Observed sampling fraction: the proportion of notified rifampicin-resistant
# cases in the study period that were sequenced. Set this for your own cohort.
OBS_LAMBDA <- as.numeric(Sys.getenv("OBS_LAMBDA", unset = "0.25"))

# Sampling-fraction grid. The TB-SIMEX reference pipeline spans lambda = 0.05
# up to 0.05 below the observed sampling fraction, in N_LEVELS equal steps.
# At the 0.64 coverage of Spies et al. this rule retains 8-92% of sequenced
# isolates; at a coverage of 0.25 it retains 20-80%.
N_LEVELS <- as.integer(Sys.getenv("N_LEVELS", unset = "10"))
N_REPS   <- as.integer(Sys.getenv("N_REPS",   unset = "5"))
LAMBDAS  <- seq(0.05, OBS_LAMBDA - 0.05, length.out = N_LEVELS)

# Random seeds. Subsampling and bootstrap draw from separate streams so that
# either can be rerun without disturbing the other.
SEED_SUBSAMPLE <- as.integer(Sys.getenv("SEED_SUBSAMPLE", unset = "20261110"))
SEED_BOOTSTRAP <- as.integer(Sys.getenv("SEED_BOOTSTRAP", unset = "20261111"))
N_BOOTSTRAP    <- as.integer(Sys.getenv("N_BOOTSTRAP",    unset = "100"))

# ------------------------------------------------------- resistance definition
# Locus tag of the resistance gene in the reference annotation (rpoB in H37Rv).
RESISTANCE_GENE <- Sys.getenv("RESISTANCE_GENE", unset = "Rv0667")

# Rifampicin resistance-determining region: codon 170, codons 426-452 and
# codon 491, in Mycobacterium tuberculosis numbering.
is_rrdr <- function(codon) {
  codon == 170L | codon == 491L | (codon >= 426L & codon <= 452L)
}

# --------------------------------------------------------- generic utilities

# Strip a trailing collection-year suffix (_23, _24, ...) from isolate labels.
# Tree tip labels and metadata identifiers must agree exactly; silently
# mismatched names would be coerced to "untreated" during classification.
normalise_id <- function(x) sub("_(2[0-9])$", "", x)

# TRUE for SNPPar internal-node labels (N1, N2, ...), FALSE for tip labels.
is_internal_node <- function(x) grepl("^N[0-9]+$", x)

# Drop alignment columns that are invariant within the given set of taxa.
# Not doing this inside each subsample lets the number of informative sites
# grow with lambda, which confounds the SIMEX trend, and leaves constant sites
# in an alignment analysed under an ascertainment-bias correction.
drop_invariant <- function(m) {
  keep <- apply(m, 2, function(column) {
    length(unique(column[column %in% c("A", "C", "G", "T")])) > 1
  })
  m[, keep, drop = FALSE]
}

write_fasta <- function(m, path) {
  connection <- file(path, "w")
  on.exit(close(connection))
  for (i in seq_len(nrow(m))) {
    cat(">", rownames(m)[i], "\n", sep = "", file = connection)
    cat(paste(m[i, ], collapse = ""), "\n", sep = "", file = connection)
  }
}

read_alignment <- function(path = file.path(INPUT, "alignment.fasta")) {
  aln <- ape::read.dna(path, format = "fasta", as.character = TRUE)
  aln <- toupper(aln)
  rownames(aln) <- sub("\\s.*$", "", rownames(aln))
  aln
}

# Map every internal node label to the tip labels descending from it.
descendants_by_node <- function(tree) {
  n_tips <- ape::Ntip(tree)
  n_nodes <- tree$Nnode
  descendants <- vector("list", n_tips + n_nodes)
  for (i in seq_len(n_tips)) descendants[[i]] <- tree$tip.label[i]
  for (node in (n_tips + n_nodes):(n_tips + 1)) {
    children <- tree$edge[tree$edge[, 1] == node, 2]
    descendants[[node]] <- unlist(descendants[children])
  }
  setNames(descendants[(n_tips + 1):(n_tips + n_nodes)], tree$node.label)
}

# ------------------------------------------------------------------ metadata
# input/metadata.tsv must contain one row per isolate and these columns:
#   ID                   isolate identifier, matching the alignment labels
#   resistance evidence  TB-Profiler resistance calls, "-" when susceptible
#   prev_tb              1 if the patient had previous TB treatment, else 0
#
# Returns a list with:
#   resistant  identifiers carrying any resistance-associated rpoB mutation
#   mappable   the subset carrying a non-synonymous RRDR SNV, the only form of
#              mutation parsimony can place on the tree; indel-only isolates
#              and minority-frequency alleles produce no mappable event
#   previous   named integer vector of previous-treatment status
load_metadata <- function(path = file.path(INPUT, "metadata.tsv")) {
  meta <- readr::read_tsv(path, show_col_types = FALSE)
  meta <- dplyr::rename(meta, id = ID, evidence = `resistance evidence`)
  meta$evidence <- tidyr::replace_na(meta$evidence, "-")
  stopifnot(!any(duplicated(meta$id)), !any(is.na(meta$prev_tb)))

  rrdr_codons <- function(evidence) {
    matches <- stringr::str_match_all(
      evidence, "rpoB p\\.[A-Za-z]{3}(\\d+)[A-Za-z]{3}")[[1]]
    if (!nrow(matches)) return(integer(0))
    as.integer(matches[, 2])
  }
  has_mappable <- function(evidence) {
    codons <- rrdr_codons(evidence)
    length(codons) > 0 && any(is_rrdr(codons))
  }

  previous <- setNames(as.integer(meta$prev_tb), normalise_id(meta$id))
  stopifnot(!any(duplicated(names(previous))))

  list(
    resistant = normalise_id(meta$id[meta$evidence != "-"]),
    mappable  = normalise_id(meta$id[vapply(meta$evidence, has_mappable, logical(1))]),
    previous  = previous
  )
}

dir.create(OUTPUT, showWarnings = FALSE, recursive = TRUE)

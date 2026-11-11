# Long-term-transmission-chains-fueling-the-drug-resistant-tuberculosis-burden-in-Southern-Brazil
Supplementar methods and results to the manusctipt Long-term transmission chains fueling the drug-resistant tuberculosis burden in Southern Brazil


# TB-SIMEX analysis — transmitted versus acquired rifampicin resistance

Scripts and parameter settings used to estimate the relative contribution of
transmitted and *de novo* acquired rifampicin resistance in a population-based
collection of *Mycobacterium tuberculosis* genomes, following the TB-SIMEX
framework of Spies et al. https://doi.org/10.1016/j.lanwpc.2026.101906

## Software versions

| Tool | Version | Role |
|---|---|---|
| IQ-TREE | 3.0.1 | maximum-likelihood phylogenies |
| SNPPar | 1.0 | parsimony mapping of mutations onto trees |
| R | ≥ 4.3 | subsampling, classification, GAM extrapolation, figures |
| Python | ≥ 3.9 | SNP table construction and per-run table subsetting |

R packages: `ape`, `phytools`, `dplyr`, `readr`, `tidyr`, `stringr`, `purrr`,
`mgcv`, `ggplot2`, `scales`. On macOS, `showtext` and `sysfonts` are used for
font rendering when available.

## Analysis parameters

| Parameter | Value | Set in |
|---|---|---|
| Observed sampling fraction (λ) | 0.25 | `OBS_LAMBDA` |
| Sampling-fraction grid | λ = 0.05 to 0.20, ten evenly spaced levels | `LAMBDAS` |
| Equivalent retained fraction of sequenced isolates | 20–80% | derived |
| Replicates per sampling fraction | 5 (50 simulated datasets) | `N_REPS` |
| Bootstrap iterations | 100 | `N_BOOTSTRAP` |
| Random seed, subsampling | 20261110 | `SEED_SUBSAMPLE` |
| Random seed, bootstrap | 20261111 | `SEED_BOOTSTRAP` |
| Substitution model | TVM+F+ASC | `MODEL` |
| Resistance gene | Rv0667 (*rpoB*) | `RESISTANCE_GENE` |
| RRDR codons | 170, 426–452, 491 | `is_rrdr()` |
| Minimum read depth / genotype quality | DP ≥ 10, GQ ≥ 20 | `01_build_snp_table.py` |
| GAM | binomial, thin-plate spline on λ, k = 10, REML | `06_…R`, `09_figure.R` |

## Classification rules

The denominator is the set of rifampicin-resistant cases, defined from each
isolate's own genotype — any resistance-associated *rpoB* mutation reported by
TB-Profiler — independently of the phylogeny. 

A resistance mutation placed on a terminal branch in a previously treated
patient is a *de novo* acquisition.

Where a resistance mutation is placed on an internal node, the descendant in
which resistance first emerged cannot be identified, so two bounding rules are
applied:

- **N** — all descendants are counted as transmitted. Lower bound for *de novo*
  acquisition.
- **N−1** — exactly one *de novo* acquisition is assigned per such node when at
  least one descendant patient had a history of previous TB treatment,
  **irrespective of how many descendants were previously treated**; the
  remaining descendants are counted as transmitted. Upper bound.

The total number of acquisitions is capped at the number of resistant cases.

Resistant isolates carrying only an indel, or a minority-frequency allele,
produce no event that parsimony can place on the tree. They remain in the
denominator and are counted as transmitted by default;
`06_classify_and_extrapolate.R` reports how many isolates this affects.

Sites invariant within a subsample are removed from that subsample before tree
reconstruction, because the ascertainment-bias correction (+ASC) assumes an
alignment of variable sites, and because retaining them lets the number of
informative sites grow with λ and confounds the SIMEX trend.

Every tree, including the one for the complete cohort, is built by the same
procedure and rooted on a single common reference group, derived once from the
root bipartition of the complete-cohort tree. Rooting each replicate
independently by midpoint would alter the inferred direction of mutation.

## Uncertainty

Two intervals are reported.

The **bootstrap** interval applies to the estimate at the observed sampling
fraction. Each of 100 iterations resamples the cohort with replacement,
rebuilds the phylogeny from the unique isolates drawn, remaps mutations and
reclassifies using the resampling weights. Confidence intervals are the 2.5th
and 97.5th percentiles across iterations, so they carry uncertainty from both
sampling variability and phylogenetic inference.

The **GAM** interval applies to the extrapolated estimate at λ = 1.0 and
reflects only the uncertainty of the fitted curve.
`06_classify_and_extrapolate.R` also reports the significance of the linear
trend, the effective degrees of freedom of the spline, and the quasibinomial
dispersion, so that the information supporting the extrapolation can be judged.

## Input files

Place these in `input/`:

| File | Description |
|---|---|
| `alignment.fasta` | variable-site alignment, one sequence per isolate |
| `snp_table.csv` | SNPPar table: `Pos` column, then one column per isolate |
| `metadata.tsv` | isolate metadata, see below |
| `NC_000962.3.gbk` | H37Rv reference GenBank record, for SNPPar annotation |

`alignment.fasta` and `snp_table.csv` are produced by
`scripts/01_build_snp_table.py` from annotated per-isolate VCFs, a reference
genome and a BED file of PE/PPE regions to mask.

`metadata.tsv` is tab-separated with one row per isolate and these columns:

| Column | Description |
|---|---|
| `ID` | isolate identifier, matching the alignment labels |
| `resistance evidence` | TB-Profiler resistance calls; `-` when susceptible |
| `prev_tb` | `1` if the patient had previous TB treatment, `0` otherwise |


## Running

```bash
# once, from the variant calls
python3 scripts/01_build_snp_table.py VCF_DIR reference.fasta pe_ppe.bed input/

# the full pipeline
bash run_all.sh
```

Or step by step:

```bash
export SIMEX_ROOT=.
Rscript scripts/02_subsample.R
bash    scripts/03_build_trees.sh    output
Rscript scripts/04_root_trees.R      output
bash    scripts/05_map_mutations.sh  output
Rscript scripts/06_classify_and_extrapolate.R
Rscript scripts/07_bootstrap_prepare.R
bash    scripts/03_build_trees.sh    output/bootstrap
Rscript scripts/04_root_trees.R      output/bootstrap
bash    scripts/05_map_mutations.sh  output/bootstrap
Rscript scripts/08_bootstrap_analyse.R
Rscript scripts/09_figure.R
```

Set `SNPPAR` if the SNPPar executable is not on `PATH`, and `PARALLEL_JOBS`,
`SNPPAR_JOBS` and `IQTREE_THREADS` to match the available hardware. Building
the 51 subsample trees and the 100 bootstrap trees is the limiting step.

## Outputs

| File | Contents |
|---|---|
| `output/manifest.csv` | one row per simulated dataset, with λ and site counts |
| `output/results/simulations.csv` | per-run classification counts, all scenarios |
| `output/results/extrapolation.csv` | observed and extrapolated estimates, with CIs and GAM diagnostics |
| `output/results/bootstrap.csv` | per-iteration bootstrap estimates |
| `output/results/figure_simex.pdf` | main figure |

Besides the primary scenario, `extrapolation.csv` reports three sensitivity
scenarios that vary the denominator (resistant cases, the mappable subset, or
the whole cohort) and the gene filter.

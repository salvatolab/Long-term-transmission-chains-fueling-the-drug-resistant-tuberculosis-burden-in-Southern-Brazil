#!/usr/bin/env bash
# Reconstruct one maximum-likelihood phylogeny per run listed in the manifest.
#
# Every tree, including the one for the complete cohort, is built by the same
# procedure and under the same substitution model, so that observed and
# simulated estimates are not separated by differences in model or in site
# filtering. The ascertainment-bias correction (+ASC) requires an alignment
# with no constant sites, which 02_subsample.R guarantees.
#
# Usage:  SIMEX_ROOT=/path/to/project bash 03_build_trees.sh [SUBDIR]
#   SUBDIR  optional output subdirectory (default: output).
#           Pass "output/bootstrap" to build the bootstrap replicates.
#
# Environment:
#   MODEL           substitution model (default: TVM+F+ASC)
#   IQTREE_THREADS  threads per IQ-TREE process (default: 2)
#   PARALLEL_JOBS   concurrent IQ-TREE processes (default: 4)
set -uo pipefail

ROOT="${SIMEX_ROOT:-.}"
SUB="${1:-$ROOT/output}"
MODEL="${MODEL:-TVM+F+ASC}"
IQTREE_THREADS="${IQTREE_THREADS:-2}"
PARALLEL_JOBS="${PARALLEL_JOBS:-4}"

[[ -f "$SUB/manifest.csv" ]] || { echo "no manifest.csv in $SUB" >&2; exit 1; }
mkdir -p "$SUB/trees" "$SUB/logs"

echo "[$(date +%T)] building trees in $SUB (model $MODEL)"

build_tree() {
    local run="$1" sub="$2" model="$3" threads="$4"
    [[ -f "$sub/trees/${run}.treefile" ]] && return 0
    iqtree -s "$sub/fasta/${run}.fa" -m "$model" -nt "$threads" \
        -pre "$sub/trees/$run" -redo \
        > "$sub/logs/${run}_iqtree.log" 2>&1 || true
}
export -f build_tree

tail -n +2 "$SUB/manifest.csv" | cut -d, -f1 \
  | xargs -P "$PARALLEL_JOBS" -I{} \
      bash -c 'build_tree "$0" "$1" "$2" "$3"' {} "$SUB" "$MODEL" "$IQTREE_THREADS"

echo "[$(date +%T)] trees built: $(ls "$SUB"/trees/*.treefile 2>/dev/null | wc -l)"

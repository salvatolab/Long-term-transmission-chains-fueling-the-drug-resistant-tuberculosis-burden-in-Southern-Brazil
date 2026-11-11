#!/usr/bin/env bash
# Map mutations onto each rooted tree by maximum parsimony, using SNPPar.
#
# For each run the master SNP table is reduced to the columns of that run's
# taxa, and rows that no longer show at least two distinct bases are dropped.
# SNPPar counts only A/C/G/T as character states, so a row left with a single
# observed base would be discarded by SNPPar itself; filtering it here keeps
# the input table consistent with the alignment the tree was built from.
#
# Usage:  SIMEX_ROOT=/path/to/project bash 05_map_mutations.sh [SUBDIR]
#
# Environment:
#   SNPPAR        path to the snppar executable (default: snppar on PATH)
#   GENBANK       reference GenBank record (default: input/NC_000962.3.gbk)
#   SNPPAR_JOBS   concurrent SNPPar processes (default: 6)
set -uo pipefail

ROOT="${SIMEX_ROOT:-.}"
SUB="${1:-$ROOT/output}"
SNPPAR="${SNPPAR:-snppar}"
GENBANK="${GENBANK:-$ROOT/input/NC_000962.3.gbk}"
SNPPAR_JOBS="${SNPPAR_JOBS:-6}"
SNP_TABLE="${SNP_TABLE:-$ROOT/input/snp_table.csv}"

[[ -f "$SUB/manifest.csv" ]] || { echo "no manifest.csv in $SUB" >&2; exit 1; }
[[ -f "$GENBANK" ]] || { echo "GenBank record not found: $GENBANK" >&2; exit 1; }
[[ -f "$SNP_TABLE" ]] || { echo "SNP table not found: $SNP_TABLE" >&2; exit 1; }
mkdir -p "$SUB/snppar_runs" "$SUB/logs"

echo "[$(date +%T)] mapping mutations in $SUB"

map_one() {
    local run="$1" sub="$2" snppar="$3" genbank="$4" table="$5"
    local dir="$sub/snppar_runs/$run"
    local tree="$sub/trees_rooted/${run}.treefile"
    [[ -f "$dir/all_mutation_events.tsv.gz" ]] && return 0
    [[ -f "$tree" ]] || return 0
    mkdir -p "$dir"

    local subset="$dir/subset.csv"
    python3 -I - "$sub/fasta/${run}.fa" "$table" "$subset" <<'PYEOF'
import sys, csv
fasta_path, table_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
tips = {line[1:].strip() for line in open(fasta_path) if line.startswith(">")}
with open(table_path, newline="") as src, open(out_path, "w", newline="\n") as dst:
    reader, writer = csv.reader(src), csv.writer(dst, lineterminator="\n")
    header = next(reader)
    keep = [0] + [i for i, name in enumerate(header)
                  if i > 0 and name.strip() in tips]
    writer.writerow([header[i] for i in keep])
    for row in reader:
        values = [row[i] for i in keep]
        if len({v for v in values[1:]} & {"A", "C", "G", "T"}) > 1:
            writer.writerow(values)
PYEOF

    "$snppar" -s "$subset" -t "$tree" -g "$genbank" -d "$dir" \
        > "$sub/logs/${run}_snppar.log" 2>&1 || true
    rm -f "$subset"
    [[ -f "$dir/all_mutation_events.tsv" ]] && gzip -f "$dir/all_mutation_events.tsv"
    return 0
}
export -f map_one

tail -n +2 "$SUB/manifest.csv" | cut -d, -f1 \
  | xargs -P "$SNPPAR_JOBS" -I{} \
      bash -c 'map_one "$0" "$1" "$2" "$3" "$4"' \
      {} "$SUB" "$SNPPAR" "$GENBANK" "$SNP_TABLE"

echo "[$(date +%T)] runs completed: $(ls -d "$SUB"/snppar_runs/*/all_mutation_events.tsv.gz 2>/dev/null | wc -l)"

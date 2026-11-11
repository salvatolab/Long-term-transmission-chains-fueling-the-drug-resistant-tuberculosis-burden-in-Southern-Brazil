#!/usr/bin/env python3
"""Build the SNP table and whole-genome SNP alignment from annotated VCFs.

The SNP table is written in the format SNPPar expects: one row per variable
position, one column per isolate, and an explicit base call in every cell.
Cells are filled with the reference base where an isolate has no variant call,
because SNPPar treats only A/C/G/T as character states and silently discards
any site that does not show at least two distinct bases.

Filters applied, in order:
  1. indels and multi-allelic records are dropped (single-base REF and ALT only)
  2. positions inside the PE/PPE mask are dropped
  3. records whose REF does not match the reference genome are dropped
  4. calls below the depth or genotype-quality thresholds are dropped

Usage:
    python3 01_build_snp_table.py VCF_DIR REFERENCE.fasta PE_PPE.bed OUTPUT_DIR

Arguments:
    VCF_DIR         directory of annotated per-isolate VCFs
    REFERENCE.fasta H37Rv reference genome (NC_000962.3)
    PE_PPE.bed      BED file of PE/PPE regions to mask
    OUTPUT_DIR      destination for snp_table.csv and alignment.fasta

Outputs:
    OUTPUT_DIR/snp_table.csv    SNPPar input table
    OUTPUT_DIR/alignment.fasta  variable-site alignment, same sites and order
"""
import sys
import os
import glob
import csv

# Minimum read depth and genotype quality for a variant call to be retained.
MIN_DP, MIN_GQ = 10, 20

# Glob pattern matching the VCF files inside VCF_DIR. The isolate identifier is
# taken from the file name: everything between VCF_PREFIX and the first dot.
VCF_PREFIX = "corrigido_"
VCF_GLOB = VCF_PREFIX + "*.concordant.eff.vcf"


def load_reference(path):
    """Return the reference sequence as an uppercase string."""
    chunks = []
    with open(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                chunks.append(line.strip())
    return "".join(chunks).upper()


def load_mask(path, genome_length):
    """Return a 1-based byte mask flagging positions to exclude."""
    mask = bytearray(genome_length + 1)
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 3:
                continue
            start, end = int(fields[1]), int(fields[2])
            for position in range(start + 1, min(end, genome_length) + 1):
                mask[position] = 1
    return mask


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    vcf_dir, reference_path, bed_path, output_dir = sys.argv[1:5]

    reference = load_reference(reference_path)
    print("reference: %d bp" % len(reference), flush=True)

    mask = load_mask(bed_path, len(reference))
    print("masked positions (PE/PPE): %d" % sum(mask), flush=True)

    files = sorted(glob.glob(os.path.join(vcf_dir, VCF_GLOB)))
    if not files:
        sys.exit("no VCF files matched %s in %s" % (VCF_GLOB, vcf_dir))

    samples = []
    calls = {}
    n_indel = n_masked = n_lowqual = n_refmismatch = 0

    for path in files:
        sample = os.path.basename(path)[len(VCF_PREFIX):].split(".")[0]
        samples.append(sample)
        with open(path) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10:
                    continue
                position = int(fields[1])
                ref_base, alt_base = fields[3].upper(), fields[4].upper()

                if len(ref_base) != 1 or len(alt_base) != 1 or "," in alt_base:
                    n_indel += 1
                    continue
                if mask[position]:
                    n_masked += 1
                    continue
                if reference[position - 1] != ref_base:
                    n_refmismatch += 1
                    continue

                sample_fields = dict(zip(fields[8].split(":"), fields[9].split(":")))
                try:
                    depth = int(sample_fields.get("DP", 0))
                    quality = int(sample_fields.get("GQ", 0))
                    if depth < MIN_DP or quality < MIN_GQ:
                        n_lowqual += 1
                        continue
                except ValueError:
                    pass  # keep the call when DP/GQ are absent or non-numeric

                calls.setdefault(position, {})[sample] = alt_base

    print("isolates: %d" % len(samples))
    print("discarded -> indel/multi-allelic: %d, PE/PPE: %d, "
          "DP<%d or GQ<%d: %d, REF mismatch: %d"
          % (n_indel, n_masked, MIN_DP, MIN_GQ, n_lowqual, n_refmismatch))

    # Keep only positions showing at least two distinct bases across isolates,
    # counting the reference base for isolates without a variant call.
    kept = [
        position for position in sorted(calls)
        if len({calls[position].get(s, reference[position - 1]) for s in samples}
               & set("ACGT")) > 1
    ]
    print("variant positions: %d -> usable by SNPPar: %d" % (len(calls), len(kept)))

    os.makedirs(output_dir, exist_ok=True)

    table_path = os.path.join(output_dir, "snp_table.csv")
    with open(table_path, "w", newline="\n") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["Pos"] + samples)
        for position in kept:
            row = calls[position]
            writer.writerow(
                [position] + [row.get(s, reference[position - 1]) for s in samples]
            )

    alignment_path = os.path.join(output_dir, "alignment.fasta")
    with open(alignment_path, "w") as handle:
        for sample in samples:
            handle.write(">%s\n" % sample)
            handle.write(
                "".join(calls[p].get(sample, reference[p - 1]) for p in kept) + "\n"
            )

    print("\nwrote: %s" % table_path)
    print("wrote: %s  (%d taxa x %d sites)"
          % (alignment_path, len(samples), len(kept)))


if __name__ == "__main__":
    main()

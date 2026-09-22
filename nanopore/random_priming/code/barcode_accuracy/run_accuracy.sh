#!/bin/bash
# Reproduce the random-priming barcode-accuracy TABLE and FIGURE from the deposited raw FASTQs.
#
# This drives the main species-mixing barcode-accuracy program (speciesmix_bobcode_metrics.py),
# which aligns each barcode's reads with STARlong to the combined human+mouse genome and compares
# each read's bobcode-declared species with the species its cDNA aligns to. Random priming uses the
# identical bob-v1 3'TSO construct, so only a run-rules row is needed (this folder provides it).
#
# Prerequisites (see the main pipeline README for building the reference):
#   PIPE        = the nanopore/transfer_screen directory of this repository (has speciesmix_bobcode_metrics.py,
#                 bobcode_general_rules.md, rdna_loci.bed)
#   STARLONG    = STARlong 2.7.11b (official Dobin build)
#   STAR_INDEX  = combined GRCh38+GRCm39 STAR index (HUMAN_/MOUSE_ contig prefixes)
#   GTF         = combined_genome.gtf (HUMAN_/MOUSE_ prefixed, Ensembl 113)
# macOS Apple-silicon: run STARlong under Rosetta, e.g.
#   export DYLD_FALLBACK_LIBRARY_PATH=/opt/homebrew/Caskroom/miniforge/base/envs/star_x86/lib
#
# Reads are capped at 10,000 per sample (max_reads in the run rules) to match the deposited table.
set -euo pipefail
: "${PIPE:?set PIPE to the nanopore/transfer_screen directory of the repository}"
: "${STARLONG:?set STARLONG to STARlong 2.7.11b}"
: "${STAR_INDEX:?set STAR_INDEX to the combined STAR index dir}"
: "${GTF:?set GTF to combined_genome.gtf}"
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(cd "$HERE/../.." && pwd)"

# 1) align (one STARlong at a time) + score barcodes -> results/bobcode_metrics.tsv
python3 "$PIPE/speciesmix_bobcode_metrics.py" \
  --run-rules "$HERE/bobcode_run_rules_randompriming.tsv" \
  --fastq-dir "$ROOT/raw_fastq" \
  --out-dir "$ROOT/results_barcode_accuracy" \
  --starlong "$STARLONG" --star-index "$STAR_INDEX" --gtf "$GTF" \
  --rdna-bed "$PIPE/rdna_loci.bed" --threads 1

# 2) draw the accuracy figure (paper style; random-priming panel)
python3 "$HERE/make_accuracy_figure_random.py" \
  --results "$ROOT/results_barcode_accuracy/bobcode_metrics.tsv" \
  --run-rules "$HERE/bobcode_run_rules_randompriming.tsv" \
  --out-dir "$ROOT/figures_reference"

echo "done: results_barcode_accuracy/bobcode_metrics.tsv + figures_reference/barcode_accuracy_all_chemistries_log.{png,pdf}"

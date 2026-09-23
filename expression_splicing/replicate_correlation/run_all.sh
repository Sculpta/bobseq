#!/usr/bin/env bash
# Regenerate both figures (Pearson r and RMSE) from a gene-count matrix + annotation.
# Usage: ./run_all.sh <gene_counts.csv[.gz]> <annotation.gtf[.gz] | annotation.tsv>
set -euo pipefail

COUNTS="${1:?path to gene-count matrix (CSV/CSV.GZ)}"
ANNOT="${2:?path to Ensembl GTF or annotation TSV}"

python3 reproduce_reproducibility_figure.py --counts "$COUNTS" --annot "$ANNOT" --metric pearson --out figure_pearson
python3 reproduce_reproducibility_figure.py --counts "$COUNTS" --annot "$ANNOT" --metric rmse    --out figure_rmse

echo "Done: figure_pearson.{svg,pdf,png,tsv} and figure_rmse.{svg,pdf,png,tsv}"

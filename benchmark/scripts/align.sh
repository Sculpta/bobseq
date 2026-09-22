#!/bin/bash
# STAR with the frozen parameters of settings.sh, one single-end FASTQ.
#   align.sh <reads.fq.gz> <outdir> [prefix]
set -euo pipefail
source "$(dirname "$0")/settings.sh"
FQ=$1; OUTDIR=$2; PREFIX=${3:-}
mkdir -p "$OUTDIR"
STAR --runThreadN "$BM_THREADS" --genomeDir "$BM_REF/STAR_index" \
     --readFilesIn "$FQ" --readFilesCommand gunzip -c \
     --outFileNamePrefix "$OUTDIR/$PREFIX" $BM_STAR_PARAMS

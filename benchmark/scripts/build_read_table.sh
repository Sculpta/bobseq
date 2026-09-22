#!/bin/bash
# Canonical read table for one arm (read_table.py), built by N parallel workers over one BAM.
#   build_read_table.sh <bam> <outdir> [nworkers]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/settings.sh"
BAM=$1; OUT=$2; N=${3:-8}
export BOBSEQ_GTF="$BM_REF/combined_genome.gtf" BOBSEQ_RDNA_BED="$BM_CONFIG/rdna_loci.bed" MPLCONFIGDIR="$OUT/.mpl"
mkdir -p "$OUT/parts" "$OUT/.mpl"
for k in $(seq 0 $((N-1))); do
  python3 "$HERE/read_table.py" "$BAM" --mod "$N" --rem "$k" > "$OUT/parts/part.$k.tsv" 2> "$OUT/parts/part.$k.log" &
done
wait
cat "$OUT"/parts/*.log | grep -v -i warn || true
cat "$OUT"/parts/part.*.tsv > "$OUT/reads.tsv"
rm -f "$OUT"/parts/part.*.tsv
echo "  reads.tsv rows: $(wc -l < "$OUT/reads.tsv")"

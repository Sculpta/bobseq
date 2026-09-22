#!/bin/bash
# Per-sample splice-junction tables for the junction panels, from the per-sample BAMs of build_per_sample_bams.sh.
#   extract_junctions.sh [nparallel]
#
# Per sample (every row of benchmark_uniform/per_sample_bams/metadata.tsv):
#   1 dedup_umi_position.py     one molecule per (UMI, fragment position): UMI + 5' end + strand for single-end reads,
#                               UMI + fragment start + template length + strand for pairs (both mates kept or dropped together)
#   2 MAPQ 255                  uniquely mapped molecules only
#   3 junctions_by_fragment.py  junction = intron coordinates of every N in the CIGAR, kept when the longest left and right
#                               anchors over all spanning records are >= 4 nt (regtools' rule); support counted per read name,
#                               so a pair spanning a junction counts once
# Output: data/junctions/<set>__<sample>_junctions.bed.gz (chrom, intron start 0-based, intron end, sample, molecules, strand '?')
# and data/junctions/sample_metadata.csv (the per-sample sheet the rarefaction scripts read: set, method, sample, in_benchmark,
# records before and after deduplication and after the MAPQ filter, junctions).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/scripts/settings.sh"
NPAR=${1:-8}
PB=$BM_UNIFORM/per_sample_bams; OUT=$BM_WORK/data/junctions; TMP=$OUT/tmp; mkdir -p "$OUT" "$TMP"
one() {   # set sample
  local set=$1 s=$2 acc="$1__$2"
  [ -s "$OUT/${acc}_junctions.bed.gz" ] && [ -s "$OUT/stats/$acc.tsv" ] && return 0
  mkdir -p "$OUT/stats"
  python3 "$HERE/scripts/dedup_umi_position.py" --set "$set" --sample "$s" --in "$PB/$set/$s.bam" --out "$TMP/$acc.dedup.bam" 2> "$TMP/$acc.dedup.log"
  python3 - "$TMP/$acc.dedup.bam" "$TMP/$acc.q255.bam" <<'PY'
import sys, pysam
b = pysam.AlignmentFile(sys.argv[1], 'rb'); o = pysam.AlignmentFile(sys.argv[2], 'wb', template=b); n = m = 0
for r in b.fetch(until_eof=True):
    n += 1
    if r.mapping_quality == 255:
        o.write(r); m += 1
o.close(); print(f'{n}\t{m}')
PY
  python3 "$HERE/scripts/junctions_by_fragment.py" --bam "$TMP/$acc.q255.bam" --sample "$acc" --out "$TMP/$acc.junctions.bed" --min-anchor 4
  gzip -c "$TMP/$acc.junctions.bed" > "$OUT/${acc}_junctions.bed.gz.tmp" && mv "$OUT/${acc}_junctions.bed.gz.tmp" "$OUT/${acc}_junctions.bed.gz"
  local n_in n_dedup n_q255 nj
  read -r _ _ _ n_in n_dedup _ < "$TMP/$acc.dedup.log"   # dedup <set> <sample> <records in> <records out> <unparsed>
  n_q255=$(python3 -c "import pysam,sys; print(sum(1 for _ in pysam.AlignmentFile(sys.argv[1],'rb').fetch(until_eof=True)))" "$TMP/$acc.q255.bam")
  nj=$(awk '$1 != "."' "$TMP/$acc.junctions.bed" | wc -l)
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$acc" "$set" "$s" "$n_in" "$n_dedup" "$n_q255" "$nj" > "$OUT/stats/$acc.tsv"
  rm -f "$TMP/$acc.dedup.bam" "$TMP/$acc.dedup.bam.bai" "$TMP/$acc.q255.bam" "$TMP/$acc.junctions.bed"
  echo "  $acc: records $n_in, deduplicated $n_dedup, MAPQ 255 $n_q255, junctions $nj"
}
export -f one; export HERE OUT TMP PB PYTHONDONTWRITEBYTECODE
python3 - "$PB/metadata.tsv" <<'PY' | xargs -P "$NPAR" -L 1 bash -c 'one "$0" "$1"'
import csv, sys
for r in csv.DictReader(open(sys.argv[1]), delimiter='\t'):
    print(r['bam'].split('/')[0], r['sample'])
PY
# the sample sheet: one row per sample-arm, joined from metadata.tsv and the per-sample counts
python3 - "$PB/metadata.tsv" "$OUT" <<'PY'
import csv, sys, os
meta, out = sys.argv[1:3]
rows = []
for r in csv.DictReader(open(meta), delimiter='\t'):
    set_dir = r['bam'].split('/')[0]; acc = f"{set_dir}__{r['sample']}"
    f = open(f'{out}/stats/{acc}.tsv').read().rstrip('\n').split('\t')
    rows.append({'acc_number': acc, 'sample_name': r['sample'], 'set': set_dir, 'method': r['method'], 'cell_line': r['cell_line'],
                 'condition': r['condition'], 'replicate': r['replicate'], 'read_layout': r['read_layout'], 'in_benchmark': r['in_benchmark'],
                 'source': r['source'], 'bobcode': r.get('bobcode', ''), 'count_unit': 'mol', 'records_in': f[3], 'records_dedup': f[4],
                 'records_q255': f[5], 'junctions': f[6]})
with open(f'{out}/sample_metadata.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(f'sample_metadata.csv: {len(rows)} sample-arms')
PY
echo "JUNCTIONS DONE"

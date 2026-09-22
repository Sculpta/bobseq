#!/bin/bash
# Per-well region composition for one arm, BENCHMARK wells only: the well rules of composition_benchmark_wells.sh
# (wells_keep.txt after well_map.tsv; species from well_species.tsv if present, else HUMAN_ contigs) plus the arm's
# units_excluded.tsv (first column = well; header skipped), so BOBseq counts the 18 benchmark wells, not the 21 human wells.
# reads.tsv columns: well/barcode, UMI, contig, pos, U|M, region (E X I G R T), gene, strand.
# Emits well <TAB> region <TAB> reads_all_mapped <TAB> reads_uniquely_mapped -> <arm>/composition_per_well.tsv (panel p3b).
D=$1; W=$D/well_species.tsv; M=$D/well_map.tsv; K=$D/wells_keep.txt; X=$D/units_excluded.tsv
mawk -F'\t' -v OFS='\t' -v hasM=$([ -s "$M" ] && echo 1 || echo 0) -v hasW=$([ -s "$W" ] && echo 1 || echo 0) -v hasK=$([ -s "$K" ] && echo 1 || echo 0) '
  FILENAME==ARGV[1] { keep[$1]=1; next }
  FILENAME==ARGV[2] { if (hasM) map[$1]=$2; next }
  FILENAME==ARGV[3] { if (hasW) sp[$1]=$2; next }
  FILENAME==ARGV[4] { if (FNR > 1) exc[$1]=1; next }
  { w = (hasM ? ($1 in map ? map[$1] : "") : $1); if (hasK && !(w in keep)) next; if (w in exc) next
    if (hasW) { if ((w in sp) && ($3 !~ sp[w])) next } else if ($3 !~ /^HUMAN_/) next
    k = w SUBSEP $6; a[k]++; if ($5=="U") u[k]++ }
  END { for (k in a) { split(k, p, SUBSEP); print p[1], p[2], a[k], (k in u ? u[k] : 0) } }' \
  "$([ -s "$K" ] && echo "$K" || echo /dev/null)" "$([ -s "$M" ] && echo "$M" || echo /dev/null)" "$([ -s "$W" ] && echo "$W" || echo /dev/null)" "$([ -s "$X" ] && echo "$X" || echo /dev/null)" "$D/reads.tsv" | LC_ALL=C sort > "$D/composition_per_well.tsv"
awk -F'\t' -v n="$(basename $D)" '{w[$1]=1; a[$1]+=$3; if($2=="I") i[$1]+=$3} END{m=0; for(k in w){m++; printf "%s\t%s\tmapped %d\tintronic %.1f%%\n", n, k, a[k], 100*i[k]/a[k]}; printf "%-24s %d wells\n", n, m}' "$D/composition_per_well.tsv" | sort > "$D/composition_per_well.log"; tail -1 "$D/composition_per_well.log"

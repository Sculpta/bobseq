#!/bin/bash
# Region composition restricted to the BENCHMARK wells (wells_keep.txt after well_map.tsv), for one arm: the plate-level
# composition_species.txt counts every well that was aligned (DRUG-seq: the whole 384-well plate in the 50-nt arm).
# Same counting as composition_by_species.sh (species from well_species.tsv if present, else HUMAN_ contigs; all mapped + unique).
# Emits region <TAB> reads_all_mapped <TAB> reads_uniquely_mapped -> <arm>/composition_benchmark_wells.txt
D=$1; W=$D/well_species.tsv; M=$D/well_map.tsv; K=$D/wells_keep.txt
[ -s "$K" ] || { cp "$D/composition_species.txt" "$D/composition_benchmark_wells.txt"; echo "$(basename $D): no wells_keep, copied plate composition"; exit 0; }
mawk -F'\t' -v OFS='\t' -v hasM=$([ -s "$M" ] && echo 1 || echo 0) -v hasW=$([ -s "$W" ] && echo 1 || echo 0) '
  FILENAME==ARGV[1] { keep[$1]=1; next }
  FILENAME==ARGV[2] { if (hasM) map[$1]=$2; next }
  FILENAME==ARGV[3] { if (hasW) sp[$1]=$2; next }
  { w = (hasM ? ($1 in map ? map[$1] : "") : $1); if (!(w in keep)) next
    if (hasW) { if ((w in sp) && ($3 !~ sp[w])) next } else if ($3 !~ /^HUMAN_/) next
    a[$6]++; if ($5=="U") u[$6]++ }
  END { for (k in a) print k, a[k], (k in u ? u[k] : 0) }' "$K" "$([ -s "$M" ] && echo "$M" || echo /dev/null)" "$([ -s "$W" ] && echo "$W" || echo /dev/null)" "$D/reads.tsv" | LC_ALL=C sort > "$D/composition_benchmark_wells.txt"
awk -F'\t' -v n="$(basename $D)" '{a+=$2; if($1=="E") e=$2; if($1=="R") r=$2} END{printf "%-24s mapped %11d  exonic-pc %5.1f%%  rRNA %5.1f%%\n", n, a, 100*e/a, 100*r/a}' "$D/composition_benchmark_wells.txt"

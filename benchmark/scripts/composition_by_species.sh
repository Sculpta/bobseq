#!/bin/bash
# Region composition for one arm, each unit counted on its declared species (well_species.tsv) or,
# without that file, on HUMAN_ contigs (the competitor arms are human-only studies).
# Emits region <TAB> reads_all_mapped <TAB> reads_uniquely_mapped.
# BOTH are needed: rRNA is repetitive, so ~99% of rRNA reads are multimappers and a unique-only
# composition reports ~0.2% rRNA where the real figure is ~20-27%.
# The UCI and gene metrics use the unique+exonic subset; composition does not.
D=$1; W=$D/well_species.tsv
if [ -s "$W" ]; then
  mawk -F'\t' 'NR==FNR{sp[$1]=$2; next} (!($1 in sp)) || ($3 ~ sp[$1]) {a[$6]++; if($5=="U") u[$6]++} END{for(k in a) print k"\t"a[k]"\t"(k in u ? u[k] : 0)}' "$W" "$D/reads.tsv"
else
  mawk -F'\t' '$3 ~ /^HUMAN_/ {a[$6]++; if($5=="U") u[$6]++} END{for(k in a) print k"\t"a[k]"\t"(k in u ? u[k] : 0)}' "$D/reads.tsv"
fi | LC_ALL=C sort > "$D/composition_species.txt"
awk -F'\t' -v n="$(basename $D)" '{a+=$2; u+=$3; if($1=="R") r=$2} END{printf "  %-20s mapped %11d  rRNA %5.1f%%  unique %4.1f%% of mapped\n", n, a, 100*r/a, 100*u/a}' "$D/composition_species.txt"

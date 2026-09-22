#!/bin/bash
# The slice: canonical read table -> wupg.tsv (well umi contig pos gene) for one
# (regions, nh) choice, with the per-arm corrections applied HERE and nowhere else:
#   species     <arm>/arm.env may set SPECIES (regex on contig); default human only
#   well_map    <arm>/well_map.tsv  observed barcode -> whitelist barcode; others dropped
#   wells_keep  <arm>/wells_keep.txt one well per line; only these survive
#   slice.sh <arm_dir> <REG> <NH>
set -euo pipefail
D=$1; REG=${2:-E}; NH=${3:-U}
S="$D/stats_${REG}_${NH}"; mkdir -p "$S"
SPECIES='^HUMAN_'
[ -f "$D/arm.env" ] && source "$D/arm.env"
case "$REG" in
  E)  RX='$6=="E"' ;;
  EI) RX='($6=="E"||$6=="I") && $7!="."' ;;
  *)  echo "slice: bad regions $REG" >&2; exit 1 ;;
esac
case "$NH" in
  U)   NX='$5=="U"' ;;
  ALL) NX='1' ;;
  *)   echo "slice: bad nh $NH" >&2; exit 1 ;;
esac
TMP="$S/wupg.tmp"
mawk -F'\t' -v OFS='\t' -v sp="$SPECIES" "($RX) && ($NX) && (\$3 ~ sp) { print \$1, \$2, \$3, \$4, \$7, \$8 }" \
  "$D/reads.tsv" > "$TMP"
n0=$(wc -l < "$TMP")
if [ -s "$D/well_map.tsv" ]; then
  mawk -F'\t' -v OFS='\t' 'NR==FNR{m[$1]=$2; next} ($1 in m){$1=m[$1]; print}' "$D/well_map.tsv" "$TMP" > "$TMP.2"
  mv "$TMP.2" "$TMP"
fi
n1=$(wc -l < "$TMP")
if [ -s "$D/wells_keep.txt" ]; then
  mawk -F'\t' 'NR==FNR{k[$1]=1; next} ($1 in k)' "$D/wells_keep.txt" "$TMP" > "$TMP.2"
  mv "$TMP.2" "$TMP"
fi
n2=$(wc -l < "$TMP")
# Per-well species (BOBseq arms): each bobcode is counted on the species it was declared for,
# so a mixed library does not count human + mouse genes as one gene set and a single-species
# library does not gain cross-mapped genes of the other species. Wells without a row pass.
if [ -s "$D/well_species.tsv" ]; then
  mawk -F'\t' 'NR==FNR{sp[$1]=$2; next} (!($1 in sp)) || ($3 ~ sp[$1])' "$D/well_species.tsv" "$TMP" > "$TMP.2"
  mv "$TMP.2" "$TMP"
fi
n3=$(wc -l < "$TMP")
mv "$TMP" "$S/wupg.tsv"
{ echo "species_regex	$SPECIES"; echo "after_region_nh_species	$n0"; echo "after_well_map	$n1"; echo "after_wells_keep	$n2"; echo "after_well_species	$n3"; echo "well_species	$([ -s "$D/well_species.tsv" ] && echo yes || echo no)"
  echo "well_map	$([ -s "$D/well_map.tsv" ] && echo yes || echo no)"; echo "wells_keep	$([ -s "$D/wells_keep.txt" ] && echo "$(wc -l < "$D/wells_keep.txt") wells" || echo no)"
} > "$S/slice_provenance.tsv"
echo "  slice $REG/$NH: $n0 -> well_map $n1 -> wells_keep $n2 -> well_species $n3 rows  (species $SPECIES)"

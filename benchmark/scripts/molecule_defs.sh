#!/bin/bash
# Unique molecules under every definition in circulation, for one arm.
#
#   A  (well, UMI)                          no position at all
#   B  (well, UMI, gene)                    zUMIs / STARsolo convention
#   C  (well, UMI, contig, SAME position)   the agreed rule, exact UMI
#   D  (well, UMI ed<=MAXED, contig, SAME pos) the agreed rule plus sequencing error; MAXED per arm (arm.env)
#   B_corr  B with the skew-aware collision correction (umi_b.py)  [THE COMMON CURRENCY of the benchmark]
#   E  (well, UMI ed<=1, contig, +/-1 kb)   the retired window, kept to show its effect
#
#   molecule_defs.sh <arm_dir> <REG> <NH>
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/settings.sh"
D=$1; REG=${2:-E}; NH=${3:-U}
# an arm may override the UMI edit distance (arm.env: UMI_MAXED=0 for a 6 nt UMI, where
# ed<=1 over-merges at hotspots: 88% vs 97% planted-molecule recovery on the synthetic 24-plex)
[ -f "$D/arm.env" ] && source "$D/arm.env"
BM_UMI_MAXED=${UMI_MAXED:-$BM_UMI_MAXED}
S="$D/stats_${REG}_${NH}"; W="$S/basic"; mkdir -p "$W"
bash "$HERE/slice.sh" "$D" "$REG" "$NH"
if ! mkdir "$W/.lock" 2>/dev/null; then echo "  locked, another run is here" >&2; exit 3; fi
trap 'rmdir "$W/.lock" 2>/dev/null' EXIT
SORT="LC_ALL=C sort -S 12G --parallel=8 -T $W"

R=$(wc -l < "$S/wupg.tsv"); echo "  reads in slice: $R"
echo "  [A] (well, UMI)                        $(date -u +%H:%M:%S)"
A=$(eval $SORT -u -k1,1 -k2,2 "$S/wupg.tsv" | wc -l)
echo "  [B] (well, UMI, gene)                  $(date -u +%H:%M:%S)"
python3 "$HERE/umi_b.py" "$D" "$REG" "$NH" > "$W/summary_B.tsv"
B=$(mawk -F'	' '$1=="mol_B_raw"{print $2}' "$W/summary_B.tsv")
echo "  sorting by position                    $(date -u +%H:%M:%S)"
eval $SORT -k1,1 -k3,3 -k6,6 -k4,4n "$S/wupg.tsv" > "$W/by_pos.tsv"
echo "  [C] same position, exact UMI           $(date -u +%H:%M:%S)"
C=$("$HERE/dedup" 0 0 /dev/null < "$W/by_pos.tsv" 2>"$W/dedup_C.log" | wc -l)
echo "  [D] same position, UMI ed<=$BM_UMI_MAXED (per arm)  $(date -u +%H:%M:%S)"
"$HERE/dedup" "$BM_UMI_MAXED" "$BM_POS_TOL" "$W/perwell_D.tsv" "$W/hist_D.tsv" < "$W/by_pos.tsv" > "$W/mol_D.tsv" 2>"$W/dedup_D.log"
Dn=$(wc -l < "$W/mol_D.tsv")
echo "  [E] +/-1 kb window, UMI ed<=1 (retired) $(date -u +%H:%M:%S)"
E=$("$HERE/dedup" "$BM_UMI_MAXED" 1000 /dev/null < "$W/by_pos.tsv" 2>"$W/dedup_E.log" | wc -l)

echo "  per-well reads and genes               $(date -u +%H:%M:%S)"
mawk -F'\t' '{c[$1]++} END{for(w in c) print w"\t"c[w]}' "$S/wupg.tsv" \
  | LC_ALL=C sort -k1,1 > "$W/reads_per_well.tsv"
eval $SORT -u -k1,1 -k2,2 "$W/mol_D.tsv" \
  | mawk -F'\t' '{g[$1]++} END{for(w in g) print w"\t"g[w]}' \
  | LC_ALL=C sort -k1,1 > "$W/genes_per_well.tsv"

{ echo -e "metric\tvalue"
  echo -e "mol_A_well_umi\t$A"
  echo -e "mol_B_well_umi_gene\t$B"
  echo -e "mol_C_samepos_exact\t$C"
  echo -e "mol_D_samepos_ed1	$Dn"
  echo -e "mol_D_samepos_ed	$Dn"
  echo -e "umi_maxed	$BM_UMI_MAXED"
  echo -e "pos_tol	$BM_POS_TOL"
  cat "$W/summary_B.tsv"
  echo -e "mol_E_1kb_ed1\t$E"
  echo -e "slice_reads\t$R"
  echo -e "n_wells\t$(wc -l < "$W/reads_per_well.tsv")"
} > "$W/summary.tsv"
rm -f "$W/by_pos.tsv"
# ordering sanity: the definitions nest; a violation is a bug, not biology
Bc=$(mawk -F'	' '$1=="mol_B_corr"{print $2}' "$W/summary_B.tsv")
[ "$A" -le "$B" ]  || echo "  SANITY: A ($A) > B ($B)" >&2
Bk=$(mawk -F'\t' '$1=="mol_B_kept"{print $2}' "$W/summary_B.tsv")
[ "$Bc" -ge "$Bk" ] || echo "  SANITY: B_corr ($Bc) < B_kept ($Bk)" >&2   # correction applies AFTER failed-call exclusion
[ "$C" -ge "$Dn" ] || echo "  SANITY: C ($C) < D ($Dn)" >&2
[ "$Dn" -ge "$E" ] || echo "  SANITY: D ($Dn) < E ($E)" >&2
sed 's/^/    /' "$W/summary.tsv"
echo "  done: $W/summary.tsv"

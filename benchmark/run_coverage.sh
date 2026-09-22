#!/bin/bash
# Per-sample coverage inputs of the panels, from the per-sample BAM store: canonical-transcript read positions, per-sample
# junction tables at matched depth, read-length and insert-length metrics, the Picard refFlats, and the wells outside the
# benchmark set (the three risdiplam 500 nM wells, the three mouse wells) through the same path.
#   run_coverage.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/scripts/settings.sh"
U=$BM_UNIFORM; PB=$U/per_sample_bams; A=$BM_COVERAGE; DATA=$BM_WORK/data; L=$BM_WORK/logs/coverage; mkdir -p "$L" "$A" "$BM_WORK/results/library_metrics/insert_per_well" "$BM_WORK/results/supplement_21"
log() { echo "[$(date '+%F %T')] $*"; }
export BM_COVSET

log "refFlats for Picard: protein-coding transcripts without ribosomal-protein and mitochondrial genes, bare contig names"
[ -s "$A/picard_refflat_pc_noRP_noMT_bare.txt" ] || python3 "$HERE/scripts/build_refflat_filtered.py" HUMAN_ "$A/picard_refflat_pc_noRP_noMT_bare.txt"
[ -s "$A/picard_refflat_pc_noRP_noMT_bare_mouse.txt" ] || python3 "$HERE/scripts/build_refflat_filtered.py" MOUSE_ "$A/picard_refflat_pc_noRP_noMT_bare_mouse.txt"

for SET in native 50nt; do
  BM_COVSET=$SET; SFX=$([ "$SET" = native ] && echo "" || echo "_50nt")
  POS=$A/positions_canonical$SFX; JX=$A/junctions$SFX; mkdir -p "$POS" "$JX"
  N=$(wc -l < "$A/samples$SFX.tsv")
  log "[$SET] read positions on canonical transcripts, $N samples"
  cut -f1,2 "$A/samples$SFX.tsv" | xargs -P "$BM_THREADS" -L 1 bash -c '[ -s "'$POS'/$0.npz" ] || python3 "'$HERE'/scripts/build_positions.py" "$0" "'$PB'/$1" "'$POS'/$0.npz"' 2>&1 | grep -v Warn || true
  [ "$(ls "$POS"/*.npz | wc -l)" -ge "$N" ] || { log "[$SET] positions incomplete"; exit 1; }
  log "[$SET] junctions per sample at matched depth"
  cut -f1,2 "$A/samples$SFX.tsv" | xargs -P "$BM_THREADS" -L 1 bash -c '[ -s "'$JX'/$0.tsv" ] || python3 "'$HERE'/scripts/junctions_per_sample.py" "$0" "'$PB'/$1" "'$JX'/$0.tsv"' 2>&1 | grep -v Warn || true
  log "[$SET] aligned read length per method and sample"
  python3 "$HERE/scripts/library_metrics.py" readlen "$SET" > "$L/readlen_$SET.log" 2>&1
done

log "insert length per well, deduplicated pair BAMs of the pipeline"
O=$BM_WORK/results/library_metrics/insert_per_well
ls "$DATA/bams/pipeline_dedup_pe/dedup_pairs_6NHH"/*_star_pairs.bam | xargs -P "$BM_THREADS" -I{} bash -c 's=$(basename {} _star_pairs.bam); [ -s "'$O'/$s.json" ] || python3 "'$HERE'/scripts/library_metrics.py" insert "$s" {} "'$O'/$s.json" > "'$O'/$s.log" 2>&1'
grep -h "insert median" "$O"/*.log | sort | cut -c1-170

log "Ris 500 nM wells (input failures) through the arm path, native set"
python3 "$HERE/scripts/supp21_build_ris500.py" > "$L/supp21_build.log" 2>&1
for w in Ris_dose_500mM-1 Ris_dose_500mM-2 Ris_dose_500mM-3; do bash "$HERE/scripts/build_read_table.sh" "$U/bob57_24plex_pe_native/supp21/$w/Aligned.out.bam" "$U/bob57_24plex_pe_native/supp21/$w" 4 2>&1 | tail -n 1; done
python3 "$HERE/scripts/supp21_tables.py" native > "$L/supp21_tables.log" 2>&1; tail -n 2 "$L/supp21_tables.log"
BM_COVSET=native python3 "$HERE/scripts/supp21_native_bams.py" > "$L/supp21_native_bams.log" 2>&1; tail -n 1 "$L/supp21_native_bams.log"

log "mouse wells (RAW 264.7) through the arm path"
python3 "$HERE/scripts/mouse3_build_wells.py" > "$L/mouse3_build.log" 2>&1; grep RAW "$L/mouse3_build.log"
for w in RAW_control_1 RAW_control_2 RAW_control_3; do bash "$HERE/scripts/build_read_table.sh" "$U/bob57_24plex_pe_native/mouse3/$w/Aligned.out.bam" "$U/bob57_24plex_pe_native/mouse3/$w" 4 2>&1 | tail -n 1; done
python3 "$HERE/scripts/mouse3_tables.py" > "$L/mouse3_tables.log" 2>&1; grep RAW "$L/mouse3_tables.log"
python3 "$HERE/scripts/rarefy.py" "$U/bob57_24plex_pe_native" E_U_mouse3 > "$L/mouse3_rarefy.log" 2>&1; tail -n 1 "$L/mouse3_rarefy.log"   # matched-depth table of the mouse wells (stacked QC supplement)
for s in picard composition positions insert; do python3 "$HERE/scripts/mouse3_native_bams.py" $s > "$L/mouse3_$s.log" 2>&1 || { log "mouse3 $s FAILED"; tail -3 "$L/mouse3_$s.log"; exit 1; }; done
log "COVERAGE DONE"

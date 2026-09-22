#!/bin/bash
# Barcode accuracy of the 24-plex run: the strict per-well funnel on the pipeline's deduplicated analysis BAMs.
#   run_funnel.sh <R2.fastq.gz>
#   <R2.fastq.gz>  read 2 of the lane as sequenced (the bobcode is at its start); the exact-bobcode pass reads it once
# Steps: annotation index from the combined GTF and the rDNA loci -> exact-bobcode read keys from the raw read 2 -> one funnel per
# well on <well>_star_combined.bam (deduplicated molecules) -> funnel_summary.json and BARCODE_ACCURACY_FUNNEL.md.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/scripts/settings.sh"
R2=$1
F=$BM_WORK/results/barcode_accuracy_funnel; D=$BM_WORK/data/bams/pipeline_dedup_pe/dedup_analysis_6NHH; mkdir -p "$F/inputs" "$F/per_well"
log() { echo "[$(date '+%F %T')] $*"; }
cd "$F"
[ -s inputs/annotation_index.npz ] || { log "annotation index"; python3 "$HERE/scripts/funnel_build_annotation_index.py" inputs/annotation_index.npz > inputs/build_annotation_index.log; }
[ -s inputs/exact_code_keys.npz ] || { log "exact bobcode calls from the raw read 2"; python3 "$HERE/scripts/funnel_exact_code_flags.py" "$R2" inputs/exact_code_keys.npz > inputs/exact_code_flags.log; tail -n 1 inputs/exact_code_flags.log; }
n=$(ls "$D"/*_star_combined.bam | wc -l); log "$n deduplicated analysis BAMs"; [ "$n" = 24 ] || { echo "expected 24 wells" >&2; exit 1; }
ls "$D"/*_star_combined.bam | xargs -P "$BM_THREADS" -I{} bash -c 's=$(basename {} _star_combined.bam); [ -s "'$F'/per_well/$s.json" ] || python3 "'$HERE'/scripts/funnel_per_well.py" "$s" {} "'$F'/inputs/annotation_index.npz" "'$F'/inputs/exact_code_keys.npz" "'$F'/per_well/$s.json" > "'$F'/per_well/$s.log" 2>&1'
log "per-well done: $(ls per_well/*.json | wc -l) wells"; grep -l Traceback per_well/*.log | sed 's/^/TRACEBACK /' || true
cd "$F" && python3 "$HERE/scripts/funnel_report.py"
log "FUNNEL DONE: $F/BARCODE_ACCURACY_FUNNEL.md"

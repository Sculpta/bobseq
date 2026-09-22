#!/bin/bash
# The whole benchmark, in order, from the raw data to the figure tree. Every stage is a script of this directory that can be
# run on its own; each skips outputs that already exist.
#   run_benchmark.sh <pipeline_out_dir> <R1.fastq.gz> <R2.fastq.gz>
#
#   <pipeline_out_dir>  the output of pipeline/run_illumina.sh on the 24-plex library (config/run_24plex.json)
#   <R1> <R2>           the same library as sequenced (the lane files the pipeline was run on)
# Settings: benchmark/settings.sh (BOBSEQ_WORK for the working directory, BOBSEQ_REF_DIR for the reference, BOBSEQ_THREADS).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/scripts/settings.sh"
PIPE=$1; R1=$2; R2=$3
mkdir -p "$BM_WORK/logs"; L=$BM_WORK/logs
log() { echo "[$(date '+%F %T')] $*"; }
stage() { local name=$1; shift; log "==== $name"; "$@" > "$L/$name.log" 2>&1 || { log "$name FAILED, see $L/$name.log"; tail -n 20 "$L/$name.log"; exit 1; }; tail -n 1 "$L/$name.log"; }

[ -x "$HERE/scripts/dedup" ] || { log "compiling dedup.c"; gcc -O2 -o "$HERE/scripts/dedup" "$HERE/scripts/dedup.c"; }
stage 01_fetch_competitors    bash "$HERE/fetch_competitors.sh"
stage 02_competitor_arms      bash "$HERE/build_arms.sh"
stage 03_bobseq_arms          bash "$HERE/build_bobseq_arms.sh" "$PIPE" "$R1" "$R2"
stage 04_per_sample_bams      bash "$HERE/build_per_sample_bams.sh"
stage 05_junctions            bash "$HERE/extract_junctions.sh" "$BM_THREADS"
stage 06_coverage             bash "$HERE/run_coverage.sh"
stage 07_barcode_accuracy     bash "$HERE/run_funnel.sh" "$R2"
stage 08_panels               bash "$HERE/run_panels.sh"
log "BENCHMARK DONE: figures in $BM_WORK/figures, values tables beside every panel"

#!/bin/bash
# The figure panels, the values tables behind them, the curated figure tree and the validators.
#   run_panels.sh
# Runs after build_arms.sh, build_bobseq_arms.sh, build_per_sample_bams.sh, extract_junctions.sh, run_coverage.sh and run_funnel.sh.
# BM_COVSET selects the set: native (each method at its own read length, the preprint's figures) or 50nt (every method reduced to
# one 50-nt read; only its poly(A)-distance panel is used, as a control). Generator output goes to preprint_figures_<set>/;
# organize_figures.py then rebuilds figures/{main,supplement,...} from it.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/scripts/settings.sh"
L=$BM_WORK/logs/panels; mkdir -p "$L" "$BM_WORK/figures/_generator_output"
ln -sfn "$BM_WORK/preprint_figures_native" "$BM_WORK/figures/_generator_output/native"
ln -sfn "$BM_WORK/preprint_figures_50nt" "$BM_WORK/figures/_generator_output/50nt"
log() { echo "[$(date '+%F %T')] $*"; }
run() {   # <set> <script> [args]: one generator, log per set
  local s=$1 sc=$2; shift 2
  BM_COVSET=$s python3 "$HERE/scripts/$sc" "$@" > "$L/${s}_$(basename "$sc" .py).log" 2>&1 || { log "[$s] FAILED $sc (see $L)"; tail -n 5 "$L/${s}_$(basename "$sc" .py).log"; exit 1; }
  log "[$s] $sc"
}
GENES="RELA MLF2 BAG6 TUBG1 EEF2 MCM7 GIGYF1 GRB2 IRS4 FUS PPP5C"

for s in native 50nt; do
  run $s composition_rules.py --par "$BM_THREADS"     # read composition by the preprint's rules, per sample (p3 values)
  run $s picard_profile.py                             # Picard CollectRnaSeqMetrics per sample, 5'-3' balance (p6, p6b)
  SFX=$([ "$s" = native ] && echo "" || echo "_$s")
  # coverage architecture along the mature transcript (sections A-K, J, V, P, Q) from the position tables of run_coverage.sh;
  # main_figure_panels.py reads A_breadth_vs_depth.json and V2_fraction_covered_vs_depth.json, validate_main_panels.py Q_UNSCALED.md
  run $s coverage_architecture.py "$BM_COVERAGE/samples$SFX.tsv" "$BM_COVERAGE/positions_canonical$SFX" "$BM_WORK/preprint_figures_$s/coverage_architecture"
  cp "$BM_COVERAGE/samples$SFX.tsv" "$BM_WORK/preprint_figures_$s/coverage_architecture/"
  run $s polya_region_genes.py
  run $s coverage_v2.py
  run $s picard_custom.py
  run $s picard_unscaled.py
  for sc in sample_table_panel.py read_fate_panels.py main_figure_panels.py composition_panels.py threshold_panels.py library_panels.py; do run $s $sc; done
done
run native panel_legends.py
log "junction panels"
for sc in junction_rarefaction.py junction_rarefaction_annotated.py junction_rarefaction_pooled.py; do run native $sc; done
# tidy the main panels and write their values tables now: the supplements below read them (p0a per-sample read fate, dedup per sample)
for s in native 50nt; do run $s tidy_panels.py; run $s panel_values.py; run $s panel_readme.py; done
log "per-sample supplement: the 21 human wells, the mouse wells, the stacked figure"
run native supplement_per_sample.py
run native supplement_per_sample_mouse.py
run native supplement_polyA_24.py
run native stacked_supplement_24plex.py                 # the two separate stacked figures (24-plex, benchmark) kept for reference
run native stacked_supplement_benchmark.py
run native stacked_supplement_merged.py
log "pooled depth-matched tracks and the gene coverage examples"
for step in hash write finish; do run native ucsc_pooled_tracks.py $step; done
run native gene_coverage_tracks.py $GENES
log "prime-seq peak check (methods text)"
run native gene_peakiness_scan.py
run native primeseq_peak_sequences.py
run native primeseq_peak_umis.py
log "values tables, panel notes, sources"
run native panel_sources.py
log "curated figure tree"
run native organize_figures.py
log "validation"
# independent recount of the genes per well without ribosomal-protein genes (rebuilt reservoir, regex-free rule); the main-panel validator compares it
log "independent gene recount per arm"
echo "drugseq:E_U primeseq:E_U bob57_24plex_pe_hs:E_U drugseq_native:E_U primeseq_native:E_U bob57_24plex_pe_native:E_U bob57_24plex_pe_native:E_U_supp21" | tr ' ' '\n' | \
  xargs -P 7 -I{} bash -c 'j={}; a=${j%%:*}; s=${j##*:}; python3 "'$HERE'/scripts/check_genes_noRP.py" "'$BM_UNIFORM'/$a" $s > "'$L'/check_genes_noRP_$a.$s.log" 2>&1; echo "  $a $s exit $?"'
for s in native 50nt; do run $s report_numbers.py; done              # per-arm numbers the validator checks the panels against
for s in native 50nt; do BM_COVSET=$s python3 "$HERE/scripts/validate_main_panels.py" > "$L/${s}_validate_main.log" 2>&1 || true; log "[$s] main panels: $(grep 'FAILED of' "$L/${s}_validate_main.log")"; done
BM_COVSET=native python3 "$HERE/scripts/validate_library_panels.py" > "$L/validate_library.log" 2>&1 || true; log "library panels: $(grep -c '^PASS' "$L/validate_library.log") PASS, $(grep -c '^FAIL' "$L/validate_library.log") FAIL"
BM_COVSET=native python3 "$HERE/scripts/validate_new_outputs.py" > "$L/validate_new_outputs.log" 2>&1 || true; log "new outputs: $(grep -c '^PASS' "$L/validate_new_outputs.log") PASS, $(grep -c '^FAIL' "$L/validate_new_outputs.log") FAIL"
log "PANELS DONE: $BM_WORK/figures"

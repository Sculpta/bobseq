#!/bin/bash
# Smoke test of the BOBseq Illumina processing on the bundled example (the first 200,000 read pairs of the 24-plex library).
#   bash tests/run_smoke_test.sh            prep + demultiplexing (Python only, about a minute); counts compared with expected_counts.tsv
#   bash tests/run_smoke_test.sh --full <reference_dir>   the whole pipeline (STAR per sample, deduplication, tables, report); needs the
#                                            reference of reference/build_reference.sh and the conda environment (environment.yml)
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(dirname "$HERE"); OUT=${SMOKE_OUT:-$HERE/smoke_out}
R1=$HERE/example_data/example_200k_R1.fastq.gz; R2=$HERE/example_data/example_200k_R2.fastq.gz; RUN_JSON=$ROOT/config/run_24plex.json
export PYTHONDONTWRITEBYTECODE=1 BOBCODE_REF_DIR=$ROOT/config BOBSEQ_NEB_INDEX_FILE=$ROOT/config/nebnext_e7600_dual_index_primers.tsv
mkdir -p "$OUT"
if [ "${1:-}" = "--full" ]; then
  bash "$ROOT/pipeline/run_illumina.sh" "$R1" "$R2" "$RUN_JSON" "$2" "$OUT" "${3:-8}"
else
  RUN=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['run'])" "$RUN_JSON")
  python3 "$ROOT/pipeline/scripts/make_guide.py" "$RUN_JSON" -o "$OUT" > "$OUT/make_guide.log"
  python3 "$ROOT/pipeline/scripts/prep_reads_illumina.py" "$R1" "$R2" "$OUT/prepped.fastq.gz" --r1-out "$OUT/prepped_R1.fastq.gz" --guide "$OUT/machine_readable_guide/${RUN}_analysisguide.txt" --stats "$OUT/prep_stats.json" > "$OUT/prep.log" 2>&1
  python3 "$ROOT/pipeline/scripts/demux_qname_illumina.py" "$OUT/prepped.fastq.gz" - "$OUT/demux" --run-json "$RUN_JSON" --r1 "$OUT/prepped_R1.fastq.gz" > "$OUT/demux.log" 2>&1
fi
python3 - "$OUT" "$HERE/expected_counts.tsv" <<'PY'
import json, sys
out, expected = sys.argv[1:3]
s = json.load(open(f'{out}/prep_stats.json')); d = json.load(open(f'{out}/demux/library_demux_stats.json'))
got = {('prep', 'pairs_in'): s['n'], ('prep', 'bobcode_called'): s['bobcode_called'], ('prep', 'exact_bobcode'): s['mm0'], ('prep', 'one_mismatch_bobcode'): s['mm1'],
       ('prep', 'pairs_written'): s['written'], ('demux', 'reads_assigned'): d['reads_assigned'], ('demux', 'reads_unclassified'): d['reads_unclassified']}
for code, v in d['per_bobcode'].items():
    got[('demux', 'reads_' + v['sample'].replace(' ', '_').replace('/', '_'))] = v['reads']
fails = 0
for l in open(expected).read().splitlines()[1:]:
    step, metric, exp = l.split('\t')
    g = got.get((step, metric))
    ok = g is not None and int(g) == int(exp)
    fails += not ok
    print(f"{'PASS' if ok else 'FAIL'} {step} {metric}: {g} (expected {exp})")
print(f'{fails} FAILED of {len(got)} checks')
sys.exit(1 if fails else 0)
PY

#!/usr/bin/env python3
"""Supplement (21 human wells): Picard CollectRnaSeqMetrics and the rule-based read composition for the three Ris 500 nM
wells on their native per-sample BAMs (per_sample_bams/bobseq_pe_native), exactly as picard_profile.py and composition_rules.py do
for the 18 benchmark wells. Picard metrics go to the same per-sample cache (picard_metrics/native/<sample>.RNA_Metrics.txt); the
composition counts to results/supplement_21/composition_rules_ris500_native.json."""

from settings import WORK, PROJECT, COVERAGE, REFFLAT
import os, sys, json, subprocess

H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, H)
P = WORK
PB = f'{P}/benchmark_uniform/per_sample_bams'
MET = f'{P}/preprint_figures_native/coverage_architecture/picard_metrics/native'
os.makedirs(MET, exist_ok=True)   # run_coverage.sh runs before the panel stage that otherwise creates the metrics cache
REFFLAT = f'{COVERAGE}/picard_refflat_pc_noRP_noMT_bare.txt'
PICARD = 'picard'
WELLS = ['Ris_dose_500mM-1', 'Ris_dose_500mM-2', 'Ris_dose_500mM-3']
for n in WELLS:
    out = f'{MET}/{n}.RNA_Metrics.txt'
    if os.path.exists(out) and '## HISTOGRAM' in open(out).read():
        print(n, 'picard cached', flush=True)
        continue
    cmd = [
        PICARD,
        'CollectRnaSeqMetrics',
        '-I',
        f'{PB}/bobseq_pe_native/{n}.bam',
        '-O',
        f'{MET}/{n}.RNA_Metrics.txt',
        '--REF_FLAT',
        REFFLAT,
        '--STRAND_SPECIFICITY',
        'NONE',
        '--MINIMUM_LENGTH',
        '1000',
        '--VALIDATION_STRINGENCY',
        'LENIENT',
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    open(f'{MET}/{n}.log', 'w').write(r.stdout + r.stderr)
    print(n, 'picard', 'done' if r.returncode == 0 else 'FAILED', flush=True)
os.environ['BM_COVSET'] = 'native'
import composition_rules as cr

res = {}
for n in WELLS:
    m, name, cnt = cr.classify_bam(('BOBseq', n, f'bobseq_pe_native/{n}.bam', True))
    res[name] = cnt
    print(name, cnt, flush=True)
json.dump(res, open(f'{PROJECT}/results/supplement_21/composition_rules_ris500_native.json', 'w'), indent=1)

# canonical-transcript read positions of the three wells (native set), read by the per-sample supplement and the 24-well
# poly(A) profiles as positions_canonical/bobseq-ris-500mM-<k>.npz, exactly as mouse3_native_bams.py does for the mouse wells
os.makedirs(f'{COVERAGE}/positions_canonical', exist_ok=True)
procs = []
for k, n in enumerate(WELLS, 1):
    out = f'{COVERAGE}/positions_canonical/bobseq-ris-500mM-{k}.npz'
    if os.path.exists(out):
        print('positions', k, 'cached', flush=True)
        continue
    procs.append((k, subprocess.Popen(['python3', f'{H}/build_positions.py', f'bobseq-ris-500mM-{k}', f'{PB}/bobseq_pe_native/{n}.bam', out],
                                      stdout=open(f'{COVERAGE}/positions_canonical/bobseq-ris-500mM-{k}.log', 'w'), stderr=subprocess.STDOUT)))
for k, pr in procs:
    print('positions', k, 'exit', pr.wait(), open(f'{COVERAGE}/positions_canonical/bobseq-ris-500mM-{k}.log').read().strip()[-160:], flush=True)
print('SUPP21 NATIVE BAMS DONE')

#!/usr/bin/env python3
"""Assemble preprint_figures_<set>/report_data/numbers.json from the benchmark tables: per arm the STAR counts, the slice, the
molecule definitions, the plate composition, the rarefaction medians per depth and the native per-well medians. Read by
validate_main_panels.py (slice reads, sample numbers, medians per depth). BM_COVSET selects the set: native (each method at
its own read length) or 50nt (every method as one 50-nt read)."""
from settings import WORK
import csv, json, os, statistics
import numpy as np

SET = os.environ.get('BM_COVSET', 'native')
U = f'{WORK}/benchmark_uniform'
FIG = 'preprint_figures_native' if SET == 'native' else f'preprint_figures_{SET}'
if SET == 'native':
    ARMS = [('drugseq_native', 'DRUG-seq'), ('primeseq_native', 'prime-seq'), ('bob57_24plex_pe_native', 'BOBseq')]
else:
    ARMS = [('drugseq', 'DRUG-seq'), ('primeseq', 'prime-seq'), ('bob57_24plex_pe_hs', 'BOBseq')]


def kv(p):
    d = {}
    for l in open(p):
        f = l.rstrip('\n').split('\t')
        if len(f) >= 2:
            d[f[0]] = f[1]
    return d


def logf(p):
    d = {}
    for l in open(p):
        if '|' in l:
            a, b = l.split('|', 1)
            d[a.strip()] = b.strip()
    return d


out = {'arms': []}
for k, lab in ARMS:
    s = kv(f'{U}/{k}/stats_E_U/basic/summary.tsv')
    # the selected-samples log when STAR ran on more samples than the benchmark uses
    L = logf(f'{U}/{k}/Log.final.selected.out') if os.path.exists(f'{U}/{k}/Log.final.selected.out') else logf(f'{U}/{k}/Log.final.out')
    comp = {}
    for l in open(f'{U}/{k}/composition_species.txt'):
        f = l.rstrip('\n').split('\t')
        if len(f) == 3 and len(f[0]) == 1:
            comp[f[0]] = (int(f[1]), int(f[2]))
    tot = sum(v[0] for v in comp.values())
    rows = [r for r in csv.DictReader(open(f'{U}/{k}/stats_E_U/rarefied.tsv'), delimiter='\t') if r.get('at_native', '0') == '0']
    depths = {}   # grid rows only, as the figures
    for dep in sorted({int(float(r['depth'])) for r in rows}):
        rr = [r for r in rows if int(float(r['depth'])) == dep]
        depths[dep] = {'n': len(rr), 'uci': statistics.median(float(r['B_corr']) for r in rr),
                       'D': statistics.median(float(r['molecules']) for r in rr), 'genes': statistics.median(float(r['genes_B']) for r in rr)}
    pw = {r['well']: r for r in csv.DictReader(open(f'{U}/{k}/stats_E_U/basic/perwell_B.tsv'), delimiter='\t')}
    rpw = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(f'{U}/{k}/stats_E_U/basic/reads_per_well.tsv')}
    native = {'n': len(pw), 'reads_med': statistics.median(rpw.values()), 'uci_med': statistics.median(float(v['B_corr']) for v in pw.values()),
              'genes_med': statistics.median(float(v['genes_B']) for v in pw.values()),
              'reads_iqr': (float(np.percentile(list(rpw.values()), 25)), float(np.percentile(list(rpw.values()), 75)))}
    out['arms'].append({
        'key': k, 'label': lab, 'reads_in': int(L['Number of input reads']),
        'mapped': int(L['Uniquely mapped reads number']) + int(L['Number of reads mapped to multiple loci']),
        'unique': int(L['Uniquely mapped reads number']), 'slice_reads': int(s['slice_reads']),
        'mol': {m: int(s[m]) for m in ('mol_A_well_umi', 'mol_B_well_umi_gene', 'mol_B_raw', 'mol_B_kept', 'mol_B_corr', 'mol_C_samepos_exact', 'mol_D_samepos_ed', 'mol_E_1kb_ed1') if m in s},
        'log': {'too_short_pct': float(L.get('% of reads unmapped: too short', '0%').rstrip('%'))},
        'umi_excluded_pct': float(s['umi_rows_excluded_pct']), 'n_samples': int(s['n_wells']),
        'composition': {c: round(100 * comp[c][0] / tot, 1) for c in 'EXIGRT' if c in comp},
        'composition_n': {c: comp[c][0] for c in 'EXIGRT' if c in comp}, 'mapped_all': tot, 'depths': depths, 'native': native})
excl = f'{U}/' + ('bob57_24plex_pe_native' if SET == 'native' else 'bob57_24plex_pe_hs') + '/units_excluded.tsv'
out['excluded'] = list(csv.DictReader(open(excl), delimiter='\t')) if os.path.exists(excl) else []
os.makedirs(f'{WORK}/{FIG}/report_data', exist_ok=True)
json.dump(out, open(f'{WORK}/{FIG}/report_data/numbers.json', 'w'), indent=1)
print('wrote numbers.json', SET)

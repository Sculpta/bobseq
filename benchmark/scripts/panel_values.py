#!/usr/bin/env python3
"""Raw values behind the kept main-figure panels, one TSV per panel in <fig set>/main_figure_panels/values/, for replotting from
the numbers. The threshold and Picard generators write their own TSVs there; this script covers p0a, p1, p2, p3 and p3b from the panel
sidecars and the rarefaction tables. BM_COVSET=native|50nt."""

from settings import WORK
import os, csv, json

P = WORK
U = f'{P}/benchmark_uniform'
SET = os.environ.get('BM_COVSET', 'native')
TAG = 'native' if SET == 'native' else 'matched'
D = f'{P}/{"preprint_figures_native" if SET == "native" else f"preprint_figures_{SET}"}/main_figure_panels'
V = f'{D}/values'
os.makedirs(V, exist_ok=True)
import os as _os, sys as _sys

_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from palette import METH

ARMS = {
    'native': {
        'DRUG-seq': 'drugseq_native',
        'BRB-seq': 'brbseq_native',
        'prime-seq': 'primeseq_native',
        'BOBseq': 'bob57_24plex_pe_native',
    },
    'matched': {'DRUG-seq': 'drugseq', 'BRB-seq': 'brbseq', 'prime-seq': 'primeseq', 'BOBseq': 'bob57_24plex_pe_hs'},
}[TAG]


def load(name):
    return json.load(open(f'{D}/{name}'))


# p0a read fate: per method summary + per well
p0 = load(f'p0a_read_fate_{TAG}.json')
pw = p0.pop('_perwell', {})
p0.pop('_arms', None)
with open(f'{V}/p0a_read_fate_{TAG}.tsv', 'w') as f:
    keys = [k for k in p0[METH[0]] if k != 'n_samples']
    f.write('method\tn_samples\t' + '\t'.join(keys) + '\n')
    for m in METH:
        f.write(f'{m}\t{p0[m].get("n_samples", "")}\t' + '\t'.join(str(p0[m].get(k, '')) for k in keys) + '\n')
if pw:
    with open(f'{V}/p0a_read_fate_{TAG}_per_sample.tsv', 'w') as f:
        first = next(iter(next(iter(pw.values())).values()))
        keys = list(first) if isinstance(first, dict) else ['value']
        f.write('method\tsample\t' + '\t'.join(keys) + '\n')
        for m in METH:
            for w, v in pw.get(m, {}).items():
                f.write(
                    f'{m}\t{w}\t'
                    + ('\t'.join(str(v.get(k, '')) for k in keys) if isinstance(v, dict) else str(v))
                    + '\n'
                )
# p1 / p2: molecules (UCI) and genes per sample vs subsampled depth, as plotted (at_native == 0 rows)
with open(f'{V}/p1_molecules_vs_depth_{TAG}.tsv', 'w') as f1, open(f'{V}/p2_genes_vs_depth_{TAG}.tsv', 'w') as f2:
    f1.write('method\tsample\tdepth_filtered_reads\tmolecules_UCI\n')
    f2.write('method\tsample\tdepth_filtered_reads\tgenes\tgenes_incl_ribosomal_protein\n')
    for m in METH:
        for r in csv.DictReader(open(f'{U}/{ARMS[m]}/stats_E_U/rarefied.tsv'), delimiter='\t'):
            if r.get('B_corr') and r['at_native'] == '0':
                f1.write(f'{m}\t{r["well"]}\t{r["depth"]}\t{r["B_corr"]}\n')
                f2.write(f'{m}\t{r["well"]}\t{r["depth"]}\t{r["genes_B_noRP"]}\t{r["genes_B"]}\n')
# p1 / p2 summary per depth: what the thick line and the ribbon show (mean over wells, 95% t-interval of the mean, depths with >= 3 wells)
import numpy as np
from scipy import stats as _st

for name, ycol in (('p1_molecules_vs_depth', 'molecules_UCI'), ('p2_genes_vs_depth', 'genes')):
    rows = list(csv.DictReader(open(f'{V}/{name}_{TAG}.tsv'), delimiter='\t'))
    with open(f'{V}/{name}_{TAG}_mean_ci95.tsv', 'w') as f:
        f.write(f'method\tdepth_filtered_reads\tn_samples\tmean\tsd\tci95_halfwidth\tci95_lo\tci95_hi\n')
        for m in METH:
            for d in sorted({int(r['depth_filtered_reads']) for r in rows if r['method'] == m}):
                y = np.array([float(r[ycol]) for r in rows if r['method'] == m and int(r['depth_filtered_reads']) == d])
                k = len(y)
                if k < 3:
                    continue
                ci = _st.t.ppf(0.975, k - 1) * y.std(ddof=1) / np.sqrt(k)
                f.write(
                    f'{m}\t{d}\t{k}\t{y.mean():.1f}\t{y.std(ddof=1):.1f}\t{ci:.1f}\t{y.mean() - ci:.1f}\t{y.mean() + ci:.1f}\n'
                )
# p3 / p3b values are written by composition_panels.py (rule-based composition)
print(SET, 'values written:', sorted(os.listdir(V)))

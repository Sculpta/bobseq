#!/usr/bin/env python3
"""Stacked per-sample figure of the benchmark, native set, same layout as stacked_supplement_24plex.py: one row per metric, the 50 benchmark
samples on a shared x axis in three method blocks (DRUG-seq 24, prime-seq 8, BOBseq 18; block width = number of samples), dots on stems in the method colour.
Within a block samples are ordered by filtered reads as sequenced (ascending), so the same x position is the same sample in every row. Matched-depth rows use
250k filtered reads (the depth of the duplicate-rate panel); samples below it stay empty. Sources: main_figure_panels/values/*.tsv of the native set (the values
behind p0a, p0c/p0d2, p1, p2, p3, p6b, p7e). Output: preprint_figures_native/supplement_stacked/stacked_benchmark_native_per_sample_qc.(svg|png) + values/."""
import os, re, csv, json, ast, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); import sys; sys.path.insert(0, H); from palette import COL, METH
from settings import WORK, COVERAGE, RUN_JSON
P = WORK; V = f'{P}/preprint_figures_native/main_figure_panels/values'; OUT = f'{P}/preprint_figures_native/supplement_stacked'; os.makedirs(f'{OUT}/values', exist_ok=True); D = 250000
RUN = json.load(open(RUN_JSON)); CODE2NAME = {w: l.replace(' ', '_').replace('/', '_') for w, l in RUN['bobcode_labels'].items()}
LAB2NAME = {l.split('\t')[0]: os.path.basename(l.split('\t')[1])[:-4] for l in open(f'{COVERAGE}/samples.tsv') if l.strip()}
def T(f): return list(csv.DictReader(open(f'{V}/{f}'), delimiter='\t'))
names = {m: sorted({r['sample'] for r in T('p6b_picard_balance_native.tsv') if r['method'] == m}) for m in METH}; assert [len(names[m]) for m in METH] == [24, 8, 18], {m: len(v) for m, v in names.items()}
def key(m, s):
    """any table's sample id -> the sample name of the per-sample BAM"""
    if s in names[m]: return s
    if s in LAB2NAME: return LAB2NAME[s]
    if m == 'BOBseq' and s in CODE2NAME: return CODE2NAME[s]
    c = [n for n in names[m] if n.endswith('_' + s)]; assert len(c) == 1, (m, s); return c[0]
R = {m: {n: {} for n in names[m]} for m in METH}
for r in T('p0a_read_fate_native_per_sample.tsv'): v = ast.literal_eval(r['value']); R[r['method']][key(r['method'], r['sample'])].update(map_pct=v[1], uniq_pct=v[2], filt_pct=v[3])
for r in T('p0c_p0d_dedup_per_sample_native.tsv'): R[r['method']][key(r['method'], r['sample'])].update(filtered=float(r['filtered_reads']), mol=float(r['molecules_UCI']), dup250=float(r['duplicate_rate_pct_at_250k']) if r['duplicate_rate_pct_at_250k'] else np.nan)
for f, k, c in (('p1_molecules_vs_depth_native.tsv', 'mol250', 'molecules_UCI'), ('p2_genes_vs_depth_native.tsv', 'genes250', 'genes')):
    for r in T(f):
        if int(r['depth_filtered_reads']) == D: R[r['method']][key(r['method'], r['sample'])][k] = float(r[c])
CL = ['rRNA', 'mitochondrial', 'ribosomal-protein', 'mRNA', 'exonic other biotype', 'intronic', 'intergenic']
for r in T('p3_composition_rules_native_per_sample.tsv'): t = sum(float(r[c]) for c in CL); R[r['method']][key(r['method'], r['sample'])].update({'pct_' + c: 100 * float(r[c]) / t for c in CL})
for r in T('p6b_picard_balance_native.tsv'): R[r['method']][key(r['method'], r['sample'])]['balance'] = float(r['balance_2x_centroid'])
for r in T('p7e_p7f_aligned_bases_per_sample_native.tsv'): R[r['method']][key(r['method'], r['sample'])].update(alen=float(r['mean_aligned_length']), bases=float(r['total_aligned_bases']) / 1e9)
for m in METH:
    for n, d in R[m].items(): assert all(k in d for k in ('filt_pct', 'filtered', 'pct_mRNA', 'balance', 'bases')), (m, n, sorted(d)); d['reads_in'] = d['filtered'] * 100 / d['filt_pct']
ORDER = {m: sorted(names[m], key=lambda n: R[m][n]['filtered']) for m in METH}
METRICS = [('reads into STAR', 'reads_in', True, (1e5, 3e7), None), ('mapping rate\n(% of reads)', 'map_pct', False, (0, 100), None), ('uniquely mapped\n(% of reads)', 'uniq_pct', False, (0, 100), None), ('filtered reads\n(% of reads)', 'filt_pct', False, (0, 100), None),
 ('filtered reads after\nde-duplication (UMI)', 'mol', True, (3e4, 3e6), None), (f'duplicate rate (%)\nat {D//1000}k filtered reads', 'dup250', False, (0, 100), None), (f'molecules (UMI)\nat {D//1000}k filtered reads', 'mol250', False, (0, 260000), None),
 (f'genes at {D//1000}k filtered reads\n(no ribosomal-protein genes)', 'genes250', False, (0, 13000), None), ('rRNA content\n(% of mapped)', 'pct_rRNA', False, (0, 100), None), ('mRNA fraction\n(% of mapped)', 'pct_mRNA', False, (0, 100), None),
 ('intronic reads\n(% of mapped)', 'pct_intronic', False, (0, 60), None), ('mitochondrial reads\n(% of mapped)', 'pct_mitochondrial', False, (0, 20), None), ('mean aligned length\nper read (nt)', 'alen', False, (0, 160), None), ('aligned bases (Gb)', 'bases', True, (0.003, 3), None), ("5'-3' balance\n(2 x centroid)", 'balance', False, (0, 2), 1.0)]
plt.rcParams.update({'svg.fonttype': 'none', 'font.family': ['Arial', 'Nimbus Sans', 'DejaVu Sans'], 'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6, 'ytick.labelsize': 6, 'axes.linewidth': 0.6, 'axes.spines.top': False, 'axes.spines.right': False, 'text.color': '#231F20', 'axes.edgecolor': '#231F20', 'axes.labelcolor': '#231F20', 'xtick.color': '#231F20', 'ytick.color': '#231F20'})
fig, axs = plt.subplots(len(METRICS), len(METH), figsize=(7.2, 0.95 * len(METRICS) + 0.6), sharey='row', sharex='col', gridspec_kw={'wspace': 0.06, 'hspace': 0.28, 'width_ratios': [len(ORDER[m]) + 1 for m in METH]}); val = []
for i, (lab, k, log, ylim, ref) in enumerate(METRICS):
    for j, m in enumerate(METH):
        ax = axs[i, j]; y = np.array([R[m][n].get(k, np.nan) for n in ORDER[m]]); x = np.arange(1, len(y) + 1); ok = ~np.isnan(y)
        assert (y[ok] >= ylim[0]).all() and (y[ok] <= ylim[1]).all(), (lab, m, float(np.nanmin(y)), float(np.nanmax(y)))
        ax.vlines(x[ok], ylim[0], y[ok], color='#BDBDBD', lw=0.6, zorder=1); ax.scatter(x[ok], y[ok], s=8, color=COL[m], lw=0, zorder=3); [val.append((lab.replace('\n', ' '), m, n, xi, v)) for n, xi, v in zip(ORDER[m], x, y)]
        if log: ax.set_yscale('log')
        ax.set_ylim(*ylim); ax.set_xlim(0.2, len(y) + 0.8); ax.set_xticks([t for t in (1, 5, 10, 15, 20, 24) if t <= len(y)] if len(y) > 8 else [1, 4, 8])
        if ref is not None: ax.axhline(ref, color='#9a9a9a', lw=0.5, ls=(0, (3, 2)), zorder=0)
        if j > 0: ax.spines['left'].set_visible(False); ax.tick_params(axis='y', left=False, which='both')
        if i < len(METRICS) - 1: ax.tick_params(axis='x', bottom=False)
        if i == 0: ax.set_title(f'{m} (n={len(y)})', fontsize=7, pad=4, color=COL[m], fontweight='bold')
    axs[i, 0].set_ylabel(lab, fontsize=6.5, rotation=0, ha='right', va='center', labelpad=6)
fig.supxlabel('benchmark samples, ordered by filtered reads within each method (same order in every row)', fontsize=7, y=0.085)
name = 'stacked_benchmark_native_per_sample_qc'
for ext, kw in (('svg', {}), ('png', {'dpi': 200})): fig.savefig(f'{OUT}/{name}.{ext}', bbox_inches='tight', pad_inches=0.05, **kw)
s = open(f'{OUT}/{name}.svg').read(); open(f'{OUT}/{name}.svg', 'w').write(re.sub(r"font-family:\s*[^;\"]+", 'font-family: Arial', s))
with open(f'{OUT}/values/{name}.tsv', 'w') as o: o.write('metric\tmethod\tsample\tx_position\tvalue\n'); [o.write(f'{a}\t{b}\t{c}\t{d}\t' + ('' if np.isnan(e) else f'{e:.6g}') + '\n') for a, b, c, d, e in val]
print('STACKED BENCHMARK DONE', len(val), 'cells;', {m: int(sum(1 for n in ORDER[m] if not np.isnan(R[m][n].get('dup250', np.nan)))) for m in METH}, f'samples at {D//1000}k')

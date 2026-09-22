#!/usr/bin/env python3
"""Stacked per-sample QC figure of the 24-plex run in the layout of the species-mixing supplement: one row per metric, the 24 samples on a
shared x axis in eight condition blocks of three replicates (mouse last), dots on stems. Only metrics that are one number per sample ("stacks directly").
Sources: supplement_per_sample_21/values/per_well_native.tsv (21 human wells), supplement_per_sample_mouse/values/per_well_native_mouse.tsv (3 mouse wells),
results/barcode_accuracy_funnel/funnel_summary.json (mouse accuracy). Output: preprint_figures_native/supplement_stacked/ (SVG + PNG, values/ with the plotted values); it replaces the single-metric per-sample panels it covers (kept for reference)."""
import os, re, csv, json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from settings import WORK, PROJECT
P = WORK; PR = PROJECT; OUT = f'{P}/preprint_figures_native/supplement_stacked'; os.makedirs(f'{OUT}/values', exist_ok=True)
hum = list(csv.DictReader(open(f'{P}/preprint_figures_native/supplement_per_sample_21/values/per_well_native.tsv'), delimiter='\t')); mus = list(csv.DictReader(open(f'{P}/preprint_figures_native/supplement_per_sample_mouse/values/per_well_native_mouse.tsv'), delimiter='\t'))
fs = json.load(open(f'{PR}/results/barcode_accuracy_funnel/funnel_summary.json'))['per_sample']
for r in mus: r['group'] = 'RAW 264.7'; r['barcode_accuracy_pct'] = f"{fs[r['sample']]['levels'][-1]['accuracy']:.3f}"
BLOCKS = [('Hek control', 'HEK293T\ncontrol'), ('CADM1 1nM', 'CADM1 ASO\n1 nM'), ('CADM1 100nM', 'CADM1 ASO\n100 nM'), ('CHX dose 1ug mL', 'CHX\n1 ug/mL'), ('CHX dose 50ug mL', 'CHX\n50 ug/mL'), ('Ris dose 25mM', 'risdiplam\n25'), ('Ris dose 500mM', 'risdiplam 500\n(input failure)'), ('RAW 264.7', 'RAW 264.7\n(mouse)')]
rows = hum + mus; by = {g: sorted([r for r in rows if r['group'] == g], key=lambda r: r['sample']) for g, _ in BLOCKS}; assert all(len(by[g]) == 3 for g, _ in BLOCKS), {g: len(v) for g, v in by.items()}
f = lambda r, k: float(r[k]) if r.get(k, '') not in ('', None) else np.nan
METRICS = [  # (label, value function(s), log, ylim, reference line)
 ('reads into STAR', [lambda r: f(r, 'reads_in')], True, (5e4, 2e7), None),
 ('mapping rate\n(% of reads)', [lambda r: 100 * f(r, 'mapped') / f(r, 'reads_in')], False, (0, 100), None),
 ('uniquely mapped\n(% of reads)', [lambda r: 100 * f(r, 'unique') / f(r, 'reads_in')], False, (0, 100), None),
 ('filtered reads\n(% of reads)', [lambda r: 100 * f(r, 'filtered_reads') / f(r, 'reads_in')], False, (0, 100), None),
 ('duplicate rate (%)', [lambda r: f(r, 'duplicate_rate_pct')], False, (0, 100), None),
 ('filtered reads after\nde-duplication (UMI)', [lambda r: f(r, 'molecules_UCI')], True, (5e4, 3e6), None),
 ('rRNA content\n(% of mapped)', [lambda r: f(r, 'pct_rRNA')], False, (0, 100), None),
 ('mRNA fraction\n(% of mapped)', [lambda r: f(r, 'pct_mRNA')], False, (0, 100), None),
 ('intronic reads\n(% of mapped)', [lambda r: f(r, 'pct_intronic')], False, (0, 20), None),
 ('median aligned length (nt)\nfilled read 2, open read 1', [lambda r: f(r, 'median_aligned_len_mate1'), lambda r: f(r, 'median_aligned_len_mate2')], False, (0, 160), None),
 ('insert length\nmedian (nt)', [lambda r: f(r, 'insert_median')], False, (0, 300), None),
 ('aligned bases (Gb)', [lambda r: f(r, 'total_aligned_bases') / 1e9], True, (0.01, 3), None),
 ("5'-3' balance\n(2 x centroid)", [lambda r: f(r, 'picard_balance')], False, (0, 2), 1.0)]   # no barcode accuracy row: shown in its own panel
plt.rcParams.update({'svg.fonttype': 'none', 'font.family': ['Arial', 'Nimbus Sans', 'DejaVu Sans'], 'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6, 'ytick.labelsize': 6, 'axes.linewidth': 0.6, 'axes.spines.top': False, 'axes.spines.right': False, 'text.color': '#231F20', 'axes.edgecolor': '#231F20', 'axes.labelcolor': '#231F20', 'xtick.color': '#231F20', 'ytick.color': '#231F20'})
fig, axs = plt.subplots(len(METRICS), len(BLOCKS), figsize=(7.2, 0.95 * len(METRICS) + 0.6), sharey='row', sharex=True, gridspec_kw={'wspace': 0.12, 'hspace': 0.28}); DOT = '#3A3A3A'; STEM = '#BDBDBD'
val = []
for i, (lab, fns, log, ylim, ref) in enumerate(METRICS):
    for j, (g, title) in enumerate(BLOCKS):
        ax = axs[i, j]; base = ylim[0]
        for x, r in enumerate(by[g]):
            for k, fn in enumerate(fns):
                v = fn(r)
                if np.isnan(v): continue
                xx = x + (0 if len(fns) == 1 else (-0.17 if k == 0 else 0.17)); ax.vlines(xx, base, v, color=STEM, lw=0.6, zorder=1)
                ax.scatter([xx], [v], s=9, zorder=3, **({'color': DOT, 'lw': 0} if k == 0 else {'facecolors': 'white', 'edgecolors': DOT, 'lw': 0.7})); val.append((lab.replace('\n', ' '), g, r['sample'], k + 1, v))
        if log: ax.set_yscale('log')
        ax.set_ylim(*ylim); ax.set_xlim(-0.6, 2.6); ax.set_xticks([0, 1, 2]); ax.set_xticklabels(['1', '2', '3'])
        if ref is not None: ax.axhline(ref, color='#9a9a9a', lw=0.5, ls=(0, (3, 2)), zorder=0)
        if j > 0: ax.spines['left'].set_visible(False); ax.tick_params(axis='y', left=False, which='both')
        if i < len(METRICS) - 1: ax.tick_params(axis='x', bottom=False)
        if i == 0: ax.set_title(title, fontsize=6.5, pad=4)
    axs[i, 0].set_ylabel(lab, fontsize=6.5, rotation=0, ha='right', va='center', labelpad=6)
fig.supxlabel('replicate, within treatment (BOBseq 24-plex, 21 human wells and 3 mouse wells)', fontsize=7, y=0.085)
for ext, kw in (('svg', {}), ('png', {'dpi': 200})): fig.savefig(f'{OUT}/stacked_24plex_per_sample_qc.{ext}', bbox_inches='tight', pad_inches=0.05, **kw)
s = open(f'{OUT}/stacked_24plex_per_sample_qc.svg').read(); open(f'{OUT}/stacked_24plex_per_sample_qc.svg', 'w').write(re.sub(r"font-family:\s*[^;\"]+", 'font-family: Arial', s))
with open(f'{OUT}/values/stacked_24plex_per_sample_qc.tsv', 'w') as o: o.write('metric\ttreatment\tsample\tseries\tvalue\n'); [o.write(f'{a}\t{b}\t{c}\t{d}\t{e:.6g}\n') for a, b, c, d, e in val]
print('STACKED DONE', len(val), 'values')

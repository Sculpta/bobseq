#!/usr/bin/env python3
"""One stacked per-sample QC figure for the supplement: the benchmark methods and all 24 BOBseq samples in one table of 14 metrics (molecules at 250k is not a row: it is the same information as the duplicate rate at 250k).
Columns: DRUG-seq (24 DMSO wells), prime-seq (8 HEK293T samples), then the BOBseq 24-plex by treatment (6 benchmark treatments x 3 replicates, the risdiplam 500
input failures, the 3 mouse wells in their own colour, set apart). Replicates of a column are stacked: one dot per sample at the column position, short horizontal
line = mean. Matched-depth rows use 250k filtered reads; a sample with fewer filtered reads is not subsampled and is marked with a small grey cross on the baseline.
Sources: main_figure_panels/values (competitors: values behind p0a, p0c/p0d2, p1, p2, p3, p6b, p7e), supplement_per_sample_21 and supplement_per_sample_mouse per-well
tables (BOBseq), rarefied.tsv of the BOBseq arm (stats_E_U, stats_E_U_supp21, stats_E_U_mouse3). Output: preprint_figures_native/supplement_stacked/stacked_qc_all_samples.*
"""

from settings import WORK, COVERAGE
import os, re, csv, ast, json, numpy as np, matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

H = os.path.dirname(os.path.abspath(__file__))
import sys

sys.path.insert(0, H)
from palette import COL

P = WORK
FN = f'{P}/preprint_figures_native'
V = f'{FN}/main_figure_panels/values'
ARM = f'{P}/benchmark_uniform/bob57_24plex_pe_native'
OUT = f'{FN}/supplement_stacked'
os.makedirs(f'{OUT}/values', exist_ok=True)
D = 250000


def T(p):
    return list(csv.DictReader(open(p), delimiter='\t'))


f = lambda r, k: float(r[k]) if r.get(k, '') not in ('', None) else np.nan
# ---- competitors (keyed by the per-sample BAM name) ----
names = {
    m: sorted({r['sample'] for r in T(f'{V}/p6b_picard_balance_native.tsv') if r['method'] == m})
    for m in ('DRUG-seq', 'prime-seq')
}
LAB2NAME = {
    l.split('\t')[0]: os.path.basename(l.split('\t')[1])[:-4] for l in open(f'{COVERAGE}/samples.tsv') if l.strip()
}


def key(m, s):
    if s in names[m]:
        return s
    if s in LAB2NAME:
        return LAB2NAME[s]
    c = [n for n in names[m] if n.endswith('_' + s)]
    assert len(c) == 1, (m, s)
    return c[0]


R = {m: {n: {} for n in names[m]} for m in names}
CL = ['rRNA', 'mitochondrial', 'ribosomal-protein', 'mRNA', 'exonic other biotype', 'intronic', 'intergenic']
for r in T(f'{V}/p0a_read_fate_native_per_sample.tsv'):
    if r['method'] in R:
        v = ast.literal_eval(r['value'])
        R[r['method']][key(r['method'], r['sample'])].update(map_pct=v[1], uniq_pct=v[2], filt_pct=v[3])
for r in T(f'{V}/p0c_p0d_dedup_per_sample_native.tsv'):
    if r['method'] in R:
        R[r['method']][key(r['method'], r['sample'])].update(
            filtered=f(r, 'filtered_reads'), mol=f(r, 'molecules_UCI'), dup250=f(r, 'duplicate_rate_pct_at_250k')
        )
for fn, k, c in (
    ('p1_molecules_vs_depth_native.tsv', 'mol250', 'molecules_UCI'),
    ('p2_genes_vs_depth_native.tsv', 'genes250', 'genes'),
):
    for r in T(f'{V}/{fn}'):
        if r['method'] in R and int(r['depth_filtered_reads']) == D:
            R[r['method']][key(r['method'], r['sample'])][k] = float(r[c])
for r in T(f'{V}/p3_composition_rules_native_per_sample.tsv'):
    if r['method'] in R:
        t = sum(float(r[c]) for c in CL)
        R[r['method']][key(r['method'], r['sample'])].update({'pct_' + c: 100 * float(r[c]) / t for c in CL})
for r in T(f'{V}/p6b_picard_balance_native.tsv'):
    if r['method'] in R:
        R[r['method']][key(r['method'], r['sample'])]['balance'] = float(r['balance_2x_centroid'])
for r in T(f'{V}/p7e_p7f_aligned_bases_per_sample_native.tsv'):
    if r['method'] in R:
        R[r['method']][key(r['method'], r['sample'])]['bases'] = float(r['total_aligned_bases']) / 1e9
for m in R:
    for n, d in R[m].items():
        d['reads_in'] = d['filtered'] * 100 / d['filt_pct']
# ---- BOBseq, all 24 wells ----
rar = {}
for s_ in ('stats_E_U', 'stats_E_U_supp21', 'stats_E_U_mouse3'):
    for r in T(f'{ARM}/{s_}/rarefied.tsv'):
        if int(r['depth']) == D and r['at_native'] == '0':
            rar[r['well']] = r
B = {}
for r in T(f'{FN}/supplement_per_sample_21/values/per_well_native.tsv') + T(
    f'{FN}/supplement_per_sample_mouse/values/per_well_native_mouse.tsv'
):
    w = r['bobcode']
    q = rar.get(w)
    g = r.get('group') or 'RAW 264.7'
    B[r['sample']] = dict(
        group=g,
        reads_in=f(r, 'reads_in'),
        map_pct=100 * f(r, 'mapped') / f(r, 'reads_in'),
        uniq_pct=100 * f(r, 'unique') / f(r, 'reads_in'),
        filt_pct=100 * f(r, 'filtered_reads') / f(r, 'reads_in'),
        filtered=f(r, 'filtered_reads'),
        mol=f(r, 'molecules_UCI'),
        dup250=100 * (1 - float(q['B_corr']) / D) if q else np.nan,
        mol250=float(q['B_corr']) if q else np.nan,
        genes250=float(q['genes_B_noRP']) if q else np.nan,
        insert=f(r, 'insert_median'),
        bases=f(r, 'total_aligned_bases') / 1e9,
        balance=f(r, 'picard_balance'),
        **{'pct_' + c: f(r, 'pct_' + c) for c in CL},
    )
assert len(B) == 24
MOUSE = '#7A4FA3'
FAIL = '#8C8C8C'
BOB = COL['BOBseq']
COLS = [
    ('DRUG-seq', 'DRUG-seq', COL['DRUG-seq'], list(R['DRUG-seq'].values())),
    ('prime-seq', 'prime-seq', COL['prime-seq'], list(R['prime-seq'].values())),
]
for g, lab, c in (
    ('Hek control', 'HEK control', BOB),
    ('CADM1 1nM', 'CADM1 1 nM', BOB),
    ('CADM1 100nM', 'CADM1 100 nM', BOB),
    ('CHX dose 1ug mL', 'CHX 1 ug/mL', BOB),
    ('CHX dose 50ug mL', 'CHX 50 ug/mL', BOB),
    ('Ris dose 25mM', 'Ris 25 nM', BOB),
    ('Ris dose 500mM', 'Ris 500 nM', FAIL),
    ('RAW 264.7', 'RAW control', MOUSE),
):
    v = [d for n, d in sorted(B.items()) if d['group'] == g]
    assert len(v) == 3, (g, len(v))
    COLS.append((g, lab, c, v))
XPOS = [0, 1] + [2.35 + i for i in range(7)] + [9.7]  # gaps: competitors | BOBseq human | mouse
METRICS = [
    ('reads into STAR (log)', 'reads_in', True, (1e5, 3e7)),
    ('mapping rate\n(% of reads)', 'map_pct', False, (0, 100)),
    ('uniquely mapped\n(% of reads)', 'uniq_pct', False, (0, 100)),
    ('filtered reads\n(% of reads)', 'filt_pct', False, (0, 100)),
    ('filtered reads after\nde-duplication (UMI, log)', 'mol', True, (3e4, 3e6)),
    (f'duplicate rate (%)\nat {D//1000}k filtered reads', 'dup250', False, (0, 100)),
    (f'genes at {D//1000}k filtered reads\n(no ribosomal-protein genes)', 'genes250', False, (0, 13000)),
    ('rRNA content\n(% of mapped)', 'pct_rRNA', False, (0, 100)),
    ('mRNA fraction\n(% of mapped)', 'pct_mRNA', False, (0, 100)),
    ('intronic reads\n(% of mapped)', 'pct_intronic', False, (0, 30)),
    ('mitochondrial reads\n(% of mapped)', 'pct_mitochondrial', False, (0, 20)),
    ('insert length\nmedian (nt)', 'insert', False, (0, 300)),
    ('aligned bases (Gb, log)', 'bases', True, (0.003, 3)),
    ("5'-3' balance\n(2 x centroid)", 'balance', False, (0, 2)),
]
plt.rcParams.update(
    {
        'svg.fonttype': 'none',
        'font.family': ['Arial', 'Nimbus Sans', 'DejaVu Sans'],
        'font.size': 7,
        'axes.linewidth': 0.6,
        'ytick.labelsize': 6,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'text.color': '#231F20',
        'axes.edgecolor': '#231F20',
        'axes.labelcolor': '#231F20',
        'xtick.color': '#231F20',
        'ytick.color': '#231F20',
    }
)
fig, axs = plt.subplots(
    len(METRICS), 1, figsize=(4.7, 0.7 * len(METRICS) + 1.7), sharex=True, gridspec_kw={'hspace': 0.3}
)
val = []


def stagger(y, log, ylim):
    """narrow placement. Triplicates: replicate 1, 2, 3 from left to right in every row (fixed, so no order is implied by the values).
    Larger columns: a compact swarm. Points are placed from the median outwards; each takes the offset closest to the column centre that does not overlap a placed point
    (distance measured in the row's axis units, log rows in log space), so the cluster is symmetric and as narrow as the data allow.
    """
    k = len(y)
    if k <= 3:
        return np.array([-0.2, 0.0, 0.2])[:k]
    t = np.log10(y) if log else np.asarray(y, float)
    lo, hi = (np.log10(ylim[0]), np.log10(ylim[1])) if log else ylim
    yn = (t - lo) / (hi - lo)
    off = np.zeros(k)
    placed = []
    cand = [0.0] + [sgn * st * 0.06 for st in range(1, 7) for sgn in (-1, 1)]
    med = np.nanmedian(yn)
    for idx in np.argsort(np.where(np.isnan(yn), np.inf, np.abs(yn - med)), kind='stable'):
        if np.isnan(yn[idx]):
            continue
        for c in cand:
            if all(((c - px) / 0.06) ** 2 + ((yn[idx] - py) / 0.06) ** 2 >= 1 for px, py in placed):
                off[idx] = c
                break
        else:
            off[idx] = max(
                cand, key=lambda c: min(((c - px) / 0.06) ** 2 + ((yn[idx] - py) / 0.06) ** 2 for px, py in placed)
            )  # column full at this height: the least crowded position
        placed.append((off[idx], yn[idx]))
    nan = np.flatnonzero(np.isnan(yn))
    off[nan] = (np.arange(len(nan)) - (len(nan) - 1) / 2) * 0.12  # samples without a value (grey crosses) side by side
    return off


for ri, (lab, k, log, ylim) in enumerate(METRICS):
    ax = axs[ri]
    for ci, (g, clab, col, samples) in enumerate(COLS):
        y = np.array([d.get(k, np.nan) for d in samples], float)
        x = XPOS[ci] + stagger(y, log, ylim)
        ok = ~np.isnan(y)
        if k == 'insert' and not ok.any():
            ax.text(
                XPOS[ci],
                np.mean(ylim),
                'single-end',
                ha='center',
                va='center',
                fontsize=5.5,
                color='#8a8a8a',
                rotation=90,
            )
            continue
        assert ((y[ok] >= ylim[0]) & (y[ok] <= ylim[1])).all(), (lab, g, float(np.nanmin(y)), float(np.nanmax(y)))
        ax.scatter(x[ok], y[ok], s=4.5 if len(y) > 3 else 8, color=col, lw=0, alpha=0.85, zorder=3)
        if ok.any():
            mu = float(np.exp(np.log(y[ok]).mean())) if log else float(y[ok].mean())
            ax.hlines(mu, XPOS[ci] - 0.38, XPOS[ci] + 0.38, color='#231F20', lw=0.8, zorder=4)
        if (~ok).any() and k in ('dup250', 'mol250', 'genes250'):
            ax.scatter(
                x[~ok], np.full((~ok).sum(), ylim[0]), marker='x', s=9, color='#9a9a9a', lw=0.6, zorder=5, clip_on=False
            )
        [val.append((lab.replace('\n', ' '), g, i + 1, '' if np.isnan(v_) else f'{v_:.6g}')) for i, v_ in enumerate(y)]
    if log:
        ax.set_yscale('log')
    ax.set_ylim(*ylim)
    ax.set_xlim(-0.7, XPOS[-1] + 0.7)
    ax.set_ylabel(lab, fontsize=6.3, rotation=0, ha='right', va='center', labelpad=6)
    ax.tick_params(axis='x', bottom=False)
    if k == 'balance':
        ax.axhline(1, color='#9a9a9a', lw=0.5, ls=(0, (3, 2)), zorder=0)
    for xs in (1.68, 9.03):
        ax.axvline(xs, color='#D0D0D0', lw=0.5, zorder=0)
axs[-1].set_xticks(XPOS)
axs[-1].set_xticklabels(
    [c[1].replace('\n', ' ') + (f' (n={len(c[3])})' if len(c[3]) > 3 else '') for c in COLS], fontsize=6.5, rotation=90
)
axs[-1].tick_params(axis='x', bottom=True, length=2)
for t, c in zip(axs[-1].get_xticklabels(), COLS):
    t.set_color(c[2])
top = axs[0]
tr = top.get_xaxis_transform()
for x0, x1, txt, c in (
    (-0.42, 1.42, 'public data', '#555555'),
    (1.93, 8.77, 'BOBseq 24-plex, human', BOB),
    (9.28, 10.12, 'mouse', MOUSE),
):
    top.plot([x0, x1], [1.18, 1.18], color=c, lw=1.4, transform=tr, clip_on=False, solid_capstyle='butt')
    top.text(
        (x0 + x1) / 2,
        1.26,
        txt,
        color=c,
        fontsize=7,
        fontweight='bold',
        ha='center',
        va='bottom',
        transform=tr,
        linespacing=0.95,
    )
fig.text(
    0.5,
    0.0,
    f'one dot per sample, line = mean (geometric mean on log rows); three replicates per BOBseq column;\ngrey cross = sample with fewer than {D//1000}k filtered reads, not subsampled; Ris 500 nM = input failures;\nmouse rRNA includes the reads the aligner places on the human rDNA copies',
    ha='center',
    va='top',
    fontsize=5.5,
    color='#555',
)
name = 'stacked_qc_all_samples'
for ext, kw in (('svg', {}), ('png', {'dpi': 200})):
    fig.savefig(f'{OUT}/{name}.{ext}', bbox_inches='tight', pad_inches=0.05, **kw)
s = open(f'{OUT}/{name}.svg').read()
open(f'{OUT}/{name}.svg', 'w').write(re.sub(r"font-family:\s*[^;\"]+", 'font-family: Arial', s))
with open(f'{OUT}/values/{name}.tsv', 'w') as o:
    o.write('metric\tcolumn\tsample_index_in_column\tvalue\n')
    [o.write('\t'.join(map(str, r)) + '\n') for r in val]
print('MERGED STACK DONE', len(val), 'cells;', {c[0]: len(c[3]) for c in COLS})

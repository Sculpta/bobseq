#!/usr/bin/env python3
"""Threshold versions of the all-reads heatmap and the mean poly(A)-anchored profile, with gene cutoffs:
transcripts >= 1 kb, ribosomal-protein genes (RPL/RPS/MRPL/MRPS/RPLP; RPS6K kinases spared) excluded, no fragment-length rule.
  p4i2_coverage_heatmaps_all_reads_1kb_noRP      one row per gene, ALL unique reads, 50 bins, each row scaled to its own peak, short to long
  p4m2_coverage_vs_nt_from_polyA_mean_1kb_noRP   coverage / gene mean at 200 reads per gene (seed 11), last 3,000 nt anchored at the poly(A)
                                                 site, mean over genes, 25-nt running mean; shorter genes contribute where they exist
  p4_thresholds_numbers.json                     gene set and profile numbers
Gene set: >= 200 unique reads in every method (pooled samples of the set), canonical transcript >= 1 kb, not a ribosomal-protein gene.
BM_COVSET=native|50nt selects the set (positions_canonical{,_50nt}, samples{,_50nt}.tsv)."""

from settings import WORK, COVERAGE
import os, re, sys, json, numpy as np, matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

P = WORK
A = COVERAGE
SET = os.environ.get('BM_COVSET', 'native')
SFX = '' if SET == 'native' else f'_{SET}'
FIG = 'preprint_figures_native' if SET == 'native' else f'preprint_figures_{SET}'
OUT = f'{P}/{FIG}/main_figure_panels'
os.makedirs(OUT, exist_ok=True)
from palette import COL, METH

RIBO = re.compile(r'^(RP[LS]|RPLP|MRP[LS])[^K]*$', re.I)
S = [l.rstrip('\n').split('\t') for l in open(f'{A}/samples{SFX}.tsv') if l.strip()]
POOL = {}
for lab, bam, m in S:
    z = np.load(f'{A}/positions_canonical{SFX}/{lab}.npz')
    for k in ('gene', 't', 'alen'):
        POOL.setdefault(m, {}).setdefault(k, []).append(z[k])
    tx_len = z['tx_len']
    genes = list(z['genes'])
for m in METH:
    POOL[m] = {k: np.concatenate(v) for k, v in POOL[m].items()}


def by_gene(d):
    g = d['gene']
    o = np.argsort(g, kind='stable')
    gs = g[o]
    st = np.flatnonzero(np.r_[True, gs[1:] != gs[:-1]])
    en = np.r_[st[1:], len(gs)]
    return {int(gs[s]): o[s:e] for s, e in zip(st, en)}


IDX = {m: by_gene(POOL[m]) for m in METH}


def gene_cov(d, ii, L):
    t = d['t'][ii].astype(int)
    e = np.minimum(t + d['alen'][ii].astype(int), L)
    diff = np.zeros(L + 1)
    np.add.at(diff, t, 1)
    np.add.at(diff, e, -1)
    return np.cumsum(diff)[:L]


def sm(y, k=25):
    y = np.asarray(y, float)
    v = np.where(np.isnan(y), 0, y)
    n = (~np.isnan(y)).astype(float)
    ker = np.ones(k)
    return np.convolve(v, ker, 'same') / np.maximum(np.convolve(n, ker, 'same'), 1)


plt.rcParams.update(
    {
        'svg.fonttype': 'none',
        'font.family': ['Arial', 'Nimbus Sans', 'DejaVu Sans'],
        'font.size': 7,
        'axes.labelsize': 7,
        'xtick.labelsize': 6.5,
        'ytick.labelsize': 6.5,
        'axes.linewidth': 0.6,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 2.5,
        'ytick.major.size': 2.5,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'axes.spines.top': False,
        'axes.spines.right': False,
        'text.color': '#231F20',
        'axes.edgecolor': '#231F20',
        'axes.labelcolor': '#231F20',
        'xtick.color': '#231F20',
        'ytick.color': '#231F20',
        'axes.grid': False,
    }
)


def psave(fig, name):
    fig.tight_layout(pad=0.4)
    fig.savefig(f'{OUT}/{name}.svg', bbox_inches='tight', pad_inches=0.03)
    fig.savefig(f'{OUT}/{name}.png', dpi=200, bbox_inches='tight', pad_inches=0.03)
    plt.close(fig)
    s = open(f'{OUT}/{name}.svg').read()
    s = re.sub(r"font-family:\s*[^;\"]+", 'font-family: Arial', s)
    open(f'{OUT}/{name}.svg', 'w').write(s)
    print('wrote', name)


D0 = 200
cand = [g for g in IDX['DRUG-seq'] if all(len(IDX[m].get(g, [])) >= D0 for m in METH) and tx_len[g] >= 1000]
common = sorted([g for g in cand if not RIBO.match(genes[g])], key=lambda g: tx_len[g])
rp_removed = sorted(genes[g] for g in cand if RIBO.match(genes[g]))
print(
    SET,
    len(common),
    'genes (>= 200 unique reads in every method, >= 1 kb, no RP);',
    len(rp_removed),
    'RP genes removed',
)
NB = 50
H = {m: np.zeros((len(common), NB)) for m in METH}
for k, g in enumerate(common):
    L = int(tx_len[g])
    idx = np.minimum((np.arange(L) * NB) // L, NB - 1)
    for m in METH:
        c = gene_cov(POOL[m], IDX[m][g], L)
        b = np.bincount(idx, weights=c, minlength=NB) / np.bincount(idx, minlength=NB)
        H[m][k] = b / b.max() if b.max() > 0 else b
fig, axs = plt.subplots(
    1, len(METH), figsize=(1.3 * len(METH), 2.4), sharey=True
)  # one panel per method
for ax, m in zip(axs, METH):
    ax.imshow(H[m], aspect='auto', cmap='Greys', vmin=0, vmax=1, interpolation='nearest')
    ax.set_title(m, fontsize=7, color=COL[m])
    ax.set_xticks([0, NB - 1])
    ax.set_xticklabels(["5'", "3'"])
    ax.tick_params(length=2)
axs[0].set_ylabel(
    f'{len(common):,} genes >= 1 kb (short to long),\nribosomal-protein genes excluded,\nall unique reads per gene'
)
axs[0].set_yticks([])
fig.text(0.5, -0.04, 'coverage along the mature transcript (each row scaled to its own peak)', ha='center', fontsize=7)
psave(fig, 'p4i2_coverage_heatmaps_all_reads_1kb_noRP')
W = 3000
rng = np.random.default_rng(11)
U = {m: np.full((len(common), W), np.nan) for m in METH}
for k, g in enumerate(common):
    L = int(tx_len[g])
    w = min(W, L)
    for m in METH:
        ii = rng.choice(IDX[m][g], D0, replace=False)
        c = gene_cov(POOL[m], ii, L)
        z = c / c.mean()
        U[m][k, :w] = z[::-1][:w]
navail = np.sum(~np.isnan(U['BOBseq']), 0)
prof = {m: sm(np.nanmean(U[m], 0)) for m in METH}
fig, ax = plt.subplots(figsize=(2.6, 2.2))
xx = np.arange(W)
for m in METH:
    ax.plot(xx, prof[m], color=COL[m], lw=1.6)
ax.axhline(1, color='#9a9a9a', lw=0.6, ls=(0, (2, 2)))
ax.set_xlim(W, 0)
ax.set_ylim(0, None)
ax.set_xlabel('distance from the poly(A) site (nt)' + ('' if SET == 'native' else '\n(all methods: one 50-nt read)'))
ax.set_ylabel(
    'coverage relative to gene mean\n(mean over genes, transcripts >= 1 kb,\nribosomal-protein genes excluded)'
)
psave(fig, 'p4m2_coverage_vs_nt_from_polyA_mean_1kb_noRP')
VAL = f'{OUT}/values'
os.makedirs(VAL, exist_ok=True)  # raw values for replotting
with open(f'{VAL}/p4i2_coverage_heatmaps_all_reads_1kb_noRP.tsv', 'w') as f:
    f.write('gene\ttranscript_length_nt\tmethod\t' + '\t'.join(f'bin{b+1:02d}' for b in range(NB)) + '\n')
    for k, g in enumerate(common):
        for m in METH:
            f.write(f'{genes[g]}\t{int(tx_len[g])}\t{m}\t' + '\t'.join(f'{v:.4f}' for v in H[m][k]) + '\n')
with open(f'{VAL}/p4m2_coverage_vs_nt_from_polyA_mean_1kb_noRP.tsv', 'w') as f:
    f.write(
        'nt_from_polyA_site\tgenes_contributing\t'
        + '\t'.join(f'{m}_mean_over_genes_smoothed25' for m in METH)
        + '\t'
        + '\t'.join(f'{m}_mean_over_genes_raw' for m in METH)
        + '\n'
    )
    for d in range(W):
        f.write(
            f'{d}\t{int(navail[d])}\t'
            + '\t'.join(f'{prof[m][d]:.4f}' for m in METH)
            + '\t'
            + '\t'.join(f'{np.nanmean(U[m][:, d]):.4f}' for m in METH)
            + '\n'
        )
json.dump(
    {
        'set': SET,
        'genes': len(common),
        'gene_names': [genes[g] for g in common],
        'rp_genes_removed': rp_removed,
        'reads_per_gene_profile': D0,
        'seed': 11,
        'window_nt': W,
        'genes_contributing_at_nt': {str(d): int(navail[d - 1]) for d in (500, 1000, 2000, 3000)},
        'profile': {
            m: {
                'peak': float(np.nanmax(prof[m])),
                'peak_nt_from_polyA': int(np.nanargmax(prof[m])),
                'mean_0_300': float(np.nanmean(prof[m][:300])),
                'mean_1000_2000': float(np.nanmean(prof[m][1000:2000])),
            }
            for m in METH
        },
        'heatmap_frac_bins_above_20pct_of_peak': {m: float(np.mean(H[m] > 0.2)) for m in METH},
    },
    open(f'{OUT}/p4_thresholds_numbers.json', 'w'),
    indent=1,
)
print('DONE')

#!/usr/bin/env python3
"""Unscaled coverage profiles in nucleotides from the poly(A) site (and from the cap): the view chosen from picard_custom for the main figure (p4m).
Genes with >= d0 unique reads in every method, each subsampled to d0 reads, per-base coverage on the canonical transcript (BOBseq both mates), equal gene weight.
Writes the variant sheet Q_unscaled_variants.png (+ numbers) and preprint-style panels p4m*, p4n*, p4o*."""
import os, sys, json, textwrap, re, logging, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt; logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
from matplotlib.lines import Line2D
from settings import WORK, COVERAGE
A = COVERAGE; P = WORK
SET = os.environ.get('BM_COVSET', 'native'); SFX = '' if SET == 'native' else f'_{SET}'; FIG = 'preprint_figures_native' if SET == 'native' else f'preprint_figures_{SET}'
BOBNOTE = 'BOBseq: both mates' if SET == 'native' else 'every method one 50-nt single read per fragment'; READNOTE = ('BOBseq both mates of the 2 x 150 pairs; competitors one native-length read per fragment' if SET == 'native' else 'read-length matched set: every method one 50-nt single read per fragment (BOBseq R2 truncated to 50 nt)')
POS = f'{A}/positions_canonical{SFX}'; OUT = f'{P}/{FIG}/coverage_architecture'; OUTP = f'{P}/{FIG}/main_figure_panels'; os.makedirs(OUT, exist_ok=True); os.makedirs(OUTP, exist_ok=True)
from palette import COL, METH
S = [l.rstrip('\n').split('\t') for l in open(f'{A}/samples{SFX}.tsv') if l.strip()]
POOL = {}
for lab, bam, m in S:
    z = np.load(f'{POS}/{lab}.npz')
    for k in ('gene', 't', 'L', 'alen', 'mate'): POOL.setdefault(m, {}).setdefault(k, []).append(z[k])
    if 'tx_len' not in globals(): tx_len = z['tx_len']; genes = z['genes']
for m in METH: POOL[m] = {k: np.concatenate(v) for k, v in POOL[m].items()}
def by_gene(d):
    g = d['gene']; o = np.argsort(g, kind='stable'); gs = g[o]; st = np.flatnonzero(np.r_[True, gs[1:] != gs[:-1]]); en = np.r_[st[1:], len(gs)]; return {int(gs[s]): o[s:e] for s, e in zip(st, en)}
IDX = {m: by_gene(POOL[m]) for m in METH}
def gene_cov(d, ii, L):
    t = d['t'][ii].astype(int); e = np.minimum(t + d['alen'][ii].astype(int), L); diff = np.zeros(L + 1); np.add.at(diff, t, 1); np.add.at(diff, e, -1); return np.cumsum(diff)[:L]
def common_genes(d0, lo=0, hi=10**9): return sorted(g for g in IDX[METH[0]] if all(len(IDX[m].get(g, [])) >= d0 for m in METH) and lo <= tx_len[g] < hi)
def covs(d0, lo, seed=11):
    r = np.random.default_rng(seed); common = common_genes(d0, lo); C = {m: [] for m in METH}
    for m in METH:
        for g in common:
            idx = IDX[m][g]; ii = r.choice(idx, d0, replace=False); C[m].append(gene_cov(POOL[m], ii, int(tx_len[g])))
    return common, C
def sm(y, k=25):
    """centred running mean, NaN-aware, for display only"""
    y = np.asarray(y, float); v = np.where(np.isnan(y), 0, y); n = (~np.isnan(y)).astype(float); ker = np.ones(k); return np.convolve(v, ker, 'same') / np.maximum(np.convolve(n, ker, 'same'), 1)
def covs_all(lo):
    """all unique reads of every gene (no subsampling) for the genes with >= 200 reads in every method"""
    common = common_genes(200, lo); C = {m: [gene_cov(POOL[m], IDX[m][g], int(tx_len[g])) for g in common] for m in METH}; return common, C
def stack(C, w, end):
    """genes x w matrix of coverage in the last (end=3) or first (end=5) w nt; genes shorter than w are NaN-padded"""
    X = np.full((len(C), w), np.nan)
    for k, c in enumerate(C):
        n = min(w, len(c)); X[k, :n] = (c[::-1] if end == 3 else c)[:n]
    return X
def rel_mean(X, C): mu = np.array([c.mean() for c in C]); return X / np.where(mu > 0, mu, 1)[:, None]
def rel_max(X, C): mx = np.array([c.max() for c in C]); return X / np.where(mx > 0, mx, 1)[:, None]
sets = {}
for d0, lo in ((200, 2000), (200, 4000), (200, 1000), (100, 2000), (400, 2000), (200, 8000)):
    sets[(d0, lo)] = covs(d0, lo); print(d0, lo, 'genes', len(sets[(d0, lo)][0]))
sets[('all', 2000)] = covs_all(2000); print('all reads >= 2 kb genes', len(sets[('all', 2000)][0]), 'median reads per gene per method', {m: int(np.median([len(IDX[m][g]) for g in sets[('all', 2000)][0]])) for m in METH})
V = []
def add(title, ylab, key, w, end, f, hline=None, ylim=None, log=False):
    common, C = sets[key]; D = {}
    for m in METH: X = stack(C[m], w, end); D[m] = f(X, C[m])
    V.append((title, ylab, D, w, end, hline, ylim, log, len(common)))
med = lambda X, C: np.nanmedian(rel_mean(X, C), 0); mean = lambda X, C: np.nanmean(rel_mean(X, C), 0)
add('1. p4m: coverage / gene mean, MEDIAN gene, transcripts >= 2 kb', 'coverage / gene mean', (200, 2000), 3000, 3, med, hline=1)
add('2. same, MEAN over genes', 'coverage / gene mean', (200, 2000), 3000, 3, mean, hline=1)
add('3. median, transcripts >= 4 kb, 4 kb window', 'coverage / gene mean', (200, 4000), 4000, 3, med, hline=1)
add('4. median, transcripts >= 8 kb, 6 kb window', 'coverage / gene mean', (200, 8000), 6000, 3, med, hline=1)
add('5. median, transcripts >= 1 kb, 1.5 kb window', 'coverage / gene mean', (200, 1000), 1500, 3, med, hline=1)
add('6. coverage / gene PEAK, median gene, >= 2 kb', 'coverage / gene peak', (200, 2000), 3000, 3, lambda X, C: np.nanmedian(rel_max(X, C), 0), ylim=(0, 1))
add('7. raw depth (reads per base at 200 reads), median gene, >= 2 kb', 'reads per base', (200, 2000), 3000, 3, lambda X, C: np.nanmedian(X, 0))
add('8. % of genes with >= 1 read at this distance, >= 2 kb', '% of genes covered', (200, 2000), 3000, 3, lambda X, C: 100 * np.nanmean(np.where(np.isnan(X), np.nan, X > 0), 0), ylim=(0, 100))
add('9. % of genes with >= 3 reads at this distance, >= 2 kb', '% of genes', (200, 2000), 3000, 3, lambda X, C: 100 * np.nanmean(np.where(np.isnan(X), np.nan, X >= 3), 0), ylim=(0, 100))
add('10. as 1 but from the CAP (5\' end), >= 2 kb', 'coverage / gene mean', (200, 2000), 3000, 5, med, hline=1)
add('11. as 1 at 100 reads per gene', 'coverage / gene mean', (100, 2000), 3000, 3, med, hline=1)
add('12. as 1 at 400 reads per gene', 'coverage / gene mean', (400, 2000), 3000, 3, med, hline=1)
fig, axs = plt.subplots(3, 4, figsize=(17, 11.8)); axs = axs.ravel()
for a, (title, ylab, D, w, end, hline, ylim, log, n) in zip(axs, V):
    xx = np.arange(w)
    for m in METH: a.plot(xx, D[m], color=COL[m], lw=1.8)
    if hline is not None: a.axhline(hline, color='#888', lw=0.7, ls=':')
    if ylim: a.set_ylim(*ylim)
    a.set_title(textwrap.fill(f'{title} ({n} genes)', 60), fontsize=8.5); a.set_ylabel(ylab, fontsize=8); a.tick_params(labelsize=7); a.grid(alpha=0.25)
    if end == 3: a.set_xlim(w, 0); a.set_xlabel('distance from the poly(A) site (nt)', fontsize=8)
    else: a.set_xlim(0, w); a.set_xlabel('distance from the cap / 5\' end (nt)', fontsize=8)
fig.legend(handles=[Line2D([], [], color=COL[m], lw=2.5) for m in METH], labels=METH, loc='lower center', ncol=4, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.045))
cap = textwrap.fill('TESTED: variations of the unscaled coverage profile (p4m). Genes with >= 200 (or 100 / 400) unique reads in every method and a canonical mature transcript at least as long as stated, each gene subsampled to exactly that many reads; '
    'per-base coverage along the canonical transcript from each read\'s 5\' position and aligned length (' + BOBNOTE + '); x = nucleotides from the poly(A) site (or from the cap in panel 10), so no scaling and no window that shrinks with '
    'length; y = coverage divided by the gene\'s own mean (1 = uniform), or by its peak, or raw reads per base, or the % of genes with coverage at that distance; the line is the median (or mean) over genes with equal weight. Genes shorter than the window contribute only where they exist.', 300)
fig.text(0.01, 0.005, cap, fontsize=7.6, va='bottom', color='#333'); fig.tight_layout(rect=(0, 0.08, 1, 1)); fig.savefig(f'{OUT}/Q_unscaled_variants.png', dpi=140); fig.savefig(f'{OUT}/Q_unscaled_variants.svg'); plt.close(fig); print('wrote Q sheet')
V2 = []
def add2(title, ylab, key, w, end, f, hline=None, ylim=None):
    common, C = sets[key]; D = {m: f(stack(C[m], w, end), C[m]) for m in METH}; V2.append((title, ylab, D, w, end, hline, ylim, len(common)))
pct1 = lambda X, C: 100 * np.nanmean(np.where(np.isnan(X), np.nan, X > 0), 0)
add2('1. p4m smoothed (25-nt running mean): coverage / gene mean, median gene, >= 2 kb', 'coverage / gene mean', (200, 2000), 3000, 3, lambda X, C: sm(med(X, C)), hline=1)
add2('2. transcripts >= 1 kb, 1.5 kb window, median gene, smoothed', 'coverage / gene mean', (200, 1000), 1500, 3, lambda X, C: sm(med(X, C)), hline=1)
add2('3. from the cap: coverage / gene mean, median gene, >= 2 kb, smoothed', 'coverage / gene mean', (200, 2000), 3000, 5, lambda X, C: sm(med(X, C)), hline=1)
add2('4. from the cap: % of genes with >= 1 read, >= 2 kb', '% of genes covered', (200, 2000), 3000, 5, pct1, ylim=(0, 100))
add2('5. ALL unique reads (no subsampling): coverage / gene mean, median gene, >= 2 kb', 'coverage / gene mean', ('all', 2000), 3000, 3, lambda X, C: sm(med(X, C)), hline=1)
add2('6. ALL unique reads: % of genes with >= 1 read, >= 2 kb', '% of genes covered', ('all', 2000), 3000, 3, pct1, ylim=(0, 100))
add2('7. ALL unique reads, from the cap: % of genes with >= 1 read, >= 2 kb', '% of genes covered', ('all', 2000), 3000, 5, pct1, ylim=(0, 100))
add2('8. ALL unique reads: % of genes with >= 3 reads, >= 2 kb', '% of genes', ('all', 2000), 3000, 3, lambda X, C: 100 * np.nanmean(np.where(np.isnan(X), np.nan, X >= 3), 0), ylim=(0, 100))
fig, axs = plt.subplots(2, 4, figsize=(17, 8.2)); axs = axs.ravel()
for a, (title, ylab, D, w, end, hline, ylim, n) in zip(axs, V2):
    xx = np.arange(w)
    for m in METH: a.plot(xx, D[m], color=COL[m], lw=1.8)
    if hline is not None: a.axhline(hline, color='#888', lw=0.7, ls=':')
    if ylim: a.set_ylim(*ylim)
    a.set_title(textwrap.fill(f'{title} ({n} genes)', 58), fontsize=8); a.set_ylabel(ylab, fontsize=8); a.tick_params(labelsize=7); a.grid(alpha=0.25)
    if end == 3: a.set_xlim(w, 0); a.set_xlabel('distance from the poly(A) site (nt)', fontsize=8)
    else: a.set_xlim(0, w); a.set_xlabel('distance from the cap / 5\' end (nt)', fontsize=8)
fig.legend(handles=[Line2D([], [], color=COL[m], lw=2.5) for m in METH], labels=METH, loc='lower center', ncol=4, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.06))
cap2 = textwrap.fill('TESTED: second round of the unscaled view. Panels 1-3 are the Q-sheet profiles with a 25-nt running mean for display; panel 4 anchors at the cap and counts the % of genes with any read at that distance; panels 5-8 use ALL unique reads of every gene '
    '(no subsampling; genes with >= 200 reads in every method and a canonical transcript >= 2 kb) to show what the browser shows at native depth: if the 3\' methods reach the transcript body at all, more reads would show it here.', 300)
fig.text(0.01, 0.005, cap2, fontsize=7.6, va='bottom', color='#333'); fig.tight_layout(rect=(0, 0.1, 1, 1)); fig.savefig(f'{OUT}/Q2_unscaled_variants.png', dpi=140); fig.savefig(f'{OUT}/Q2_unscaled_variants.svg'); plt.close(fig); print('wrote Q2 sheet')
with open(f'{OUT}/Q_UNSCALED.md', 'w') as f:
    f.write('## Q. Unscaled profiles (nt from the poly(A) site or the cap), equal gene weight\n\n| variant | genes | method | at 100 nt | at 500 nt | at 1,000 nt | at 2,000 nt | at 3,000 nt (or window end) | max |\n|---|---|---|---|---|---|---|---|---|\n')
    for title, ylab, D, w, end, *_, n in V + [('Q2 ' + t, yl, D, w, e, h, yl2, n) for t, yl, D, w, e, h, yl2, n in V2]:
        for m in METH:
            y = D[m]; pk = lambda i: f'{y[min(i, w - 1)]:.2f}'
            f.write(f'| {title} | {n} | {m} | {pk(100)} | {pk(500)} | {pk(1000)} | {pk(2000)} | {pk(w - 1)} | {np.nanmax(y):.2f} |\n')
print(open(f'{OUT}/Q_UNSCALED.md').read())
# ---------------- preprint-style panels ----------------
plt.rcParams.update({'svg.fonttype': 'none', 'font.family': ['Arial', 'Nimbus Sans', 'DejaVu Sans'], 'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6.5,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.direction': 'out', 'ytick.direction': 'out', 'axes.spines.top': False, 'axes.spines.right': False,
    'text.color': '#231F20', 'axes.edgecolor': '#231F20', 'axes.labelcolor': '#231F20', 'xtick.color': '#231F20', 'ytick.color': '#231F20', 'axes.grid': False})
def psave(fig, name):
    fig.tight_layout(pad=0.4); fig.savefig(f'{OUTP}/{name}.svg', bbox_inches='tight', pad_inches=0.03); fig.savefig(f'{OUTP}/{name}.png', dpi=200, bbox_inches='tight', pad_inches=0.03); plt.close(fig)
    s = open(f'{OUTP}/{name}.svg').read(); s = re.sub(r"font-family:\s*[^;\"]+", 'font-family: Arial', s); open(f'{OUTP}/{name}.svg', 'w').write(s); print('wrote', name)
LEGEND = os.environ.get('BM_PANEL_LEGEND', '0') == '1'   # legends off by default; one shared legend SVG (panel_legends.py) instead
def note(ax, *lines, loc='tl'):
    """small grey in-axes note stating the set and the matching rule (no titles)"""
    if os.environ.get('BM_PANEL_NOTES', '0') != '1': return   # no notes; matching lives in the axis labels only
    x, y, ha, va = {'tl': (0.03, 0.97, 'left', 'top'), 'tr': (0.97, 0.97, 'right', 'top'), 'bl': (0.03, 0.04, 'left', 'bottom'), 'br': (0.97, 0.04, 'right', 'bottom'), 'tc': (0.5, 0.97, 'center', 'top'), 'ml': (0.03, 0.42, 'left', 'center'), 'mr': (0.97, 0.42, 'right', 'center')}[loc]
    ax.text(x, y, '\n'.join(lines), transform=ax.transAxes, fontsize=5.5, color='#6B6B6B', ha=ha, va=va, linespacing=1.3, zorder=5)
SETLBL = 'native read length (BOBseq 2 x 150 PE)' if SET == 'native' else '50-nt single reads, all methods'; NCOV = f"{sum(1 for s_ in S if s_[2] == 'BOBseq')} BOBseq, 3 per competitor"
def plegend(fig, ncol=4):
    if not LEGEND: return
    fig.legend(handles=[Line2D([], [], color=COL[m], lw=1.6) for m in METH], labels=METH, loc='lower center', frameon=False, ncol=ncol, handlelength=1.4, columnspacing=1.0, bbox_to_anchor=(0.5, -0.06))
def panel(name, D, w, end, ylab, hline=1, ylim=None, size=(2.6, 2.2), n=None, reads='200 unique reads each', loc=None):
    fig, ax = plt.subplots(figsize=size); xx = np.arange(w)
    for m in METH: ax.plot(xx, D[m], color=COL[m], lw=1.6)
    if hline is not None: ax.axhline(hline, color='#9a9a9a', lw=0.6, ls=(0, (2, 2)))
    if end == 3: ax.set_xlim(w, 0); ax.set_xlabel('distance from the poly(A) site (nt)')
    else: ax.set_xlim(0, w); ax.set_xlabel("distance from the 5' end (nt)")
    ax.set_ylim(*(ylim or (0, None))); ax.set_ylabel(ylab); note(ax, SETLBL, f'{n} genes, {reads}', NCOV, loc=loc or ('tl' if end == 3 else 'tr')); plegend(fig); psave(fig, name)
byname = {t.split('.')[0]: (D, w, end) for t, _, D, w, end, *_ in V}; by2 = {t.split('.')[0]: (D, w, end) for t, _, D, w, end, *_ in V2}
def smoothed(D): return {m: sm(D[m]) for m in D}
N2, N1, NA = len(sets[(200, 2000)][0]), len(sets[(200, 1000)][0]), len(sets[('all', 2000)][0])
panel('p4m_coverage_vs_nt_from_polyA', smoothed(byname['1'][0]), 3000, 3, 'coverage relative to gene mean\n(median gene, transcripts >= 2 kb)', n=N2)
panel('p4m2_coverage_vs_nt_from_polyA_mean', smoothed(byname['2'][0]), 3000, 3, 'coverage relative to gene mean\n(mean over genes, transcripts >= 2 kb)', n=N2)
panel('p4m3_coverage_vs_nt_from_polyA_1kb', smoothed(byname['5'][0]), 1500, 3, 'coverage relative to gene mean\n(median gene, transcripts >= 1 kb)', n=N1)
panel('p4n_genes_covered_vs_nt_from_polyA', *byname['8'], 'genes covered at this distance (%)\n(200 reads per gene, transcripts >= 2 kb)', hline=None, ylim=(0, 100), n=N2)
panel('p4n2_depth_vs_nt_from_polyA', smoothed(byname['7'][0]), 3000, 3, 'reads per base, median gene\n(200 reads, transcripts >= 2 kb)', hline=None, n=N2)
panel('p4o_coverage_vs_nt_from_cap', smoothed(byname['10'][0]), 3000, 5, 'coverage relative to gene mean\n(median gene, transcripts >= 2 kb)', n=N2, ylim=(0, 1.3))
panel('p4o2_genes_covered_vs_nt_from_cap', *by2['4'], 'genes covered at this distance (%)\n(200 reads per gene, transcripts >= 2 kb)', hline=None, ylim=(0, 100), n=N2, loc='br')
panel('p4q_allreads_coverage_vs_nt_from_polyA', *by2['5'], 'coverage relative to gene mean\n(median gene, all reads, transcripts >= 2 kb)', n=NA, reads='all unique reads (no subsampling)')
panel('p4q2_allreads_genes_covered_vs_nt_from_polyA', *by2['6'], 'genes covered at this distance (%)\n(all reads, transcripts >= 2 kb)', hline=None, ylim=(0, 100), n=NA, reads='all unique reads (no subsampling)', loc='ml')
panel('p4q3_allreads_genes_covered_vs_nt_from_cap', *by2['7'], 'genes covered at this distance (%)\n(all reads, transcripts >= 2 kb)', hline=None, ylim=(0, 100), n=NA, reads='all unique reads (no subsampling)', loc='ml')
# both ends in one panel: first 1 kb from the cap | last 1 kb to the poly(A) site, transcripts >= 2 kb (no overlap)
def both_ends(name, key, w, f, ylab, hline=1, ylim=None):
    common, C = sets[key]; fig, (a1, a2) = plt.subplots(1, 2, figsize=(3.2, 2.35), sharey=True, gridspec_kw={'wspace': 0.12})
    for m in METH: a1.plot(np.arange(w), f(stack(C[m], w, 5), C[m]), color=COL[m], lw=1.6); a2.plot(np.arange(w), f(stack(C[m], w, 3), C[m]), color=COL[m], lw=1.6)
    for a in (a1, a2):
        if hline is not None: a.axhline(hline, color='#9a9a9a', lw=0.6, ls=(0, (2, 2)))
        a.set_ylim(*(ylim or (0, None)))
    a1.set_xlim(0, w); a2.set_xlim(w, 0); a1.set_xticks([0, w // 2, w]); a2.set_xticks([w, w // 2, 0]); a1.set_xticklabels(['0', f'{w//2}', '']); a2.set_xticklabels(['', f'{w//2}', '0'])
    a1.set_xlabel("nt from the 5' end", labelpad=2); a2.set_xlabel('nt to the poly(A) site', labelpad=2); a1.set_ylabel(ylab); a2.spines['left'].set_visible(False); a2.tick_params(axis='y', length=0)
    note(a1, SETLBL.split(' (')[0], f'{len(common)} genes, ' + ('all unique reads' if key[0] == 'all' else '200 reads each'), NCOV, loc='tr'); plegend(fig); psave(fig, name)
both_ends('p4p_coverage_both_ends', (200, 2000), 1000, lambda X, C: sm(med(X, C)), 'coverage relative to gene mean\n(median gene, transcripts >= 2 kb)')
both_ends('p4p2_genes_covered_both_ends', (200, 2000), 1000, pct1, 'genes covered at this distance (%)\n(200 reads per gene, transcripts >= 2 kb)', hline=None, ylim=(0, 100))
both_ends('p4p3_allreads_genes_covered_both_ends', ('all', 2000), 1000, pct1, 'genes covered at this distance (%)\n(all reads, transcripts >= 2 kb)', hline=None, ylim=(0, 100))
# ---- all reads, no subsampling: p4m2 counterpart (mean over genes) and p4i counterpart (heatmap) ----
common_a, Ca = sets[('all', 2000)]
panel('p4q4_allreads_coverage_vs_nt_from_polyA_mean', smoothed({m: np.nanmean(rel_mean(stack(Ca[m], 3000, 3), Ca[m]), 0) for m in METH}), 3000, 3, 'coverage relative to gene mean\n(mean over genes, all reads, >= 2 kb)', n=NA, reads='all unique reads (no subsampling)')
common_h = sorted(common_genes(200), key=lambda g: tx_len[g]); nb = 50; H = {m: np.zeros((len(common_h), nb)) for m in METH}; nreads = {}
for m in METH:
    nreads[m] = int(np.median([len(IDX[m][g]) for g in common_h]))
    for r, g in enumerate(common_h):
        L = int(tx_len[g]); cov = gene_cov(POOL[m], IDX[m][g], L); idx = np.minimum((np.arange(L) * nb) // L, nb - 1); b = np.bincount(idx, weights=cov, minlength=nb) / np.bincount(idx, minlength=nb); H[m][r] = b / b.max() if b.max() > 0 else b
fig, ax = plt.subplots(1, len(METH), figsize=(1.3 * len(METH), 2.4), sharey=True)   # one panel per method
for a, m in zip(ax, METH):
    a.imshow(H[m], aspect='auto', cmap='Greys', vmin=0, vmax=1, interpolation='nearest'); a.set_title(m, fontsize=7, color=COL[m]); a.set_xticks([0, nb - 1]); a.set_xticklabels(["5'", "3'"]); a.tick_params(length=2)
    for sp in a.spines.values(): sp.set_visible(True); sp.set_linewidth(0.5)
ax[0].set_ylabel(f'{len(common_h):,} genes (short to long)\nall unique reads per gene'); ax[0].set_yticks([]); fig.text(0.5, -0.04, 'coverage along the mature transcript (each row scaled to its own peak)', ha='center', va='top', fontsize=7); psave(fig, 'p4i2_coverage_heatmaps_all_reads')
json.dump({'genes': len(common_h), 'median_unique_reads_per_gene': nreads, 'frac_bins_above_20pct_of_peak': {m: float(np.mean(H[m] > 0.2)) for m in METH}}, open(f'{OUTP}/p4i2_all_reads_numbers.json', 'w'), indent=1); print('p4i2', nreads)
print('DONE')

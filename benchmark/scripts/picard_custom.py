#!/usr/bin/env python3
"""Custom Picard-style gene-body profile built from base coverage at matched reads per gene (canonical transcript, unique reads, BOBseq both mates).
Unlike Picard: one canonical transcript per gene, every gene subsampled to the same number of unique reads, per-gene normalization, equal gene weight.
Writes a variant sheet (preprint_figures_native/coverage_architecture/P_picard_custom_variants.png) + numbers, and preprint-style panels p4j/p4k/p4l."""
import os, sys, json, textwrap, re, logging, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt; logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
from matplotlib.lines import Line2D
from settings import WORK, COVERAGE
A = COVERAGE; P = WORK
SET = os.environ.get('BM_COVSET', 'native'); SFX = '' if SET == 'native' else f'_{SET}'; FIG = 'preprint_figures_native' if SET == 'native' else f'preprint_figures_{SET}'
BOBNOTE = 'BOBseq: both mates' if SET == 'native' else 'every method one 50-nt single read per fragment'; READNOTE = ('BOBseq both mates of the 2 x 150 pairs; competitors one native-length read per fragment' if SET == 'native' else 'read-length matched set: every method one 50-nt single read per fragment (BOBseq R2 truncated to 50 nt)')
POS = f'{A}/positions_canonical{SFX}'; OUT = f'{P}/{FIG}/coverage_architecture'; OUTP = f'{P}/{FIG}/main_figure_panels'; os.makedirs(OUT, exist_ok=True); os.makedirs(OUTP, exist_ok=True)
from palette import COL, METH
S = [l.rstrip('\n').split('\t') for l in open(f'{A}/samples{SFX}.tsv') if l.strip()]; rng = np.random.default_rng(11)
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
NB = 100
def binned(cov, nb=NB):
    L = len(cov); idx = np.minimum((np.arange(L) * nb) // L, nb - 1); return np.bincount(idx, weights=cov, minlength=nb) / np.bincount(idx, minlength=nb)
def profiles(d0, lo=0, hi=10**9, seed=11):
    """per method: matrix genes x NB of binned base coverage at d0 reads per gene (raw depth), plus 3'-anchored unscaled coverage (0..3000 nt from the 3' end)"""
    r = np.random.default_rng(seed); common = common_genes(d0, lo, hi); M = {m: np.zeros((len(common), NB)) for m in METH}; U = {m: np.full((len(common), 3000), np.nan) for m in METH}
    for m in METH:
        for k, g in enumerate(common):
            idx = IDX[m][g]; ii = r.choice(idx, d0, replace=False); L = int(tx_len[g]); cov = gene_cov(POOL[m], ii, L); M[m][k] = binned(cov); w = min(3000, L); U[m][k, :w] = cov[::-1][:w]
    return common, M, U
def norm_mean(M): mu = M.mean(1, keepdims=True); return M / np.where(mu > 0, mu, 1)
def norm_max(M): mx = M.max(1, keepdims=True); return M / np.where(mx > 0, mx, 1)
x = (np.arange(NB) + 0.5)
# ---------------- variant sheet ----------------
common, M, U = profiles(200); n = len(common); print('genes', n)
commonL, ML, UL = profiles(200, lo=2000); print('genes >= 2 kb', len(commonL))
common4, M4, _ = profiles(400); print('genes at 400', len(common4))
common1, M1, _ = profiles(100); print('genes at 100', len(common1))
V = []   # (title, ylabel, dict method -> (y, lo, hi) , kwargs)
def add(title, ylab, f, M_=None, hline=None, ylim=None, xlab=None, xs=None):
    V.append((title, ylab, {m: f(M_[m] if M_ is not None else M[m]) for m in METH}, hline, ylim, xlab, xs))
add(f'1. Picard-like: per-gene coverage / gene mean, MEAN over {n} genes', 'coverage / gene mean', lambda A_: (norm_mean(A_).mean(0), None, None), hline=1)
add(f'2. per-gene coverage / gene mean, MEDIAN over genes', 'coverage / gene mean', lambda A_: (np.median(norm_mean(A_), 0), None, None), hline=1)
add(f'3. median with interquartile band over genes', 'coverage / gene mean', lambda A_: (np.median(norm_mean(A_), 0), np.percentile(norm_mean(A_), 25, 0), np.percentile(norm_mean(A_), 75, 0)), hline=1)
add(f'4. per-gene coverage / gene MAXIMUM, mean over genes', 'coverage / gene peak', lambda A_: (norm_max(A_).mean(0), None, None), ylim=(0, 1))
add(f'5. % of genes with any coverage at this position', '% of genes covered', lambda A_: (100 * (A_ > 0).mean(0), None, None), ylim=(0, 100))
add(f'6. % of genes with coverage >= 20% of their peak at this position', '% of genes', lambda A_: (100 * (norm_max(A_) >= 0.2).mean(0), None, None), ylim=(0, 100))
add(f'7. % of genes with coverage >= half of their mean', '% of genes', lambda A_: (100 * (norm_mean(A_) >= 0.5).mean(0), None, None), ylim=(0, 100))
add(f'8. log2(coverage / gene mean), median over genes', 'log2 coverage / gene mean', lambda A_: (np.median(np.log2(np.maximum(norm_mean(A_), 1 / 64)), 0), None, None), hline=0)
add(f'9. as 1 (mean), transcripts >= 2 kb only ({len(commonL)} genes)', 'coverage / gene mean', lambda A_: (norm_mean(A_).mean(0), None, None), M_=ML, hline=1)
add(f'10. as 1 (mean), at 100 reads per gene ({len(common1)} genes)', 'coverage / gene mean', lambda A_: (norm_mean(A_).mean(0), None, None), M_=M1, hline=1)
add(f'11. as 1 (mean), at 400 reads per gene ({len(common4)} genes)', 'coverage / gene mean', lambda A_: (norm_mean(A_).mean(0), None, None), M_=M4, hline=1)
# unscaled: coverage / gene mean by distance from the 3' end in nt (genes >= 2 kb so every gene spans the window)
def unscaled(Uu, Mm, agg=np.nanmedian):
    mu = Mm.mean(1); Z = Uu / np.where(mu > 0, mu, 1)[:, None]; return agg(Z, 0), None, None
add(f'12. unscaled: coverage / gene mean vs nt from the 3\' end, MEAN, transcripts >= 2 kb', 'coverage / gene mean', lambda A_: unscaled(UL[A_], ML[A_], np.nanmean), M_={m: m for m in METH}, hline=1, xlab='distance from the 3\' end (nt)', xs=np.arange(3000))
fig, axs = plt.subplots(3, 4, figsize=(17, 11)); axs = axs.ravel()
for a, (title, ylab, D, hline, ylim, xlab, xs) in zip(axs, V):
    for m in METH:
        y, lo, hi = D[m]; xx = xs if xs is not None else x; a.plot(xx, y, color=COL[m], lw=1.8)
        if lo is not None: a.fill_between(xx, lo, hi, color=COL[m], alpha=0.12, lw=0)
    if hline is not None: a.axhline(hline, color='#888', lw=0.7, ls=':')
    if ylim: a.set_ylim(*ylim)
    a.set_title(title, fontsize=8.5); a.set_ylabel(ylab, fontsize=8); a.set_xlabel(xlab or "position along the mature transcript, 5' to 3' (%)", fontsize=8); a.tick_params(labelsize=7); a.grid(alpha=0.25)
    if xs is not None: a.invert_xaxis()
fig.legend(handles=[Line2D([], [], color=COL[m], lw=2.5) for m in METH], labels=METH, loc='lower center', ncol=4, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.045))
cap = textwrap.fill(f'TESTED: a Picard-style gene-body profile that keeps what made the 200-read plots work. Genes with >= 200 unique reads in every method ({n} genes; same genes in every panel unless stated), each gene subsampled to exactly 200 reads, '
    'per-base coverage along the Ensembl canonical mature transcript from each read\'s 5\' position and aligned length (' + BOBNOTE + '), 100 bins from 5\' to 3\'. Panels differ only in the per-gene normalization (gene mean as Picard, or gene peak), '
    'the aggregation over genes (mean, median, band, or the % of genes covered at that position), the gene set (length, depth) and the x axis (scaled % or nucleotides from the 3\' end). Equal gene weight everywhere; no isoform averaging; unique reads (MAPQ 255).', 300)
fig.text(0.01, 0.005, cap, fontsize=7.6, va='bottom', color='#333'); fig.tight_layout(rect=(0, 0.08, 1, 1)); fig.savefig(f'{OUT}/P_picard_custom_variants.png', dpi=140); fig.savefig(f'{OUT}/P_picard_custom_variants.svg'); plt.close(fig); print('wrote variants')
# numbers
num = {}
for title, ylab, D, *_ in V:
    num[title] = {m: {'5p_10pct': float(np.mean(D[m][0][:10])), 'mid': float(np.mean(D[m][0][45:55])), '3p_10pct': float(np.mean(D[m][0][-10:])), 'max': float(np.max(D[m][0])), 'min': float(np.min(D[m][0]))} if len(D[m][0]) == NB else {'max': float(np.nanmax(D[m][0]))} for m in METH}
json.dump(num, open(f'{OUT}/P_picard_custom_numbers.json', 'w'), indent=1)
with open(f'{OUT}/P_PICARD_CUSTOM.md', 'w') as f:
    f.write(f'## P. Custom Picard-style profiles, {n} genes at 200 reads (canonical transcript, base coverage, equal gene weight)\n\n| variant | method | 5\'-most 10% | middle 10% | 3\'-most 10% | max | min |\n|---|---|---|---|---|---|---|\n')
    for title, d in num.items():
        for m in METH:
            v = d[m]
            if '5p_10pct' in v: f.write(f'| {title} | {m} | {v["5p_10pct"]:.2f} | {v["mid"]:.2f} | {v["3p_10pct"]:.2f} | {v["max"]:.2f} | {v["min"]:.2f} |\n')
print(open(f'{OUT}/P_PICARD_CUSTOM.md').read())
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
def plegend(fig):
    if not LEGEND: return
    fig.legend(handles=[Line2D([], [], color=COL[m], lw=1.6) for m in METH], labels=METH, loc='lower center', frameon=False, ncol=4, handlelength=1.4, columnspacing=1.0, bbox_to_anchor=(0.5, -0.06))
XT = [0, 25, 50, 75, 100]; XL = ["5'", '25', '50', '75', "3'"]
# p4j: median coverage / gene mean (Picard axis)
fig, ax = plt.subplots(figsize=(2.6, 2.2))
for m in METH: ax.plot(x, norm_mean(M[m]).mean(0), color=COL[m], lw=1.6)
ax.axhline(1, color='#9a9a9a', lw=0.6, ls=(0, (2, 2))); ax.set_xlim(0, 100); ax.set_ylim(bottom=0); ax.set_xticks(XT); ax.set_xticklabels(XL); ax.set_xlabel('position along the mature transcript (%)'); ax.set_ylabel('coverage relative to gene mean\n(mean over genes, 200 reads per gene)'); note(ax, SETLBL, f'{n} genes at 200 unique reads each', NCOV, loc='tc'); plegend(fig); psave(fig, 'p4j_coverage_relative_to_gene_mean')
# p4k: % of genes covered at this position
fig, ax = plt.subplots(figsize=(2.6, 2.2))
for m in METH: ax.plot(x, 100 * (M[m] > 0).mean(0), color=COL[m], lw=1.6)
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.set_xticks(XT); ax.set_xticklabels(XL); ax.set_xlabel('position along the mature transcript (%)'); ax.set_ylabel('genes covered at this position (%)\n(200 reads per gene)'); note(ax, SETLBL.split(' (')[0], f'{n} genes at 200 reads each', NCOV, loc='br'); plegend(fig); psave(fig, 'p4k_genes_covered_along_transcript')
# p4l: median with IQR band, coverage / gene mean
fig, ax = plt.subplots(figsize=(2.6, 2.2))
for m in METH:
    Z = norm_mean(M[m]); ax.fill_between(x, np.percentile(Z, 25, 0), np.percentile(Z, 75, 0), color=COL[m], alpha=0.15, lw=0); ax.plot(x, np.median(Z, 0), color=COL[m], lw=1.6)
ax.axhline(1, color='#9a9a9a', lw=0.6, ls=(0, (2, 2))); ax.set_xlim(0, 100); ax.set_ylim(bottom=0); ax.set_xticks(XT); ax.set_xticklabels(XL); ax.set_xlabel('position along the mature transcript (%)'); ax.set_ylabel('coverage relative to gene mean\n(median and IQR over genes)'); note(ax, SETLBL, f'{n} genes at 200 unique reads each', NCOV, loc='tc'); plegend(fig); psave(fig, 'p4l_coverage_relative_to_gene_mean_iqr')
print('DONE')

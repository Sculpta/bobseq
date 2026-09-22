#!/usr/bin/env python3
"""Main-figure panels for Illustrator, one SVG per panel: Arial 7 pt, thin black axes, no grid, no titles, no notes.
Data: matched (50 nt) and native rarefaction tables, native composition, coverage_architecture positions (canonical transcripts),
per-sample junctions at matched depth. Output: preprint_figures_native/main_figure_panels/<panel>.svg (+ .png preview).
"""

from settings import WORK, COVERAGE
import os, csv, json, glob, collections, numpy as np, matplotlib
from scipy import stats as _st

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

P = WORK
U = f'{P}/benchmark_uniform'
A = COVERAGE
SET = os.environ.get('BM_COVSET', 'native')
SFX = '' if SET == 'native' else f'_{SET}'
FIG = 'preprint_figures_native' if SET == 'native' else f'preprint_figures_{SET}'
BOBNOTE = 'BOBseq: both mates' if SET == 'native' else 'every method one 50-nt single read per fragment'
READNOTE = (
    'BOBseq both mates of the 2 x 150 pairs; competitors one native-length read per fragment'
    if SET == 'native'
    else 'read-length matched set: every method one 50-nt single read per fragment (BOBseq R2 truncated to 50 nt)'
)
OUT = f'{P}/{FIG}/main_figure_panels'
CV = f'{P}/{FIG}/coverage_architecture'
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update(
    {
        'svg.fonttype': 'none',
        'font.family': ['Arial', 'Nimbus Sans', 'DejaVu Sans'],
        'font.size': 7,
        'axes.labelsize': 7,
        'xtick.labelsize': 6.5,
        'ytick.labelsize': 6.5,
        'legend.fontsize': 6.5,
        'axes.linewidth': 0.6,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 2.5,
        'ytick.major.size': 2.5,
        'xtick.minor.size': 1.5,
        'ytick.minor.size': 1.5,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'axes.spines.top': False,
        'axes.spines.right': False,
        'text.color': '#231F20',
        'axes.edgecolor': '#231F20',
        'axes.labelcolor': '#231F20',
        'xtick.color': '#231F20',
        'ytick.color': '#231F20',
    }
)
import os as _os, sys as _sys

_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from palette import COL, METH  # one palette shared by every panel

SETS = {
    'matched': {'DRUG-seq': 'drugseq', 'BRB-seq': 'brbseq', 'prime-seq': 'primeseq', 'BOBseq': 'bob57_24plex_pe_hs'},
    'native': {
        'DRUG-seq': 'drugseq_native',
        'BRB-seq': 'brbseq_native',
        'prime-seq': 'primeseq_native',
        'BOBseq': 'bob57_24plex_pe_native',
    },
}
SZ = (2.6, 2.2)
COMPSET = 'native' if SET == 'native' else 'matched'
if SET != 'native':
    SETS = {'matched': SETS['matched']}  # the 50-nt directory holds only read-length-matched panels


def save(fig, name):
    """figure panels (panel_keep.py) go to the panel directory, everything else straight to extra/, so no legacy panels appear in the main dir"""
    OUT = globals()['OUT'] if name in KEEP else globals()['OUT'] + '/extra'
    fig.tight_layout(pad=0.4)
    fig.savefig(f'{OUT}/{name}.svg', bbox_inches='tight', pad_inches=0.03)
    fig.savefig(f'{OUT}/{name}.png', dpi=200, bbox_inches='tight', pad_inches=0.03)
    plt.close(fig)
    s = open(f'{OUT}/{name}.svg').read()
    import re

    s = re.sub(r"font-family:\s*[^;\"]+", 'font-family: Arial', s)
    open(f'{OUT}/{name}.svg', 'w').write(s)
    print('wrote', name)


LEGEND = (
    os.environ.get('BM_PANEL_LEGEND', '0') == '1'
)  # legends off by default; one shared legend SVG (panel_legends.py) is drawn instead


def note(ax, *lines, loc='tl'):
    """small grey in-axes note stating the set and the matching rule (off by default: no titles, no notes)"""
    if os.environ.get('BM_PANEL_NOTES', '0') != '1':
        return  # notes off: the matching rule lives in the axis labels only
    x, y, ha, va = {
        'tl': (0.03, 0.97, 'left', 'top'),
        'tr': (0.97, 0.97, 'right', 'top'),
        'bl': (0.03, 0.04, 'left', 'bottom'),
        'br': (0.97, 0.04, 'right', 'bottom'),
        'tc': (0.5, 0.97, 'center', 'top'),
        'ml': (0.03, 0.42, 'left', 'center'),
        'mr': (0.97, 0.42, 'right', 'center'),
    }[loc]
    ax.text(
        x,
        y,
        '\n'.join(lines),
        transform=ax.transAxes,
        fontsize=5.5,
        color='#6B6B6B',
        ha=ha,
        va=va,
        linespacing=1.3,
        zorder=5,
    )


SETLBL = {'native': 'native read length (BOBseq 2 x 150 PE)', 'matched': '50-nt single reads, all methods'}


def legend(ax, loc=None, ncol=4):
    if not LEGEND:
        return
    ax.figure.legend(
        handles=[Line2D([], [], color=COL[m], lw=1.6) for m in METH],
        labels=METH,
        loc='lower center',
        frameon=False,
        ncol=ncol,
        handlelength=1.4,
        columnspacing=1.0,
        bbox_to_anchor=(0.5, -0.06),
    )  # legend below the axes; one shared legend is placed per figure in Illustrator


def rare(arm):
    return [
        r
        for r in csv.DictReader(open(f'{U}/{arm}/stats_E_U/rarefied.tsv'), delimiter='\t')
        if r.get('B_corr') and r['at_native'] == '0'
    ]


DEPTHS = [50000, 100000, 250000, 500000, 1000000, 2000000]


def band(
    ax, x, y, col, lw=1.3, alpha=0.3, dots=True
):  # thin line, stronger ribbon so the 95% interval reads
    """replicate dots (one per sample at each depth, light, small log-x jitter so they do not stack) under a thick mean line
    with the 95% t-interval over samples at each depth (depths with >= 3 samples)"""
    if dots:
        jit = np.exp(np.random.default_rng(7).normal(0, 0.03, len(x)))
        ax.scatter(x * jit, y, s=5, color=col, alpha=0.35, lw=0, zorder=1)
    ds = [d for d in sorted(set(x)) if (x == d).sum() >= 3]
    mu = np.array([y[x == d].mean() for d in ds])
    n = np.array([(x == d).sum() for d in ds])
    ci = np.array([_st.t.ppf(0.975, k - 1) * y[x == d].std(ddof=1) / np.sqrt(k) for d, k in zip(ds, n)])
    ax.fill_between(ds, mu - ci, mu + ci, color=col, alpha=alpha, lw=0, zorder=2)
    ax.plot(ds, mu, color=col, lw=lw, zorder=3)


from panel_keep import keep as _keep

KEEP = _keep(COMPSET)
os.makedirs(f'{OUT}/extra', exist_ok=True)


def nlegend(ax, arms, loc):
    """in-axes legend with the number of wells per method"""
    ax.legend(
        handles=[Line2D([], [], color=COL[m], lw=1.3) for m in METH],
        labels=[f"{m} (n={len({r['well'] for r in rare(arms[m])})})" for m in METH],
        loc=loc,
        frameon=False,
        fontsize=6.5,
        handlelength=1.6,
        borderaxespad=0.3,
        labelspacing=0.3,
    )


# ---------- 1. molecules (UCI) per sample vs filtered reads ----------
for setname, arms in SETS.items():
    fig, ax = plt.subplots(figsize=SZ)
    for (
        m
    ) in (
        METH
    ):  # thick line = mean over wells, ribbon = 95% t-interval of the mean, linear axis, no dots, in-axes legend with n
        rows = rare(arms[m])
        x = np.array([int(r['depth']) for r in rows])
        y = np.array([float(r['B_corr']) for r in rows])
        band(ax, x, y, COL[m], dots=False)
    ax.set_xscale('log')
    ax.set_xlabel('filtered reads per sample (subsampled)')
    ax.set_ylabel('molecules per sample (UMI)')
    ax.set_xticks(DEPTHS)
    ax.set_xticklabels(['50k', '100k', '250k', '500k', '1M', '2M'])
    ax.set_ylim(bottom=0)
    ax.set_yticks([0, 250000, 500000, 750000, 1000000])
    ax.set_yticklabels(['0', '250k', '500k', '750k', '1M'])
    nlegend(ax, arms, 'upper left')
    ns = ' / '.join(str(len({r['well'] for r in rare(arms[m])})) for m in METH)
    note(ax, SETLBL[setname], 'every sample subsampled to each depth', f'n = {ns} samples', loc='tl')
    save(fig, f'p1_molecules_vs_depth_{setname}')
    fig, ax = plt.subplots(figsize=SZ)  # log-axis version of the same panel (both versions are kept)
    for m in METH:
        rows = rare(arms[m])
        x = np.array([int(r['depth']) for r in rows])
        y = np.array([float(r['B_corr']) for r in rows])
        band(ax, x, y, COL[m], dots=False)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('filtered reads per sample (subsampled)')
    ax.set_ylabel('molecules per sample (UMI)')
    ax.set_xticks(DEPTHS)
    ax.set_xticklabels(['50k', '100k', '250k', '500k', '1M', '2M'])
    nlegend(ax, arms, 'upper left')
    save(fig, f'p1_molecules_vs_depth_log_{setname}')


# ---------- 2. genes per sample at matched depth ----------
def strip(ax, data, ylabel, log=False, xlabel=None):
    for i, m in enumerate(METH):
        v = np.array(data[m])
        jit = np.random.default_rng(1).normal(0, 0.07, len(v))
        ax.scatter(i + jit, v, s=7, color=COL[m], alpha=0.8, lw=0, zorder=2)
        ax.hlines(np.median(v), i - 0.28, i + 0.28, color='#231F20', lw=1.0, zorder=3)
    ax.set_xticks(range(len(METH)))
    ax.set_xticklabels(METH, rotation=30, ha='right')
    ax.set_ylabel(ylabel)
    ax.set_xlim(-0.6, len(METH) - 0.4)
    if xlabel and os.environ.get('BM_PANEL_NOTES', '0') == '1':
        ax.set_xlabel(xlabel, fontsize=5.5, color='#6B6B6B', labelpad=2)
    if log:
        ax.set_yscale('log')


for setname, arms in SETS.items():
    for d in (250000, 500000):
        data = {m: [float(r['genes_B_noRP']) for r in rare(arms[m]) if int(r['depth']) == d] for m in METH}
        fig, ax = plt.subplots(figsize=(2.2, 2.2))
        strip(
            ax,
            data,
            f'genes per sample at {d//1000}k filtered reads\n(ribosomal-protein genes excluded)',
            xlabel=SETLBL[setname],
        )
        ax.set_ylim(bottom=0)
        note(ax, f'each sample subsampled to {d//1000}k', 'filtered reads; (n) = samples', loc='bl')
        save(fig, f'p2_genes_at_{d//1000}k_{setname}')
        ax = None
    fig, ax = plt.subplots(figsize=SZ)
    for m in METH:
        rows = rare(arms[m])
        x = np.array([int(r['depth']) for r in rows])
        y = np.array([float(r['genes_B_noRP']) for r in rows])
        band(ax, x, y, COL[m], dots=False)  # mean line and 95% ribbon, as p1
    ax.set_xscale('log')
    ax.set_ylabel('genes per sample\n(ribosomal-protein genes excluded)')
    ax.set_xlabel('filtered reads per sample (subsampled)')
    ax.set_xticks(DEPTHS)
    ax.set_xticklabels(['50k', '100k', '250k', '500k', '1M', '2M'])
    ax.set_ylim(bottom=0)
    nlegend(ax, arms, 'lower right')
    ns = ' / '.join(str(len({r['well'] for r in rare(arms[m])})) for m in METH)
    note(ax, SETLBL[setname], 'every sample subsampled to each depth', f'n = {ns} samples', loc='br')
    save(fig, f'p2_genes_vs_depth_{setname}')
    # UCI at matched depth strip
    for d in (250000, 500000):
        data = {m: [float(r['B_corr']) for r in rare(arms[m]) if int(r['depth']) == d] for m in METH}
        fig, ax = plt.subplots(figsize=(2.2, 2.2))
        strip(ax, data, f'molecules (UMI) at {d//1000}k filtered reads', xlabel=SETLBL[setname])
        ax.set_ylim(bottom=0)
        note(ax, f'each sample subsampled to {d//1000}k', 'filtered reads; (n) = samples', loc='bl')
        save(fig, f'p1_molecules_at_{d//1000}k_{setname}')
# ---------- 3. read composition (native arms) ----------
CLS = [
    ('E', 'exonic, protein-coding', '#0B6E63'),
    ('X', 'exonic, other biotype', '#7FB8B0'),
    ('I', 'intronic', '#B9C7C3'),
    ('G', 'intergenic', '#DDE3E1'),
    ('R', 'rRNA', '#A5342B'),
    ('T', 'mitochondrial', '#E0A59F'),
]
fig, ax = plt.subplots(figsize=(3.0, 1.9))
comp = {}
for m in METH:
    d = {}
    cf = f'{U}/{SETS[COMPSET][m]}/composition_benchmark_wells.txt'
    cf = (
        cf if os.path.exists(cf) else f'{U}/{SETS[COMPSET][m]}/composition_species.txt'
    )  # benchmark wells only (the plate-level file counts every aligned well)
    for ln in open(cf):
        p = ln.rstrip('\n').split('\t')
        d[p[0]] = int(p[1])
    tot = sum(d.values())
    comp[m] = {k: 100 * d.get(k, 0) / tot for k, _, _ in CLS}
left = np.zeros(len(METH))
for k, lab, c in CLS:
    v = np.array([comp[m][k] for m in METH])
    ax.barh(range(len(METH)), v, left=left, color=c, edgecolor='white', lw=0.4, label=lab, height=0.7)
    left += v
ax.set_yticks(range(len(METH)))
ax.set_yticklabels(METH)
ax.invert_yaxis()
ax.set_xlabel('% of mapped reads (multimappers included)')
ax.set_xlim(0, 100)
if LEGEND:
    ax.legend(
        frameon=False, ncol=2, loc='upper center', bbox_to_anchor=(0.5, -0.38), handlelength=1.0, columnspacing=1.0
    )
save(fig, f'p3_read_composition_readtable_{COMPSET}')
json.dump(comp, open(f'{OUT}/p3_read_composition_readtable_{COMPSET}.json', 'w'), indent=1)
# ---------- 3b. intronic reads per sample (the p3 intronic class per well, one dot per benchmark sample) ----------
pw = {}
for m in METH:
    d = collections.defaultdict(dict)
    for ln in open(
        f'{U}/{SETS[COMPSET][m]}/composition_per_well.tsv'
    ):  # composition_per_well.sh: benchmark wells minus units_excluded (BOBseq n = 18)
        w, reg, a, u = ln.rstrip('\n').split('\t')
        d[w][reg] = int(a)
    pw[m] = {w: 100 * r.get('I', 0) / sum(r.values()) for w, r in d.items()}
fig, ax = plt.subplots(figsize=(2.2, 2.2))
strip(
    ax,
    {m: list(pw[m].values()) for m in METH},
    'intronic reads (% of mapped reads)\nper sample',
    xlabel=SETLBL[COMPSET],
)
ax.set_ylim(bottom=0)
save(fig, f'p3b_intronic_per_sample_readtable_{COMPSET}')
json.dump(pw, open(f'{OUT}/p3b_intronic_per_sample_readtable_{COMPSET}.json', 'w'), indent=1)
# ---------- 4. coverage along the transcript (canonical transcripts, unique reads, BOBseq both mates) ----------
S = [l.rstrip('\n').split('\t') for l in open(f'{A}/samples{SFX}.tsv') if l.strip()]
NCOV = f"{sum(1 for s_ in S if s_[2] == 'BOBseq')} BOBseq, 3 per competitor"
POOL = {}
for lab, bam, m in S:
    z = np.load(f'{A}/positions_canonical{SFX}/{lab}.npz')
    for k in ('gene', 't', 'L'):
        POOL.setdefault(m, {}).setdefault(k, []).append(z[k])
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
rng = np.random.default_rng(7)
# 4a breadth vs reads per gene + 4b 5' reached (from the A json if present, else recompute)
Aj = json.load(open(f'{CV}/A_breadth_vs_depth.json'))
depths = sorted(int(d) for d in Aj['DRUG-seq'])
for key, ylab, name in (
    ('occ_mean', 'transcript bins covered (of 20)', 'p4a_breadth_vs_reads_per_gene'),
    ('hit5', "genes with 5' end covered (%)", 'p4b_5prime_reached'),
):
    fig, ax = plt.subplots(figsize=SZ)
    for m in METH:
        ax.plot(
            depths,
            [Aj[m][str(d)][key] * (100 if key == 'hit5' else 1) for d in depths],
            'o-',
            color=COL[m],
            lw=1.6,
            ms=3,
        )
    ax.set_xscale('log')
    ax.set_xticks(depths)
    ax.set_xticklabels([str(d) for d in depths])
    ax.set_xlabel('unique reads per gene (subsampled)')
    ax.set_ylabel(ylab)
    if key == 'occ_mean':
        ax.set_ylim(0, 20)
        ax.axhline(20, color='#999', lw=0.5, ls=(0, (2, 2)))
    else:
        ax.set_ylim(0, 100)
    note(ax, SETLBL[COMPSET], 'each gene subsampled to x reads', NCOV, loc='br')
    legend(ax, 'lower right')
    save(fig, name)
# 4c cumulative distance from the poly(A) site (all canonical transcripts >= 500 nt, genes >= 100 reads, equal gene weight)
xs = np.arange(0, 3001, 25)
fig, ax = plt.subplots(figsize=SZ)
cdfs = {}
for m in METH:
    d = POOL[m]
    gl = [g for g, ii in IDX[m].items() if len(ii) >= 100 and d['L'][ii[0]] >= 500]
    cdf = np.zeros(len(xs))
    for g in gl:
        ii = IDX[m][g]
        dist = d['L'][ii] - 1 - d['t'][ii]
        cdf += np.searchsorted(np.sort(dist), xs, side='right') / len(ii)
    cdf /= len(gl)
    cdfs[m] = cdf
    ax.plot(xs, 100 * cdf, color=COL[m], lw=1.6)
ax.set_xlabel("distance from the poly(A) site (nt)")
ax.set_ylabel('read starts within distance (%)')
ax.set_ylim(0, 100)
ax.set_xlim(0, 3000)
note(ax, SETLBL[COMPSET], 'genes >= 100 reads, equal gene weight', NCOV, loc='br')
legend(ax, 'lower right')
save(fig, 'p4c_distance_from_polyA')
# 4d scaled profile at 200 reads per gene (equal gene weight), 20 bins
d0 = 200
common = [g for g in IDX[METH[0]] if all(len(IDX[m].get(g, [])) >= d0 for m in METH)]
fig, ax = plt.subplots(figsize=SZ)
x = np.arange(20) * 5 + 2.5
for m in METH:
    Pm = np.zeros(20)
    for g in common:
        ii = rng.choice(IDX[m][g], d0, replace=False)
        b = np.minimum((POOL[m]['t'][ii] / np.maximum(POOL[m]['L'][ii], 1) * 20).astype(int), 19)
        Pm += np.bincount(b, minlength=20) / d0
    ax.plot(x, 100 * Pm / len(common), color=COL[m], lw=1.6)
ax.axhline(5, color='#999', lw=0.5, ls=(0, (2, 2)))
ax.set_xlabel("position along the mature transcript (%, 5' to 3')")
ax.set_ylabel('read starts per 5% bin (%)')
ax.set_xlim(0, 100)
ax.set_ylim(bottom=0)
note(ax, SETLBL[COMPSET], f'{len(common):,} genes at {d0} unique reads each', NCOV, loc='tc')
legend(ax)
save(fig, 'p4d_scaled_profile_200_reads_per_gene')
json.dump(
    {
        'genes_scaled_profile': len(common),
        'cdf_polyA': {m: {str(int(c)): float(np.interp(c, xs, cdfs[m])) for c in (300, 500, 1000)} for m in METH},
    },
    open(f'{OUT}/p4_numbers.json', 'w'),
    indent=1,
)
# ---------- 4e-4g. base coverage (the browser view): covered fraction vs reads per gene, peakiness, heatmaps ----------
V2j = json.load(open(f'{CV}/V2_fraction_covered_vs_depth.json'))
vd = sorted(int(d) for d in V2j['BOBseq'])
fig, ax = plt.subplots(figsize=SZ)
for m in METH:
    ax.plot(vd, [100 * V2j[m][str(d)]['frac1'] for d in vd], 'o-', color=COL[m], lw=1.6, ms=3)
ax.set_xscale('log')
ax.set_xticks(vd)
ax.set_xticklabels([str(d) for d in vd])
ax.set_xlabel('unique reads per gene (subsampled)')
ax.set_ylabel('transcript covered by >= 1 read (%)')
ax.set_ylim(0, 100)
note(ax, SETLBL[COMPSET], 'each gene subsampled to x reads', NCOV, loc='br')
legend(ax)
save(fig, 'p4e_transcript_covered_vs_reads_per_gene')
fig, ax = plt.subplots(figsize=SZ)
for m in METH:
    ax.plot(vd, [100 * V2j[m][str(d)]['peak'] for d in vd], 'o-', color=COL[m], lw=1.6, ms=3)
ax.set_xscale('log')
ax.set_xticks(vd)
ax.set_xticklabels([str(d) for d in vd])
ax.set_xlabel('unique reads per gene (subsampled)')
ax.set_ylabel("reads in the gene's densest 300 nt (%)")
ax.set_ylim(0, 100)
note(ax, SETLBL[COMPSET], 'each gene subsampled to x reads', NCOV, loc='bl')
legend(ax)
save(fig, 'p4f_peakiness_vs_reads_per_gene')
# per-gene strips at 200 reads (covered fraction, densest-300 share) and the heatmaps, recomputed from the positions
alen = {}
for lab, bam, m in S:
    z = np.load(f'{A}/positions_canonical{SFX}/{lab}.npz')
    alen.setdefault(m, []).append(z['alen'])
for m in METH:
    POOL[m]['alen'] = np.concatenate(alen[m])
tx_len = np.load(f'{A}/positions_canonical{SFX}/{S[0][0]}.npz')['tx_len']


def gene_cov(m, ii, L):
    t = POOL[m]['t'][ii].astype(int)
    e = np.minimum(t + POOL[m]['alen'][ii].astype(int), L)
    diff = np.zeros(L + 1)
    np.add.at(diff, t, 1)
    np.add.at(diff, e, -1)
    return np.cumsum(diff)[:L], t


d0 = 200
common = sorted([g for g in IDX[METH[0]] if all(len(IDX[m].get(g, [])) >= d0 for m in METH)], key=lambda g: tx_len[g])
rng2 = np.random.default_rng(11)
cov_frac = {m: [] for m in METH}
peak = {m: [] for m in METH}
nb = 50
Hm = {m: np.zeros((len(common), nb)) for m in METH}
for m in METH:
    for r, g in enumerate(common):
        ii = rng2.choice(IDX[m][g], d0, replace=False)
        L = int(tx_len[g])
        cov, t = gene_cov(m, ii, L)
        cov_frac[m].append(100 * np.mean(cov > 0))
        cs = np.cumsum(np.bincount(t, minlength=L))
        w = min(300, L)
        peak[m].append(100 * (cs[w - 1 :] - np.r_[0, cs[:-w]]).max() / d0)
        idx = np.minimum((np.arange(L) * nb) // L, nb - 1)
        b = np.bincount(idx, weights=cov, minlength=nb) / np.bincount(idx, minlength=nb)
        Hm[m][r] = b / b.max() if b.max() > 0 else b
fig, ax = plt.subplots(figsize=(2.2, 2.2))
strip(ax, cov_frac, f'transcript covered by >= 1 read (%)\nat {d0} reads per gene', xlabel=SETLBL[COMPSET])
ax.set_ylim(0, 100)
save(fig, 'p4g_transcript_covered_at_200_reads')
fig, ax = plt.subplots(figsize=(2.2, 2.2))
strip(ax, peak, f"reads in the gene's densest 300 nt (%)\nat {d0} reads per gene", xlabel=SETLBL[COMPSET])
ax.set_ylim(0, 100)
save(fig, 'p4h_peakiness_at_200_reads')
fig, ax = plt.subplots(1, len(METH), figsize=(1.3 * len(METH), 2.4), sharey=True)  # one panel per method
for a, m in zip(ax, METH):
    a.imshow(Hm[m], aspect='auto', cmap='Greys', vmin=0, vmax=1, interpolation='nearest')
    a.set_title(m, fontsize=7, color=COL[m])
    a.set_xticks([0, nb - 1])
    a.set_xticklabels(["5'", "3'"])
    a.tick_params(length=2)
    for sp in a.spines.values():
        sp.set_visible(True)
        sp.set_linewidth(0.5)
ax[0].set_ylabel(f'{len(common):,} genes (short to long)\n{d0} unique reads per gene')
ax[0].set_yticks([])
fig.text(
    0.5,
    -0.04,
    'coverage along the mature transcript (each row scaled to its own peak)',
    ha='center',
    va='top',
    fontsize=7,
)
save(fig, 'p4i_coverage_heatmaps_200_reads')
json.dump(
    {
        'genes': len(common),
        'covered_frac_median': {m: float(np.median(cov_frac[m])) for m in METH},
        'peak_median': {m: float(np.median(peak[m])) for m in METH},
    },
    open(f'{OUT}/p4_base_coverage_numbers.json', 'w'),
    indent=1,
)
# ---------- 5. junctions per sample at matched depth ----------
J = []
meth_of = {lab: m for lab, _, m in S}
for f in glob.glob(f'{A}/junctions{SFX}/*.tsv'):
    J += [
        r for r in csv.DictReader(open(f), delimiter='\t') if r['sample'] in meth_of
    ]  # only the samples of the set (stale files of dropped samples are ignored)
if J:
    for d in (250000, 500000):
        data = {
            m: [
                int(r['junctions_ge2'])
                for r in J
                if meth_of[r['sample']] == m and int(r['depth']) == d and r['at_native'] == '0'
            ]
            for m in METH
        }
        fig, ax = plt.subplots(figsize=(2.2, 2.2))
        strip(
            ax,
            data,
            f'unique splice junctions per sample\n(>= 2 reads, at {d//1000}k unique reads)',
            xlabel=SETLBL[COMPSET],
        )
        ax.set_ylim(bottom=0)
        note(ax, f'subsampled to {d//1000}k', 'unique reads per sample', '(n) = samples', loc='tl')
        save(fig, f'p5_junctions_at_{d//1000}k')
        data = {
            m: [
                100 * int(r['spliced_reads']) / int(r['reads_used'])
                for r in J
                if meth_of[r['sample']] == m and int(r['depth']) == d and r['at_native'] == '0'
            ]
            for m in METH
        }
        fig, ax = plt.subplots(figsize=(2.2, 2.2))
        strip(ax, data, f'spliced reads (%) at {d//1000}k unique reads', xlabel=SETLBL[COMPSET])
        ax.set_ylim(bottom=0)
        note(ax, f'subsampled to {d//1000}k', 'unique reads per sample', '(n) = samples', loc='bl')
        save(fig, f'p5_spliced_fraction_at_{d//1000}k')
    fig, ax = plt.subplots(figsize=SZ)
    for m in METH:
        ds = sorted({int(r['depth']) for r in J if meth_of[r['sample']] == m and r['at_native'] == '0'})
        med = [
            np.median(
                [
                    int(r['junctions_ge2'])
                    for r in J
                    if meth_of[r['sample']] == m and int(r['depth']) == d and r['at_native'] == '0'
                ]
            )
            for d in ds
        ]
        pts = [
            (int(r['depth']), int(r['junctions_ge2'])) for r in J if meth_of[r['sample']] == m and r['at_native'] == '0'
        ]
        ax.scatter(
            [d * np.exp(np.random.default_rng(0).normal(0, 0.04)) for d, _ in pts],
            [v for _, v in pts],
            s=5,
            color=COL[m],
            alpha=0.5,
            lw=0,
            zorder=2,
        )
        ax.plot(ds, med, 'o-', color=COL[m], lw=1.6, ms=3, zorder=3)
    ax.set_xscale('log')
    ax.set_xlabel('unique reads per sample')
    ax.set_ylabel('unique splice junctions per sample\n(>= 2 reads per junction)')
    ax.set_xticks([100000, 250000, 500000, 1000000, 2000000])
    ax.set_xticklabels(['100k', '250k', '500k', '1M', '2M'])
    ax.set_ylim(bottom=0)
    ax.set_xlabel('unique reads per sample (subsampled)')
    ns = ' / '.join(str(len({r['sample'] for r in J if meth_of[r['sample']] == m})) for m in METH)
    note(ax, SETLBL[COMPSET], 'every sample subsampled to each depth', f'n = {ns} samples', loc='tl')
    legend(ax)
    save(fig, 'p5_junctions_vs_depth')
    for thr in (3, 5):
        if f'junctions_ge{thr}' not in J[0]:
            continue
        fig, ax = plt.subplots(figsize=SZ)
        for m in METH:
            ds = sorted({int(r['depth']) for r in J if meth_of[r['sample']] == m and r['at_native'] == '0'})
            med = [
                np.median(
                    [
                        int(r[f'junctions_ge{thr}'])
                        for r in J
                        if meth_of[r['sample']] == m and int(r['depth']) == d and r['at_native'] == '0'
                    ]
                )
                for d in ds
            ]
            pts = [
                (int(r['depth']), int(r[f'junctions_ge{thr}']))
                for r in J
                if meth_of[r['sample']] == m and r['at_native'] == '0'
            ]
            ax.scatter(
                [d * np.exp(np.random.default_rng(0).normal(0, 0.04)) for d, _ in pts],
                [v for _, v in pts],
                s=5,
                color=COL[m],
                alpha=0.5,
                lw=0,
                zorder=2,
            )
            ax.plot(ds, med, 'o-', color=COL[m], lw=1.6, ms=3, zorder=3)
        ax.set_xscale('log')
        ax.set_xlabel('unique reads per sample (subsampled)')
        ax.set_ylabel(f'unique splice junctions per sample\n(>= {thr} reads per junction)')
        ax.set_xticks([100000, 250000, 500000, 1000000, 2000000])
        ax.set_xticklabels(['100k', '250k', '500k', '1M', '2M'])
        ax.set_ylim(bottom=0)
        legend(ax)
        save(fig, f'p5_junctions_vs_depth_ge{thr}')
        data = {
            m: [
                int(r[f'junctions_ge{thr}'])
                for r in J
                if meth_of[r['sample']] == m and int(r['depth']) == 500000 and r['at_native'] == '0'
            ]
            for m in METH
        }
        fig, ax = plt.subplots(figsize=(2.2, 2.2))
        strip(
            ax,
            data,
            f'unique splice junctions per sample\n(>= {thr} reads, at 500k unique reads)',
            xlabel=SETLBL[COMPSET],
        )
        ax.set_ylim(bottom=0)
        save(fig, f'p5_junctions_at_500k_ge{thr}')
    with open(f'{OUT}/p5_junctions_per_sample.tsv', 'w') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(
            ['sample', 'method', 'depth', 'reads_used', 'spliced_reads', 'junctions_ge1', 'junctions_ge2', 'at_native']
        )
        w.writerows(
            [
                [
                    r['sample'],
                    meth_of[r['sample']],
                    r['depth'],
                    r['reads_used'],
                    r['spliced_reads'],
                    r['junctions_ge1'],
                    r['junctions_ge2'],
                    r['at_native'],
                ]
                for r in J
            ]
        )
print('done')

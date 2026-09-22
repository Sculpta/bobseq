#!/usr/bin/env python3
"""Per-sample supplement on all 21 human BOBseq wells (replicate variation; the three risdiplam 500 nM wells, input failures excluded
from every main figure, are shown here marked). BM_COVSET=native|50nt -> <fig set>/supplement_per_sample_21/ (separate from the figure
panels). Same sources and rules as the main panels for the 18 benchmark wells; the Ris 500 nM wells come from supp21_* (arm path).
Panels (one SVG + PNG each, values/ TSV beside): S01 read fate, S02 filtered reads and molecules, S03 duplicate rate, S04 molecules vs depth,
S05 genes vs depth, S06 composition, S07 intronic %, S08 read length (per mate), S09 insert length (native), S10 aligned bases, S11 Picard
profile, S12 Picard balance, (barcode accuracy per well only in the values table)."""

from settings import WORK, PROJECT, COVERAGE, RUN_JSON
import os, re, sys, csv, json, glob, collections, numpy as np, matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.ticker
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, H)
from palette import COL, CLS

P = WORK
U = f'{P}/benchmark_uniform'
PR = PROJECT
A = COVERAGE
SET = os.environ.get('BM_COVSET', 'native')
assert SET == 'native', 'the per-sample supplement is native only'
TAG = 'native'
FIG = 'preprint_figures_native' if SET == 'native' else f'preprint_figures_{SET}'
ARM = f'{U}/bob57_24plex_pe_native' if SET == 'native' else f'{U}/bob57_24plex_pe_hs'
OUT = f'{P}/{FIG}/supplement_per_sample_21'
V = f'{OUT}/values'
os.makedirs(V, exist_ok=True)
MAIN = f'{P}/{FIG}/main_figure_panels'
sfx = '' if SET == 'native' else '_50nt'
MET = f'{P}/{FIG}/coverage_architecture/picard_metrics/{TAG}'
RUN = json.load(open(RUN_JSON))
LAB = {
    w: l.replace(' ', '_').replace('/', '_')
    for w, l in RUN['bobcode_labels'].items()
    if RUN['tso_species_map'][w] == 'human'
}
# the 21 human wells of the run in configuration order (18 benchmark wells and the three risdiplam 500 nM input failures);
# the arm's wells_keep.txt lists only the 18 benchmark wells
WELLS = list(LAB)
assert len(WELLS) == 21, len(WELLS)


def group(n):
    return re.sub(r'[-_]\d+$', '', n).replace('_', ' ')


GROUPS = [
    'Hek control',
    'CADM1 1nM',
    'CADM1 100nM',
    'CHX dose 1ug mL',
    'CHX dose 50ug mL',
    'Ris dose 25mM',
    'Ris dose 500mM',
]
GCOL = {
    'Hek control': '#2E6B32',
    'CADM1 1nM': '#4C78A8',
    'CADM1 100nM': '#1F4E79',
    'CHX dose 1ug mL': '#E6B422',
    'CHX dose 50ug mL': '#B8860B',
    'Ris dose 25mM': '#B5379B',
    'Ris dose 500mM': '#8C8C8C',
}


def gkey(w):
    g = group(LAB[w])
    return GROUPS.index(g) if g in GROUPS else 99


WELLS = sorted(WELLS, key=lambda w: (gkey(w), LAB[w]))
NAMES = [LAB[w] for w in WELLS]
RIS500 = {w for w in WELLS if group(LAB[w]) == 'Ris dose 500mM'}
EXC = {l.split('\t')[0] for i, l in enumerate(open(f'{ARM}/units_excluded.tsv')) if i > 0}


def col(w):
    return GCOL.get(group(LAB[w]), '#444444')


def short(n):
    return (
        n.replace('_dose', '')
        .replace('_mL', '')
        .replace('Hek_control', 'Hek')
        .replace('CHX_', 'CHX ')
        .replace('CADM1_', 'CADM1 ')
        .replace('Ris_', 'Ris ')
    )


plt.rcParams.update(
    {
        'svg.fonttype': 'none',
        'font.family': ['Arial', 'Nimbus Sans', 'DejaVu Sans'],
        'font.size': 7,
        'axes.labelsize': 7,
        'xtick.labelsize': 6,
        'ytick.labelsize': 6.5,
        'legend.fontsize': 6,
        'axes.linewidth': 0.6,
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


def save(fig, name):
    os.makedirs(os.path.dirname(f'{OUT}/{name}'), exist_ok=True)   # superseded/ variants live in a subdirectory
    fig.tight_layout(pad=0.4)
    fig.savefig(f'{OUT}/{name}.svg', bbox_inches='tight', pad_inches=0.03)
    fig.savefig(f'{OUT}/{name}.png', dpi=200, bbox_inches='tight', pad_inches=0.03)
    plt.close(fig)
    s = open(f'{OUT}/{name}.svg').read()
    s = re.sub(r"font-family:\s*[^;\"]+", 'font-family: Arial', s)
    open(f'{OUT}/{name}.svg', 'w').write(s)
    print('wrote', name, flush=True)


def glegend(ax, loc='upper right', ncol=1):
    ax.legend(
        handles=[Patch(color=GCOL[g]) for g in GROUPS],
        labels=[g + (' (input failure)' if g == 'Ris dose 500mM' else '') for g in GROUPS],
        loc=loc,
        ncol=ncol,
        frameon=False,
        fontsize=5.5,
        handlelength=1.0,
        borderaxespad=0.3,
        labelspacing=0.25,
    )


def bars(ax, vals, ylabel, log=False, fmt=None):
    x = np.arange(len(WELLS))
    ax.bar(x, [vals[w] for w in WELLS], color=[col(w) for w in WELLS], width=0.75, hatch=None)
    [ax.patches[i].set_hatch('///') for i, w in enumerate(WELLS) if w in RIS500]
    ax.set_xticks(x)
    ax.set_xticklabels([short(LAB[w]) for w in WELLS], rotation=90, fontsize=5.5)
    ax.set_ylabel(ylabel)
    ax.set_xlim(-0.7, len(WELLS) - 0.3)
    if log:
        ax.set_yscale('log')


def two(kind, f18, fsup):
    """rows of a per-well table from the benchmark file plus the supplement file"""
    rows = []
    for fn in (f18, fsup):
        if os.path.exists(fn):
            rows += (
                list(csv.DictReader(open(fn), delimiter='\t'))
                if kind == 'dict'
                else [l.rstrip('\n').split('\t') for l in open(fn) if l.strip()]
            )
    return rows


# ---- per-well core numbers ----
inp = {
    l.split('\t')[0]: int(l.split('\t')[1]) for l in open(f'{ARM}/input_counts.tsv') if l.strip()
}  # current table (18 benchmark wells after the N-UMI drop; the native file also lists the other wells before that drop)
[
    inp.setdefault(l.split('\t')[0], int(l.split('\t')[1]))
    for l in open(f'{U}/bob57_24plex_pe_native/input_counts.tsv')
    if l.strip()
]  # 50-nt set: reads into STAR are the same wells' reads (truncated read 2)
for l in open(
    f'{ARM}/supp21/n_umi_dropped_per_well.tsv'
):  # Ris 500 nM wells: apply the same N-UMI drop the 18 wells carry in input_counts.tsv
    w, d = l.rstrip('\n').split('\t')
    inp[w] = inp[w] - int(d)
pm = {
    r[0]: (int(r[1]), int(r[2]))
    for r in two('list', f'{ARM}/perwell_mapped.tsv', f'{ARM}/supp21/perwell_mapped.tsv')
    if len(r) >= 3
}  # the 18-well table (read_fate_panels.py uses the same)
rpw = {
    r[0]: int(r[1])
    for r in two(
        'list', f'{ARM}/stats_E_U/basic/reads_per_well.tsv', f'{ARM}/stats_E_U_supp21/basic/reads_per_well.tsv'
    )
}
pb = {
    r['well']: r
    for r in two('dict', f'{ARM}/stats_E_U/basic/perwell_B.tsv', f'{ARM}/stats_E_U_supp21/basic/perwell_B.tsv')
}
rar = two('dict', f'{ARM}/stats_E_U/rarefied.tsv', f'{ARM}/stats_E_U_supp21/rarefied.tsv')
missing = [LAB[w] for w in WELLS if w not in pm or w not in rpw or w not in pb]
assert not missing, f'wells without tables: {missing}'
core = {
    w: dict(
        sample=LAB[w],
        group=group(LAB[w]),
        input_failure=w in RIS500,
        reads_in=inp[w],
        mapped=pm[w][0],
        unique=pm[w][1],
        filtered=rpw[w],
        molecules_UCI=float(pb[w]['B_corr']),
        genes_at_native=int(pb[w]['genes_B']),
        duplicate_rate_pct=100 * (1 - float(pb[w]['B_corr']) / rpw[w]),
    )
    for w in WELLS
}
# S01 read fate
fig, ax = plt.subplots(figsize=(2.8, 2.3))
steps = ['reads_in', 'mapped', 'unique', 'filtered']
for w in WELLS:
    y = [100 * core[w][k] / core[w]['reads_in'] for k in steps]
    ax.plot(range(4), y, color=col(w), lw=0.9, ls=(0, (3, 2)) if w in RIS500 else '-', alpha=0.9)
ax.set_xticks(range(4))
ax.set_xticklabels(['reads\ninto STAR', 'mapped', 'uniquely\nmapped', 'filtered'], fontsize=6)
ax.set_ylim(0, 100)
ax.set_ylabel('% of reads into STAR (per well)')
glegend(ax, 'lower left')
save(fig, f'S01_read_fate_{TAG}')
# S02 filtered reads and molecules, S03 duplicate rate
fig, ax = plt.subplots(figsize=(3.4, 2.4))
bars(
    ax, {w: core[w]['molecules_UCI'] for w in WELLS}, 'filtered reads per sample\nafter de-duplication (molecules, UMI)'
)
ax.set_yscale('log')
save(fig, f'S02_reads_after_dedup_{TAG}')
fig, ax = plt.subplots(figsize=(3.4, 2.4))
bars(ax, {w: core[w]['duplicate_rate_pct'] for w in WELLS}, 'duplicate rate (%)\nas sequenced')
ax.set_ylim(0, 100)
save(fig, f'S03_duplicate_rate_{TAG}')
# S04 / S05 depth curves per well
DEP = [50000, 100000, 250000, 500000, 1000000, 2000000]
for name, key, yl in (
    ('S04_molecules_vs_depth', 'B_corr', 'molecules per well (UMI)'),
    ('S05_genes_vs_depth', 'genes_B_noRP', 'genes per well\n(ribosomal-protein genes excluded)'),
):
    fig, ax = plt.subplots(figsize=(2.8, 2.3))
    for w in WELLS:
        pts = sorted(
            (int(r['depth']), float(r[key]))
            for r in rar
            if r['well'] == w and r.get('B_corr') and r['at_native'] == '0'
        )
        if pts:
            ax.plot(
                [p[0] for p in pts],
                [p[1] for p in pts],
                color=col(w),
                lw=0.9,
                ls=(0, (3, 2)) if w in RIS500 else '-',
                marker='o',
                ms=1.8,
            )
    ax.set_xscale('log')
    ax.set_xticks(DEP)
    ax.set_xticklabels(['50k', '100k', '250k', '500k', '1M', '2M'])
    ax.set_xlabel('filtered reads per well (subsampled)')
    ax.set_ylabel(yl)
    ax.set_ylim(0, None)
    ax.yaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(
            lambda v, _: f'{v/1e6:g}M' if v >= 1e6 else (f'{v/1e3:g}k' if v >= 1e3 else f'{v:g}')
        )
    )
    glegend(ax, 'upper left' if key == 'B_corr' else 'lower right')
    save(fig, f'{name}_{TAG}')  # S05 legend bottom right, away from the curves
# S03b duplicate rate at 250k filtered reads (matched depth), S04b molecules vs depth on a log axis
d250 = {
    r['well']: 100 * (1 - float(r['B_corr']) / 250000)
    for r in rar
    if r['well'] in LAB and r.get('B_corr') and r['at_native'] == '0' and int(r['depth']) == 250000
}
fig, ax = plt.subplots(figsize=(3.4, 2.4))
bars(ax, {w: d250.get(w, np.nan) for w in WELLS}, 'duplicate rate (%)\nat 250k filtered reads')
ax.set_ylim(0, 100)
save(fig, f'S03b_duplicate_rate_at_250k_{TAG}')
fig, ax = plt.subplots(figsize=(2.8, 2.3))
for w in WELLS:
    pts = sorted(
        (int(r['depth']), float(r['B_corr']))
        for r in rar
        if r['well'] == w and r.get('B_corr') and r['at_native'] == '0'
    )
    if pts:
        ax.plot(
            [q[0] for q in pts],
            [q[1] for q in pts],
            color=col(w),
            lw=0.9,
            ls=(0, (3, 2)) if w in RIS500 else '-',
            marker='o',
            ms=1.8,
        )
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xticks(DEP)
ax.set_xticklabels(['50k', '100k', '250k', '500k', '1M', '2M'])
ax.set_xlabel('filtered reads per well (subsampled)')
ax.set_ylabel('molecules per well (UMI)')
glegend(ax, 'upper left')
save(fig, f'S04b_molecules_vs_depth_log_{TAG}')
# S06 composition, S07 intronic
comp = {
    r['sample']: r
    for r in csv.DictReader(open(f'{MAIN}/values/p3_composition_rules_{TAG}_per_sample.tsv'), delimiter='\t')
    if r['method'] == 'BOBseq'
}
sup = f'{PR}/results/supplement_21/composition_rules_ris500_{SET}.json'
if os.path.exists(sup):
    for n, c in json.load(open(sup)).items():
        comp[n] = {**{k: str(c.get(k, 0)) for k, _, _ in CLS}, 'mapped_records': str(sum(c.values()))}
CLASSES = [(k, l, c) for k, l, c in CLS]
fig, ax = plt.subplots(figsize=(3.4, 2.4))
bottom = np.zeros(len(WELLS))
x = np.arange(len(WELLS))
compv = {}
for k, l, c in CLASSES:
    v = np.array(
        [
            (
                100 * float(comp[LAB[w]][k]) / sum(float(comp[LAB[w]][kk]) for kk, _, _ in CLASSES)
                if LAB[w] in comp
                else np.nan
            )
            for w in WELLS
        ]
    )
    compv[k] = v
    ax.bar(x, v, 0.75, bottom=bottom, color=c, lw=0)
    bottom += np.nan_to_num(v)
ax.set_xticks(x)
ax.set_xticklabels([short(LAB[w]) for w in WELLS], rotation=90, fontsize=5.5)
ax.set_ylabel('% of mapped reads')
ax.set_ylim(0, 100)
ax.set_xlim(-0.7, len(WELLS) - 0.3)
save(fig, f'S06_composition_{TAG}')
fig, ax = plt.subplots(figsize=(3.4, 2.4))
bars(ax, {w: compv['intronic'][i] for i, w in enumerate(WELLS)}, 'intronic reads\n(% of mapped)')
save(fig, f'S07_intronic_{TAG}')
# S08 read length, S10 aligned bases (positions npz), S09 insert (native)
S = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open(f'{A}/samples{sfx}.tsv') if l.strip()}
lab_of = {os.path.basename(r[1])[:-4]: r[0] for r in S.values() if r[2] == 'BOBseq'}
for k, n in ((1, 'Ris_dose_500mM-1'), (2, 'Ris_dose_500mM-2'), (3, 'Ris_dose_500mM-3')):
    lab_of[n] = f'bobseq-ris-500mM-{k}'
rl = {}
ab = {}
for w in WELLS:
    n = LAB[w]
    f = f'{A}/positions_canonical{sfx}/{lab_of.get(n, "")}.npz'
    if not os.path.exists(f):
        continue
    z = np.load(f)
    al = z['alen'].astype(np.int64)
    mate = z['mate']
    nf = int((mate == 1).sum()) if (mate > 0).any() else len(al)
    rl[w] = {
        'm1': float(np.median(al[mate == 1])) if (mate == 1).any() else float(np.median(al)),
        'm2': float(np.median(al[mate == 2])) if (mate == 2).any() else None,
    }
    ab[w] = {'total': int(al.sum()), 'per_fragment': al.sum() / max(nf, 1), 'reads': len(al)}
if rl:
    fig, ax = plt.subplots(figsize=(3.4, 2.4))
    x = np.arange(len(WELLS))
    ax.bar(
        x - 0.2,
        [rl.get(w, {}).get('m1', np.nan) for w in WELLS],
        0.4,
        color=[col(w) for w in WELLS],
        label='read 2 (mate 1)' if SET == 'native' else '50-nt read',
    )
    if SET == 'native':
        ax.bar(
            x + 0.2,
            [rl.get(w, {}).get('m2') or np.nan for w in WELLS],
            0.4,
            color=[col(w) for w in WELLS],
            alpha=0.45,
            label='read 1 (mate 2)',
        )
    ax.set_xticks(x)
    ax.set_xticklabels([short(LAB[w]) for w in WELLS], rotation=90, fontsize=5.5)
    ax.set_ylabel('median aligned length (nt)\ndark read 2, light read 1')
    ax.set_xlim(-0.7, len(WELLS) - 0.3)
    save(fig, f'S08_read_length_{TAG}')
    fig, ax = plt.subplots(figsize=(3.4, 2.4))
    bars(
        ax,
        {w: ab[w]['total'] / 1e9 if w in ab else np.nan for w in WELLS},
        'aligned bases (Gb)\nunique reads, as sequenced',
        log=True,
    )
    save(fig, f'S10_aligned_bases_{TAG}')
ins = {}
if SET == 'native':
    for w in WELLS:
        f = f'{PR}/results/library_metrics/insert_per_well/{LAB[w]}.json'
        if os.path.exists(f):
            d = json.load(open(f))
            ins[w] = d
    if ins:
        fig, ax = plt.subplots(figsize=(2.8, 2.3))
        xs = np.arange(1, 601)
        for w in WELLS:
            if w not in ins:
                continue
            h = np.array(ins[w]['hist_1nt_to_3000'], float)
            pct = 100 * h / h.sum()
            y = np.convolve(pct, np.ones(5) / 5, 'same')[
                1:601
            ]  # % of the well's pairs per nt, 5-nt running mean, x to 600
            ax.plot(xs, y, color=col(w), lw=0.9, ls=(0, (3, 2)) if w in RIS500 else '-', alpha=0.9)
        ax.set_xlabel('insert length on the mature transcript (nt)')
        ax.set_ylabel('% of deduplicated pairs per nt')
        ax.set_xlim(0, 600)
        ax.set_ylim(0, None)
        glegend(ax, 'upper right')
        save(fig, f'S09_insert_length_{TAG}')
# S14 poly(A)-anchored coverage profile per well (the p4m2 gene set: >= 1 kb, ribosomal-protein genes excluded; per well all of its unique
# reads on those genes, genes with >= 20 reads in that well, coverage / gene mean anchored at the 3' end, 3-kb window, mean over genes, 25-nt smoothing)
thr = json.load(open(f'{MAIN}/p4_thresholds_numbers.json'))
GS = set(thr['gene_names'])
W3 = 3000
MINR = 20


def sm(y, k=25):
    y = np.asarray(y, float)
    v = np.where(np.isnan(y), 0, y)
    n = (~np.isnan(y)).astype(float)
    ker = np.ones(k)
    return np.convolve(v, ker, 'same') / np.maximum(np.convolve(n, ker, 'same'), 1)


prof = {}
ngen = {}
for w in WELLS:
    f = f'{A}/positions_canonical{sfx}/{lab_of.get(LAB[w], "")}.npz'
    if not os.path.exists(f):
        continue
    # arrays loaded once per well: indexing the npz handle inside the loop would decompress the whole array every time
    z = np.load(f)
    genes_ = list(z['genes'])
    tl = z['tx_len']
    g = z['gene']
    T_ = z['t']
    AL_ = z['alen']
    o = np.argsort(g, kind='stable')
    gs = g[o]
    st = np.flatnonzero(np.r_[True, gs[1:] != gs[:-1]])
    en = np.r_[st[1:], len(gs)]
    idx = {int(gs[a]): o[a:b] for a, b in zip(st, en)}
    Uw = []
    for gi, ii in idx.items():
        if genes_[gi] not in GS or len(ii) < MINR:
            continue
        L = int(tl[gi])
        t = T_[ii].astype(int)
        e = np.minimum(t + AL_[ii].astype(int), L)
        d = np.zeros(L + 1)
        np.add.at(d, t, 1)
        np.add.at(d, e, -1)
        c = np.cumsum(d)[:L]
        zz = c / c.mean()
        row = np.full(W3, np.nan)
        wl = min(W3, L)
        row[:wl] = zz[::-1][:wl]
        Uw.append(row)
    if Uw:
        Uw = np.array(Uw)
        prof[w] = sm(np.nanmean(Uw, 0))
        ngen[w] = len(Uw)
print('S14: wells with a profile', len(prof), 'of', len(WELLS), '| gene set', len(GS), flush=True)
if prof:
    fig, ax = plt.subplots(figsize=(2.8, 2.3))
    xx = np.arange(W3)
    for w in WELLS:
        if w in prof:
            ax.plot(xx, prof[w], color=col(w), lw=0.9, ls=(0, (3, 2)) if w in RIS500 else '-')
    ax.axhline(1, color='#9a9a9a', lw=0.6, ls=(0, (2, 2)))
    ax.set_xlim(W3, 0)
    ax.set_ylim(0, 8.5)
    ax.set_xlabel('distance from the poly(A) site (nt)')
    ax.set_ylabel('coverage / gene mean\n(genes >= 1 kb, RP excluded)')
    glegend(ax, 'upper left')
    save(
        fig, f'superseded/S14_human_only_coverage_vs_nt_from_polyA_{TAG}'
    )  # y 0 to 8.5 = the range of the main p4m2 panel; the all-24 S14 of supplement_polyA_24.py replaces this human-only panel
    with open(f'{V}/polyA_profiles_{TAG}.tsv', 'w') as f:
        f.write('sample\tgenes\t' + '\t'.join(f'nt{d}' for d in range(W3)) + '\n')
        [f.write(f'{LAB[w]}\t{ngen[w]}\t' + '\t'.join(f'{v:.4f}' for v in prof[w]) + '\n') for w in WELLS if w in prof]


# S11 / S12 Picard per well
def parse(path):
    rows = [l.rstrip('\n').split('\t') for l in open(path)]
    hist = []
    inh = False
    for r in rows:
        if r and r[0].startswith('## HISTOGRAM'):
            inh = True
            continue
        if inh and len(r) >= 2 and r[0].isdigit():
            hist.append(float(r[1]))
    return np.array(hist)


pic = {}
p = np.arange(101) / 100
for w in WELLS:
    f = f'{MET}/{LAB[w]}.RNA_Metrics.txt'
    if os.path.exists(f) and '## HISTOGRAM' in open(f).read():
        h = parse(f)
        if len(h) == 101 and h.sum() > 0:
            pic[w] = {'profile': h, 'balance': float(2 * (p * h).sum() / h.sum())}
if pic:
    fig, ax = plt.subplots(figsize=(2.8, 2.3))
    for w in WELLS:
        if w in pic:
            ax.plot(np.arange(101), pic[w]['profile'], color=col(w), lw=0.9, ls=(0, (3, 2)) if w in RIS500 else '-')
    ax.axhline(1, color='#7f7f7f', lw=0.7, ls=(0, (4, 3)))
    ax.set_xlabel("gene-body percentile (5' to 3')")
    ax.set_ylabel('normalized read coverage')
    ax.set_ylim(0, 7)
    glegend(ax, 'upper left')
    save(fig, f'S11_picard_profile_{TAG}')  # y to 7 like the main Picard panel
    fig, ax = plt.subplots(figsize=(3.4, 2.4))
    bars(ax, {w: pic[w]['balance'] if w in pic else np.nan for w in WELLS}, "5' to 3' balance\n(2 x centroid)")
    ax.axhline(1, color='#7f7f7f', lw=0.7, ls=(0, (4, 3)))
    ax.set_ylim(0, 2)
    save(fig, f'S12_picard_balance_{TAG}')
# S13 barcode accuracy (final funnel on the native deduplicated pair BAMs, set-independent)
fs = json.load(open(f'{PR}/results/barcode_accuracy_funnel/funnel_summary.json'))['per_sample']
acc = {}
for w in WELLS:
    d = fs.get(LAB[w])
    if d:
        l = d['levels'][-1]
        acc[w] = (l['accuracy'], l['wilson_lo'], l['wilson_hi'], l['n'])
# (no S13 panel is drawn; the per-well accuracy goes to values/per_well_*.tsv)

# S15 unique splice junctions (>= 3 molecules) vs deduplicated molecules per well (the per-sample junction BEDs of extract_junctions.sh; same
# rarefaction as the main p5 panel: binomial thinning of each junction's molecule count at p = depth / the well's dedup MAPQ-255 molecules)
import gzip
from scipy import stats as _st2

JD = f'{PR}/data/junctions'
JMETA = {r['acc_number']: r for r in csv.DictReader(open(f'{JD}/sample_metadata.csv'))}
JDEP = [25000, 50000, 100000, 250000, 500000, 1000000, 2000000]
KJ = 3
junc = {}
for w in WELLS:
    jacc = f'bobseq_pe_native__{LAB[w]}'
    if jacc not in JMETA or not os.path.exists(f'{JD}/{jacc}_junctions.bed.gz'):
        continue
    c = collections.Counter()
    for l in gzip.open(f'{JD}/{jacc}_junctions.bed.gz', 'rt'):
        f_ = l.rstrip('\n').split('\t')
        c[(f_[0], f_[1], f_[2])] += int(f_[4])
    nj = np.array(sorted(c.values()), dtype=np.int64)
    N = int(JMETA[jacc]['records_q255'])
    junc[w] = {
        'N': N,
        'total': int(len(nj)),
        'curve': [(d, float(_st2.binom.sf(KJ - 1, nj, d / N).sum())) for d in JDEP if d <= N],
    }
if junc:
    fig, ax = plt.subplots(figsize=(2.8, 2.3))
    for w in WELLS:
        if w in junc and junc[w]['curve']:
            ax.plot(
                [d for d, _ in junc[w]['curve']],
                [v for _, v in junc[w]['curve']],
                color=col(w),
                lw=0.9,
                ls=(0, (3, 2)) if w in RIS500 else '-',
                marker='o',
                ms=1.8,
            )
    ax.set_xscale('log')
    ax.set_xticks(JDEP)
    ax.set_xticklabels(['25k', '50k', '100k', '250k', '500k', '1M', '2M'])
    ax.set_xlabel('deduplicated molecules per well\n(subsampled)')
    ax.set_ylabel(f'unique junctions (>= {KJ} molecules)')
    ax.set_ylim(0, None)
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v/1e3:g}k' if v >= 1e3 else f'{v:g}'))
    glegend(ax, 'upper left')
    save(fig, f'S15_junctions_vs_depth_min{KJ}_{TAG}')
    with open(f'{V}/junctions_min{KJ}_{TAG}.tsv', 'w') as f:
        f.write('bobcode\tsample\tdedup_molecules_total\tunique_junctions_total\tdepth\tjunctions_expected\n')
        [
            f.write(f"{w}\t{LAB[w]}\t{junc[w]['N']}\t{junc[w]['total']}\t{d}\t{v:.1f}\n")
            for w in WELLS
            if w in junc
            for d, v in junc[w]['curve']
        ]
# ---- values + README ----
with open(f'{V}/per_well_{TAG}.tsv', 'w') as f:
    f.write(
        'bobcode\tsample\tgroup\tinput_failure\treads_in\tmapped\tunique\tfiltered_reads\tmolecules_UCI\tgenes_at_native_incl_ribosomal_protein\tduplicate_rate_pct\t'
        + '\t'.join(f'pct_{k}' for k, _, _ in CLASSES)
        + '\tmedian_aligned_len_mate1\tmedian_aligned_len_mate2\ttotal_aligned_bases\tbases_per_fragment\tinsert_median\tinsert_p10\tinsert_p90\tpicard_balance\tbarcode_accuracy_pct\tbarcode_accuracy_n\n'
    )
    for i, w in enumerate(WELLS):
        c = core[w]
        cells = [
            w,
            c['sample'],
            c['group'],
            str(int(c['input_failure'])),
            str(c['reads_in']),
            str(c['mapped']),
            str(c['unique']),
            str(c['filtered']),
            f"{c['molecules_UCI']:.0f}",
            str(c['genes_at_native']),
            f"{c['duplicate_rate_pct']:.2f}",
        ]
        cells += [f'{compv[k][i]:.2f}' if LAB[w] in comp else '' for k, _, _ in CLASSES]
        cells += [
            f"{rl[w]['m1']:.0f}" if w in rl else '',
            f"{rl[w]['m2']:.0f}" if w in rl and rl[w]['m2'] else '',
            str(ab[w]['total']) if w in ab else '',
            f"{ab[w]['per_fragment']:.1f}" if w in ab else '',
        ]
        cells += (
            [
                f"{ins[w]['insert_median']:.0f}",
                f"{ins[w]['insert_quantiles']['10']:.0f}",
                f"{ins[w]['insert_quantiles']['90']:.0f}",
            ]
            if w in ins
            else ['', '', '']
        )
        cells += [
            f"{pic[w]['balance']:.3f}" if w in pic else '',
            f"{acc[w][0]:.3f}" if w in acc else '',
            str(acc[w][3]) if w in acc else '',
        ]
        f.write('\t'.join(cells) + '\n')
with open(f'{V}/depth_curves_{TAG}.tsv', 'w') as f:
    f.write('bobcode\tsample\tdepth\tmolecules_UCI\tgenes\tgenes_incl_ribosomal_protein\n')
    [
        f.write(f"{r['well']}\t{LAB[r['well']]}\t{r['depth']}\t{r['B_corr']}\t{r['genes_B_noRP']}\t{r['genes_B']}\n")
        for r in rar
        if r['well'] in LAB and r.get('B_corr') and r['at_native'] == '0'
    ]
if pic:
    with open(f'{V}/picard_profiles_{TAG}.tsv', 'w') as f:
        f.write('sample\t' + '\t'.join(f'pct{i:03d}' for i in range(101)) + '\n')
        [f.write(LAB[w] + '\t' + '\t'.join(f'{v:.5f}' for v in pic[w]['profile']) + '\n') for w in WELLS if w in pic]
open(f'{OUT}/README.md', 'w').write(
    f"# Per-sample supplement, all 21 human BOBseq wells, {SET} set ({__import__('datetime').date.today()})\n\nOne panel per metric with every human well of the 24-plex run: the 18 benchmark wells plus the three risdiplam 500 nM wells (input failures, excluded from every main figure; hatched bars / dashed lines, grey). Wells are grouped by treatment (colour), Hek control 2 is the well with the T-poor UMI oligo. Sources and rules are those of the main panels (the 18 wells read the same tables; the Ris 500 nM wells were run through the same arm path: mate-1 records of their raw pair BAMs, the canonical read table, the E/U slice, rarefaction, the molecule rule; run_coverage.sh, supp21_tables.py). Native set only (per-sample means BOBseq only; the 50-nt version was dropped). Panels: S01 read fate, S02 filtered reads and molecules, S03 duplicate rate as sequenced, S03b at 250k filtered reads, S04 molecules vs depth (S04b log axis), S05 genes vs depth, S06 composition, S07 intronic, S08 read length per mate, S09 insert length, S10 aligned bases, S11 Picard profile (y 0 to 7 as in the main panel), S12 Picard balance, S14 poly(A)-anchored coverage profile, S15 unique junctions (>= 3 molecules) vs deduplicated molecules (the per-sample junction BEDs of extract_junctions.sh). The per-well barcode accuracy (final funnel) is in the values table only.\n\nvalues/per_well_{TAG}.tsv has every number per well; depth_curves_{TAG}.tsv, picard_profiles_{TAG}.tsv and polyA_profiles_{TAG}.tsv the curves.\n"
)
print('SUPPLEMENT DONE', SET)

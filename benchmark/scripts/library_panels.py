#!/usr/bin/env python3
"""Library panels: read length, insert length, base quality by cycle. BM_COVSET=native|50nt.
  p7a_read_length_<tag>          aligned read length per method (BOBseq per mate), % of unique reads on canonical transcripts (library_metrics.py readlen)
  p7b_insert_length_bobseq       native only: insert length on the mature transcript of deduplicated BOBseq pairs, 18 benchmark wells pooled,
                                 per-well medians in the json (library_metrics.py insert)
  p7c_quality_by_cycle           mean Phred per cycle: BOBseq read 1, BOBseq read 2, and a poly(dT) library read that runs through the poly(A)
                                 (read 2 of an internal poly(dT) library, quality illustration only) (quality_by_cycle.py)
  p7d_base_share_by_cycle        A and T share per cycle of the same reads: where the homopolymer sits and where the quality collapses
Values TSVs in values/, json sidecars beside the panels."""

from settings import WORK, PROJECT, COVERAGE
import os, sys, json, csv, glob, numpy as np, matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, H)
from palette import COL, METH

P = WORK
PR = PROJECT
SET = os.environ.get('BM_COVSET', 'native')
TAG = 'native' if SET == 'native' else 'matched'
FIG = 'preprint_figures_native' if SET == 'native' else f'preprint_figures_{SET}'
OUT = f'{P}/{FIG}/main_figure_panels'
V = f'{OUT}/values'
os.makedirs(V, exist_ok=True)
A = COVERAGE
LM = f'{PR}/results/library_metrics'
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
import re


def save(fig, name):
    fig.tight_layout(pad=0.4)
    fig.savefig(f'{OUT}/{name}.svg', bbox_inches='tight', pad_inches=0.03)
    fig.savefig(f'{OUT}/{name}.png', dpi=200, bbox_inches='tight', pad_inches=0.03)
    plt.close(fig)
    s = open(f'{OUT}/{name}.svg').read()
    s = re.sub(r"font-family:\s*[^;\"]+", 'font-family: Arial', s)
    open(f'{OUT}/{name}.svg', 'w').write(s)
    print('wrote', name)


def legend(ax, handles, labels, loc='upper left', anchor=None):
    ax.legend(
        handles=handles,
        labels=labels,
        loc=loc,
        bbox_to_anchor=anchor,
        frameon=False,
        fontsize=6.5,
        handlelength=1.6,
        borderaxespad=0.3,
        labelspacing=0.3,
    )


LIGHT = '#7FB07F'  # BOBseq mate 2 (read 1)
# ---- p7a read length ----
R = json.load(open(f'{A}/read_length_{SET}.json'))
fig, ax = plt.subplots(figsize=(2.6, 2.2))
H_ = []
L_ = []
series = [(m, m, COL[m], '-') for m in METH if m != 'BOBseq'] + (
    [
        ('BOBseq mate 1', 'BOBseq read 2 (mate 1)', COL['BOBseq'], '-'),
        ('BOBseq mate 2', 'BOBseq read 1 (mate 2)', LIGHT, '-'),
    ]
    if 'BOBseq mate 1' in R
    else [('BOBseq', 'BOBseq (50-nt read 2)', COL['BOBseq'], '-')]
)
with open(f'{V}/p7a_read_length_{TAG}.tsv', 'w') as f:
    f.write('series\taligned_length_nt\tpct_of_reads\n')
    for key, lab, c, ls in series:
        h = np.array(R[key]['hist_0_to_301'], float)
        pct = 100 * h / h.sum()
        x = np.arange(len(h))
        sel = x <= 160
        ax.plot(x[sel], pct[sel], color=c, lw=1.3, ls=ls)
        H_.append(Line2D([], [], color=c, lw=1.3))
        L_.append(f"{lab} (median {R[key]['median']} nt)")
        for xi, p in zip(x[sel], pct[sel]):
            f.write(f'{lab}\t{xi}\t{p:.4f}\n')
ax.set_xlabel('aligned read length (nt)')
ax.set_ylabel('% of unique reads')
ax.set_xlim(0, 160)
ax.set_ylim(0, None)
legend(ax, H_, L_, 'upper left', (0.4, 1.0) if SET == 'native' else (0.45, 1.0))
save(fig, f'p7a_read_length_{TAG}')  # legend right of the 50-nt spikes, above the BOBseq peaks
json.dump(
    {k: {kk: vv for kk, vv in v.items() if kk != 'hist_0_to_301'} for k, v in R.items()},
    open(f'{OUT}/p7a_read_length_{TAG}.json', 'w'),
    indent=1,
)
# ---- p7e / p7f aligned bases (total bases aligned = read number x mean aligned length) ----
AB = list(csv.DictReader(open(f'{A}/aligned_bases_{SET}.tsv'), delimiter='\t'))


def strip_panel(name, key, ylabel, log=False, ylim=None, yticks=None, scale=1.0):
    fig, ax = plt.subplots(figsize=(2.2, 2.2))
    nums = {}
    for i, m in enumerate(METH):
        v = np.array([float(r[key]) for r in AB if r['method'] == m]) * scale
        jit = np.random.default_rng(1).normal(0, 0.07, len(v))
        ax.scatter(i + jit, v, s=7, color=COL[m], alpha=0.8, lw=0, zorder=2)
        ax.hlines(np.median(v), i - 0.28, i + 0.28, color='#231F20', lw=1.0, zorder=3)
        nums[m] = dict(
            n=int(len(v)), median=float(np.median(v)), mean=float(v.mean()), min=float(v.min()), max=float(v.max())
        )
    if log:
        ax.set_yscale('log')
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    if yticks:
        ax.set_yticks(yticks[0])
        ax.set_yticklabels(yticks[1])
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_xticks(range(len(METH)))
    ax.set_xticklabels(METH, rotation=30, ha='right')
    ax.set_xlim(-0.6, len(METH) - 0.4)
    ax.set_ylabel(ylabel)
    save(fig, name)
    json.dump(nums, open(f'{OUT}/{name}.json', 'w'), indent=1)
    print(name, {m: round(nums[m]['median'], 2) for m in METH})


import matplotlib.ticker

D250 = 250000
for r in AB:
    r['aligned_bases_at_250k'] = (
        f"{float(r['bases_per_fragment']) * D250:.0f}"  # matched depth: 250k filtered fragments x the sample's aligned bases per fragment (depth-independent)
    )
fig, ax = plt.subplots(figsize=(2.6, 2.2))
nums = {}
off = (-0.19, 0.19)
for i, m in enumerate(METH):
    nums[m] = {}
    for k, (key, o) in enumerate(zip(('total_aligned_bases', 'aligned_bases_at_250k'), off)):
        v = np.array([float(r[key]) for r in AB if r['method'] == m]) * 1e-9
        jit = np.random.default_rng(1 + k).normal(0, 0.05, len(v))
        if k == 0:
            ax.scatter(i + o + jit, v, s=7, color=COL[m], alpha=0.8, lw=0, zorder=2)
        else:
            ax.scatter(i + o + jit, v, s=8, facecolors='white', edgecolors=COL[m], lw=0.6, alpha=0.9, zorder=2)
        ax.hlines(np.median(v), i + o - 0.15, i + o + 0.15, color='#231F20', lw=1.0, zorder=3)
        nums[m][key] = dict(
            n=int(len(v)), median=float(np.median(v)), mean=float(v.mean()), min=float(v.min()), max=float(v.max())
        )
ax.set_yscale('log')
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_yticks([0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2])
ax.set_yticklabels(['0.01', '0.02', '0.05', '0.1', '0.2', '0.5', '1', '2'])
ax.set_xticks(range(len(METH)))
ax.set_xticklabels(METH, rotation=30, ha='right')
ax.set_xlim(-0.6, len(METH) - 0.4)
ax.set_ylabel('aligned bases per sample (Gb)\nunique reads')
ax.legend(
    handles=[
        Line2D([], [], marker='o', ls='', color='#555', markersize=3.5),
        Line2D([], [], marker='o', ls='', markerfacecolor='white', markeredgecolor='#555', markersize=3.8),
    ],
    labels=['as sequenced', 'at 250k filtered reads\n(matched depth)'],
    loc='upper left',
    frameon=False,
    fontsize=6,
    handletextpad=0.3,
    borderaxespad=0.3,
    labelspacing=0.3,
)
save(fig, f'p7e_aligned_bases_total_{TAG}')
json.dump(nums, open(f'{OUT}/p7e_aligned_bases_total_{TAG}.json', 'w'), indent=1)
print('p7e', {m: {k: round(nums[m][k]['median'], 3) for k in nums[m]} for m in METH})
strip_panel(
    f'p7f_aligned_bases_per_fragment_{TAG}',
    'bases_per_fragment',
    'aligned bases per unique fragment (nt)\nBOBseq: both mates',
    ylim=(0, None),
)
with open(f'{V}/p7e_p7f_aligned_bases_per_sample_{TAG}.tsv', 'w') as f:
    f.write(
        'method\tsample\treads\tmean_aligned_length\ttotal_aligned_bases\tfragments\tbases_per_fragment\taligned_bases_at_250k\n'
    )
    [
        f.write(
            '\t'.join(
                r[k]
                for k in (
                    'method',
                    'sample',
                    'reads',
                    'mean_aligned_length',
                    'total_aligned_bases',
                    'fragments',
                    'bases_per_fragment',
                    'aligned_bases_at_250k',
                )
            )
            + '\n'
        )
        for r in AB
    ]
# ---- p7b insert length (native, BOBseq pairs) ----
if SET == 'native':
    files = sorted(glob.glob(f'{LM}/insert_per_well/*.json'))
    W = {}
    for fn in files:
        d = json.load(open(fn))
        s = d['sample']
        if s.startswith('Ris_dose_500mM') or s.startswith('RAW_control'):
            continue  # benchmark wells only (18)
        W[s] = d
    if W:
        h = np.sum([np.array(d['hist_1nt_to_3000'], float) for d in W.values()], axis=0)
        pct = 100 * h / h.sum()
        x = np.arange(len(h))
        cum = np.cumsum(pct)
        med = int(np.searchsorted(cum, 50))
        p10 = int(np.searchsorted(cum, 10))
        p90 = int(np.searchsorted(cum, 90))
        n_pairs = int(sum(d['pairs_same_transcript'] for d in W.values()))
        fig, ax = plt.subplots(figsize=(2.6, 2.2))
        sel = (x >= 1) & (x <= 600)
        ax.fill_between(x[sel], 0, pct[sel], color=COL['BOBseq'], alpha=0.35, lw=0)
        ax.plot(x[sel], pct[sel], color=COL['BOBseq'], lw=1.0)
        ax.axvline(med, color='#231F20', lw=0.7, ls=(0, (4, 3)))
        ax.text(
            med + 25,
            ax.get_ylim()[1] * 0.92 if ax.get_ylim()[1] else 0.1,
            f'median {med} nt\np10 {p10}, p90 {p90}',
            fontsize=6,
            va='top',
        )
        ax.set_xlabel('insert length on the mature transcript (nt)')
        ax.set_ylabel('% of deduplicated pairs')
        ax.set_xlim(0, 600)
        ax.set_ylim(0, None)
        ax.set_title(f'BOBseq, {len(W)} wells, {n_pairs:,} pairs', fontsize=7, loc='left')
        save(fig, 'p7b_insert_length_bobseq')
        json.dump(
            {
                'wells': {s: {k: v for k, v in d.items() if k != 'hist_1nt_to_3000'} for s, d in W.items()},
                'pooled': {'pairs': n_pairs, 'median': med, 'p10': p10, 'p90': p90},
            },
            open(f'{OUT}/p7b_insert_length_bobseq.json', 'w'),
            indent=1,
        )
        with open(f'{V}/p7b_insert_length_bobseq.tsv', 'w') as f:
            f.write('insert_length_nt\tpct_of_pairs\n')
            [f.write(f'{xi}\t{p:.5f}\n') for xi, p in zip(x, pct) if xi >= 1]
        with open(f'{V}/p7b_insert_length_bobseq_per_well.tsv', 'w') as f:
            f.write('sample\tpairs\tinsert_median\tinsert_p10\tinsert_p90\ttlen_genomic_median\n')
            [
                f.write(
                    f"{s}\t{d['pairs_same_transcript']}\t{d['insert_median']:.0f}\t{d['insert_quantiles']['10']:.0f}\t{d['insert_quantiles']['90']:.0f}\t{d['tlen_genomic_median']:.0f}\n"
                )
                for s, d in W.items()
            ]
        print('insert: wells', len(W), 'median', med, 'p10', p10, 'p90', p90)
    else:
        print('insert: no per-well json yet')
# ---- p7c / p7d quality and base share by cycle ----
Qd = f'{LM}/quality_by_cycle'
QS = [
    ('bobseq_R1', 'BOBseq read 1', COL['BOBseq'], '-'),
    ('bobseq_R2', 'BOBseq read 2', LIGHT, '-'),
    ('polydT_R2', 'poly(dT) library, read through the poly(A)', '#6E6E6E', '-'),
]
QS = [(k, lab, c, ls) for k, lab, c, ls in QS if os.path.exists(f'{Qd}/{k}.json')]
if QS and SET == 'native':
    fig, ax = plt.subplots(figsize=(2.9, 2.2))
    H_ = []
    L_ = []
    with open(f'{V}/p7c_quality_by_cycle.tsv', 'w') as f:
        f.write('series\tcycle\tmean_phred\tfrac_q30\tfrac_A\tfrac_T\n')
        for k, lab, c, ls in QS:
            d = json.load(open(f'{Qd}/{k}.json'))
            q = np.array(d['mean_q'])
            x = np.arange(1, len(q) + 1)
            ax.plot(x, q, color=c, lw=1.3, ls=ls)
            H_.append(Line2D([], [], color=c, lw=1.3))
            L_.append(lab)
            for i in range(len(q)):
                f.write(
                    f"{lab}\t{i+1}\t{q[i]:.3f}\t{d['frac_q30'][i]:.4f}\t{d['base_share']['A'][i]:.4f}\t{d['base_share']['T'][i]:.4f}\n"
                )
    ax.axhline(30, color='#7f7f7f', lw=0.7, ls=(0, (4, 3)))
    ax.set_xlabel('sequencing cycle')
    ax.set_ylabel('mean base quality (Phred)')
    ax.set_xlim(0, 151)
    ax.set_ylim(0, 42)
    legend(ax, H_, L_, 'lower left')
    save(fig, 'p7c_quality_by_cycle')
    fig, ax = plt.subplots(figsize=(2.9, 2.2))
    H_ = []
    L_ = []
    for k, lab, c, ls in QS:
        d = json.load(open(f'{Qd}/{k}.json'))
        x = np.arange(1, d['cycles'] + 1)
        ax.plot(x, 100 * np.array(d['base_share']['A']), color=c, lw=1.3)
        ax.plot(x, 100 * np.array(d['base_share']['T']), color=c, lw=1.3, ls=(0, (2, 2)))
        H_.append(Line2D([], [], color=c, lw=1.3))
        L_.append(lab)
    H_ += [Line2D([], [], color='#231F20', lw=1.0), Line2D([], [], color='#231F20', lw=1.0, ls=(0, (2, 2)))]
    L_ += ['A share', 'T share']
    ax.set_xlabel('sequencing cycle')
    ax.set_ylabel('% of bases')
    ax.set_xlim(0, 151)
    ax.set_ylim(0, 100)
    legend(ax, H_, L_, 'upper right')
    save(fig, 'p7d_base_share_by_cycle')
    json.dump(
        {
            k: {
                'label': lab,
                'reads': json.load(open(f'{Qd}/{k}.json'))['reads'],
                'cycles': json.load(open(f'{Qd}/{k}.json'))['cycles'],
                'mean_q_at': {
                    c: json.load(open(f'{Qd}/{k}.json'))['mean_q'][c - 1]
                    for c in (1, 25, 50, 75, 100, 125, 150)
                    if c <= json.load(open(f'{Qd}/{k}.json'))['cycles']
                },
            }
            for k, lab, _, _ in QS
        },
        open(f'{OUT}/p7c_quality_by_cycle.json', 'w'),
        indent=1,
    )
print('DONE')

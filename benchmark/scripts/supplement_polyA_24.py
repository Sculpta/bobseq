#!/usr/bin/env python3
"""Coverage vs distance from the poly(A) site for all 24 samples of the 24-plex in one panel: the 21 human wells and the 3 mouse wells.
One species-neutral rule for every well (the human-only S14 uses the benchmark gene set shared by the three methods, which has no mouse counterpart): canonical
protein-coding transcripts >= 1 kb, ribosomal-protein genes excluded, genes with >= 20 unique reads in that well; per gene coverage / gene mean anchored at the
3' end, 3-kb window, mean over genes, 25-nt smoothing (as S14). Human profiles are computed here; mouse profiles come from supplement_per_sample_mouse.py (same rule).
Output: supplement_per_sample_21/S14_coverage_vs_nt_from_polyA_all24_native.(svg|png) (the S14 of the supplement; the single-species versions go to superseded/) + values/polyA_profiles_all24_native.tsv
"""

from settings import WORK, COVERAGE, RUN_JSON
import os, re, sys, csv, json, numpy as np, matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from concurrent.futures import ProcessPoolExecutor

H = os.path.dirname(os.path.abspath(__file__))
A = COVERAGE
P = WORK
OUT = f'{P}/preprint_figures_native/supplement_per_sample_21'
V = f'{OUT}/values'
RIBO = re.compile(r'^(RP[LS]|RPLP|MRP[LS])[^K]*$', re.I)
W3 = 3000
MINR = 20
RUN = json.load(open(RUN_JSON))
LAB = {
    w: l.replace(' ', '_').replace('/', '_')
    for w, l in RUN['bobcode_labels'].items()
    if RUN['tso_species_map'][w] == 'human'
}


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
LABEL = {'Ris dose 25mM': 'Ris 25 nM', 'Ris dose 500mM': 'Ris 500 nM'}   # the sample labels carry a legacy unit; the risdiplam arms are 25 and 500 nM
GCOL = {
    'Hek control': '#2E6B32',
    'CADM1 1nM': '#4C78A8',
    'CADM1 100nM': '#1F4E79',
    'CHX dose 1ug mL': '#E6B422',
    'CHX dose 50ug mL': '#B8860B',
    'Ris dose 25mM': '#B5379B',
    'Ris dose 500mM': '#8C8C8C',
}
MC = '#7A4FA3'
S = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open(f'{A}/samples.tsv') if l.strip()}
lab_of = {os.path.basename(r[1])[:-4]: r[0] for r in S.values() if r[2] == 'BOBseq'}
for k in (1, 2, 3):
    lab_of[f'Ris_dose_500mM-{k}'] = f'bobseq-ris-500mM-{k}'
NAMES = sorted(LAB.values(), key=lambda n: (GROUPS.index(group(n)), n))
assert len(NAMES) == 21 and all(n in lab_of for n in NAMES)


def sm(y, k=25):
    y = np.asarray(y, float)
    v = np.where(np.isnan(y), 0, y)
    n = (~np.isnan(y)).astype(float)
    ker = np.ones(k)
    return np.convolve(v, ker, 'same') / np.maximum(np.convolve(n, ker, 'same'), 1)


def profile(n):
    z = np.load(f'{A}/positions_canonical/{lab_of[n]}.npz')
    genes = list(z['genes'])
    tl = z['tx_len']
    g = z['gene']
    T = z['t']
    AL = z['alen'].astype(int)
    o = np.argsort(g, kind='stable')
    gs = g[o]
    st = np.flatnonzero(np.r_[True, gs[1:] != gs[:-1]])
    en = np.r_[st[1:], len(gs)]
    U = []
    for a, b in zip(st, en):
        gi = int(gs[a])
        ii = o[a:b]
        L = int(tl[gi])
        if len(ii) < MINR or L < 1000 or RIBO.match(genes[gi]):
            continue
        t = T[ii].astype(int)
        e = np.minimum(t + AL[ii], L)
        d = np.zeros(L + 1)
        np.add.at(d, t, 1)
        np.add.at(d, e, -1)
        c = np.cumsum(d)[:L]
        zz = c / c.mean()
        row = np.full(W3, np.nan)
        w = min(W3, L)
        row[:w] = zz[::-1][:w]
        U.append(row)
    return n, sm(np.nanmean(np.array(U), 0)), len(U)


if __name__ == '__main__':
    with ProcessPoolExecutor(11) as ex:
        res = {n: (p, k) for n, p, k in ex.map(profile, NAMES)}
    mouse = {
        r['sample']: (np.array([float(r[f'nt{d}']) for d in range(W3)]), int(r['genes']))
        for r in csv.DictReader(
            open(f'{P}/preprint_figures_native/supplement_per_sample_mouse/values/polyA_profiles_native_mouse.tsv'),
            delimiter='\t',
        )
    }
    assert len(mouse) == 3
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
            'axes.spines.top': False,
            'axes.spines.right': False,
            'text.color': '#231F20',
            'axes.edgecolor': '#231F20',
            'axes.labelcolor': '#231F20',
            'xtick.color': '#231F20',
            'ytick.color': '#231F20',
        }
    )
    fig, ax = plt.subplots(figsize=(2.8, 2.3))
    xx = np.arange(W3)
    for n in NAMES:
        ax.plot(xx, res[n][0], color=GCOL[group(n)], lw=0.9, ls=(0, (3, 2)) if group(n) == 'Ris dose 500mM' else '-')
    for n in sorted(mouse):
        ax.plot(xx, mouse[n][0], color=MC, lw=0.9)
    ax.axhline(1, color='#9a9a9a', lw=0.6, ls=(0, (2, 2)))
    ax.set_xlim(W3, 0)
    ax.set_ylim(0, 8.5)
    ax.set_xlabel('distance from the poly(A) site (nt)')
    ax.set_ylabel('coverage / gene mean\n(genes >= 1 kb, RP excluded)')
    ax.legend(
        handles=[Patch(color=GCOL[g]) for g in GROUPS] + [Patch(color=MC)],
        labels=[LABEL.get(g, g) + (' (input failure)' if g == 'Ris dose 500mM' else '') for g in GROUPS] + ['RAW 264.7 (mouse)'],
        loc='upper left',
        frameon=False,
        fontsize=5.5,
        handlelength=1.0,
        borderaxespad=0.3,
        labelspacing=0.25,
    )
    name = 'S14_coverage_vs_nt_from_polyA_all24_native'
    fig.tight_layout(pad=0.4)
    fig.savefig(f'{OUT}/{name}.svg', bbox_inches='tight', pad_inches=0.03)
    fig.savefig(f'{OUT}/{name}.png', dpi=200, bbox_inches='tight', pad_inches=0.03)
    s = open(f'{OUT}/{name}.svg').read()
    open(f'{OUT}/{name}.svg', 'w').write(re.sub(r"font-family:\s*[^;\"]+", 'font-family: Arial', s))
    with open(f'{V}/polyA_profiles_all24_native.tsv', 'w') as f:
        f.write('sample\tspecies\tgenes\t' + '\t'.join(f'nt{d}' for d in range(W3)) + '\n')
        [f.write(f'{n}\thuman\t{res[n][1]}\t' + '\t'.join(f'{v:.4f}' for v in res[n][0]) + '\n') for n in NAMES]
        [
            f.write(f'{n}\tmouse\t{mouse[n][1]}\t' + '\t'.join(f'{v:.4f}' for v in mouse[n][0]) + '\n')
            for n in sorted(mouse)
        ]
    print(
        'wrote',
        name,
        '| genes per well: human',
        min(v[1] for v in res.values()),
        'to',
        max(v[1] for v in res.values()),
        '| mouse',
        sorted(v[1] for v in mouse.values()),
    )

# the note that accompanies S14 in the figure tree (organize_figures.py places it beside the panel)
README = "# S14: coverage vs distance from the poly(A) site, all 24 samples of the 24-plex\n\nThe 21 human wells (treatment colours; the three risdiplam 500 nM wells, input failures, dashed) and the three mouse wells (RAW 264.7, purple) in one panel, the only S14 panel; the single-species versions are in superseded/ of the generator output. One rule for every well: canonical protein-coding transcripts of at least 1 kb (human GRCh38 or mouse GRCm39, Ensembl 113), ribosomal-protein genes excluded, genes with at least 20 unique reads in that well; per gene, coverage divided by the gene mean, anchored at the 3' end, 3-kb window, mean over genes, 25-nt smoothing. The y range (0 to 8.5) is that of the benchmark panel p4m2. The benchmark panel itself uses the gene set shared by the three methods; on a HEK well the two rules differ by at most 0.1. Values: `values/polyA_profiles_all24_native.tsv` (genes per well in the table). Script: `supplement_polyA_24.py`; mouse inputs from `mouse3_chain.sh`.\n"
open(f'{OUT}/README_S14_all24.md', 'w').write(README)
print('wrote README_S14_all24.md')

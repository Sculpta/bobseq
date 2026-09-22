#!/usr/bin/env python3
"""Pooled (dataset-level) unique annotated junctions vs depth, native benchmark: all samples of a method pooled on the same x axis, because the
competitors keep re-detecting the same junctions, so the per-sample view understates the difference. Per method the molecule counts of every junction are summed over its benchmark
samples (DRUG-seq 24, prime-seq 8, BOBseq 18), N = the pooled UMI-deduplicated MAPQ-255 molecules; the pool is thinned binomially to each depth (exact in expectation:
E[junctions with >= k molecules at d] = sum_j P(Binom(n_j, d/N) >= k)); a junction counts when its intron is in the Ensembl 113 annotation (junction_rarefaction_annotated.py's
cache). Output: junctions/pooled/ under WORK (panel per threshold, overview, values, summary with fold values).
"""

from settings import WORK, PROJECT
import os, re, sys, csv, gzip, json, collections, numpy as np, matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats as st
from matplotlib.lines import Line2D

H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, H)
from palette import COL, METH

PR = PROJECT
D = f'{PR}/data/junctions'
P = WORK
OUT = f'{WORK}/junctions/pooled'
os.makedirs(OUT, exist_ok=True)
ANN = {
    tuple(l.rstrip('\n').split('\t'))
    for l in gzip.open(f'{WORK}/junctions/annotated/annotated_introns_GRCh38.113.tsv.gz', 'rt')
}
SETS = {'DRUG-seq': 'drugseq_native', 'prime-seq': 'primeseq_native', 'BOBseq': 'bobseq_pe_native'}
K = [1, 2, 3]
meta = {r['acc_number']: r for r in csv.DictReader(open(f'{D}/sample_metadata.csv'))}
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
    fig.tight_layout(pad=0.4)
    fig.savefig(f'{OUT}/{name}.svg', bbox_inches='tight', pad_inches=0.03)
    fig.savefig(f'{OUT}/{name}.png', dpi=200, bbox_inches='tight', pad_inches=0.03)
    plt.close(fig)
    s = open(f'{OUT}/{name}.svg').read()
    open(f'{OUT}/{name}.svg', 'w').write(re.sub(r"font-family:\s*[^;\"]+", 'font-family: Arial', s))


pool = {}
POOL_ALL = {}
info = {}
for m, s_ in SETS.items():
    accs = [
        a
        for a, r in meta.items()
        if r['set'] == s_ and not (m == 'BOBseq' and r['sample_name'].startswith('Ris_dose_500mM'))
    ]
    c = collections.Counter()
    call = collections.Counter()
    N = 0
    per_sample_ann = []
    for a in accs:
        cs = collections.Counter()
        for l in gzip.open(f'{D}/{a}_junctions.bed.gz', 'rt'):
            f = l.rstrip('\n').split('\t')
            k_ = (f[0], f[1], f[2])
            call[k_] += int(f[4])
            if k_ in ANN:
                cs[k_] += int(f[4])
        c.update(cs)
        N += int(meta[a]['records_q255'])
        per_sample_ann.append(len(cs))
    nj = np.array(sorted(c.values()), dtype=np.int64)
    pool[m] = (nj, N)
    POOL_ALL[m] = (np.array(sorted(call.values()), dtype=np.int64), N)
    info[m] = {
        'samples': len(accs),
        'pooled_molecules': N,
        'pooled_annotated_junctions_ge1': int(len(nj)),
        'pooled_ge3': int((nj >= 3).sum()),
        'pooled_all_junctions_ge1': int(len(call)),
        'pooled_all_ge3': int((POOL_ALL[m][0] >= 3).sum()),
        'median_per_sample_annotated_junctions_ge1': float(np.median(per_sample_ann)),
    }
    print(m, info[m], flush=True)
DEP = [25000, 50000, 100000, 250000, 500000, 1000000, 2000000, 5000000, 10000000, 20000000, 50000000, 100000000]
DL = ['25k', '50k', '100k', '250k', '500k', '1M', '2M', '5M', '10M', '20M', '50M', '100M']


def expected(nj, N, d, k):
    return np.nan if d > N else float(st.binom.sf(k - 1, nj, d / N).sum())


E = {m: {k: [expected(*pool[m], d, k) for d in DEP] for k in K} for m in METH}
fig, axs = plt.subplots(1, len(K), figsize=(2.6 * len(K), 2.5))
for ki, k in enumerate(K):
    for ax in (axs[ki],):
        for m in METH:
            y = np.array(E[m][k])
            ok = ~np.isnan(y)
            ax.plot(
                list(np.array(DEP)[ok]) + [pool[m][1]],
                list(y[ok]) + [expected(*pool[m], pool[m][1], k)],
                color=COL[m],
                lw=1.3,
                marker='o',
                ms=2.2,
            )  # the last point = the full pooled depth of the method
        ax.set_xscale('log')
        ax.set_xlabel('deduplicated molecules, all samples of a method pooled')
        ax.set_ylim(0, None)
        ax.set_title(f'>= {k} molecule{"s" if k > 1 else ""} per junction', fontsize=7)
        ax.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda v, _: f'{v/1e3:g}k' if v >= 1e3 else f'{v:g}')
        )
        if ki == 0:
            ax.set_ylabel('unique annotated junctions (pooled)')
            ax.legend(
                handles=[Line2D([], [], color=COL[m], lw=1.3, marker='o', ms=2.2) for m in METH],
                labels=[f'{m} ({info[m]["samples"]} samples)' for m in METH],
                loc='upper left',
                frameon=False,
                fontsize=6,
            )
    fig2, ax2 = plt.subplots(figsize=(2.8, 2.3))
    for m in METH:
        y = np.array(E[m][k])
        ok = ~np.isnan(y)
        ax2.plot(
            list(np.array(DEP)[ok]) + [pool[m][1]],
            list(y[ok]) + [expected(*pool[m], pool[m][1], k)],
            color=COL[m],
            lw=1.3,
            marker='o',
            ms=2.2,
        )
    ax2.set_xscale('log')
    ax2.set_xticks([1e5, 1e6, 1e7, 1e8])
    ax2.set_xticklabels(['100k', '1M', '10M', '100M'])
    ax2.set_xlabel('deduplicated molecules\n(all samples of a method pooled, subsampled)')
    ax2.set_ylabel(f'unique annotated junctions\n(>= {k} molecule{"s" if k > 1 else ""} per junction, pooled)')
    ax2.set_ylim(0, None)
    ax2.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v/1e3:g}k' if v >= 1e3 else f'{v:g}'))
    ax2.legend(
        handles=[Line2D([], [], color=COL[m], lw=1.3, marker='o', ms=2.2) for m in METH],
        labels=[f'{m} (n={info[m]["samples"]})' for m in METH],
        loc='upper left',
        frameon=False,
        fontsize=6,
        handlelength=1.8,
        borderaxespad=0.3,
        labelspacing=0.3,
    )
    save(fig2, f'p5_pooled_annotated_junctions_vs_depth_min{k}_native')
    if k == 3:  # main figure copy
        MAIN = f'{P}/preprint_figures_native/main_figure_panels'
        os.makedirs(f'{MAIN}/values', exist_ok=True)
        nm = f'p5_pooled_annotated_junctions_vs_depth_min{k}_native'
        import shutil

        [shutil.copy2(f'{OUT}/{nm}.{ext}', f'{MAIN}/{nm}.{ext}') for ext in ('svg', 'png')]
        json.dump(
            {
                'set': 'native',
                'min_molecules': k,
                'pooled': True,
                'per_method': {m: {**info[m], 'at_full_pool': expected(*pool[m], pool[m][1], k)} for m in METH},
                'source': 'as the per-sample annotated panel, all samples of a method pooled before thinning',
            },
            open(f'{MAIN}/{nm}.json', 'w'),
            indent=1,
        )
        with open(f'{MAIN}/values/{nm}.tsv', 'w') as fv:
            fv.write('method\tsamples_pooled\tpooled_molecules\tdepth\tannotated_junctions_expected\n')
            [
                fv.write(f'{m}\t{info[m]["samples"]}\t{pool[m][1]}\t{d}\t{E[m][k][i]:.1f}\n')
                for m in METH
                for i, d in enumerate(DEP)
                if not np.isnan(E[m][k][i])
            ]
            [
                fv.write(
                    f'{m}\t{info[m]["samples"]}\t{pool[m][1]}\t{pool[m][1]}\t{expected(*pool[m], pool[m][1], k):.1f}\n'
                )
                for m in METH
            ]
fig.suptitle(
    'pooled unique annotated junctions vs pooled depth, native benchmark; the last point of a curve = the method at its full pooled depth',
    fontsize=7.5,
)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(f'{OUT}/overview_pooled_junctions_native.png', dpi=170, bbox_inches='tight')
plt.close(fig)
with open(f'{OUT}/values_pooled_junction_rarefaction.tsv', 'w') as f:
    f.write('method\tsamples_pooled\tpooled_molecules\tdepth\tmin_molecules\tannotated_junctions_expected\n')
    [
        f.write(f'{m}\t{info[m]["samples"]}\t{pool[m][1]}\t{d}\t{k}\t{E[m][k][i]:.1f}\n')
        for m in METH
        for k in K
        for i, d in enumerate(DEP)
        if not np.isnan(E[m][k][i])
    ]
    [
        f.write(
            f'{m}\t{info[m]["samples"]}\t{pool[m][1]}\t{pool[m][1]}\t{k}\t{expected(*pool[m], pool[m][1], k):.1f}\n'
        )
        for m in METH
        for k in K
    ]
L = [
    '# Pooled unique annotated junctions vs depth, native benchmark',
    '',
    'All samples of a method pooled (junction molecule counts summed), thinned binomially to each depth; annotated junctions only (Ensembl 113 introns); support = UMI-deduplicated MAPQ-255 molecules.',
    '',
    '| method | samples | pooled molecules | annotated junctions >= 1 (full pool) | >= 3 (full pool) | median per-sample junctions >= 1 |',
    '|---|---|---|---|---|---|',
] + [
    f"| {m} | {info[m]['samples']} | {info[m]['pooled_molecules']:,} | {info[m]['pooled_annotated_junctions_ge1']:,} | {info[m]['pooled_ge3']:,} | {info[m]['median_per_sample_annotated_junctions_ge1']:,.0f} |"
    for m in METH
]
L += [
    '',
    '| min molecules | pooled depth | ' + ' | '.join(METH) + ' | BOBseq / DRUG-seq | BOBseq / prime-seq |',
    '|---|---|' + '---|' * (len(METH) + 2),
]
for k in K:
    for i, d in enumerate(DEP):
        v = {m: E[m][k][i] for m in METH}
        if all(np.isnan(v[m]) for m in METH) or np.isnan(v['BOBseq']):
            continue
        L.append(
            f'| {k} | {DL[i]} | '
            + ' | '.join('na' if np.isnan(v[m]) else f'{v[m]:,.0f}' for m in METH)
            + ' | '
            + ' | '.join('na' if np.isnan(v[o]) else f"{v['BOBseq'] / v[o]:.1f}x" for o in ('DRUG-seq', 'prime-seq'))
            + ' |'
        )
    L.append(
        f'| {k} | full pool of each | '
        + ' | '.join(f'{expected(*pool[m], pool[m][1], k):,.0f}' for m in METH)
        + ' | '
        + ' | '.join(
            f"{expected(*pool['BOBseq'], pool['BOBseq'][1], k) / expected(*pool[o], pool[o][1], k):.1f}x"
            for o in ('DRUG-seq', 'prime-seq')
        )
        + ' |'
    )
# ---- all junctions (annotated or not), pooled, min 3 ----
k = 3
EA = {m: [expected(*POOL_ALL[m], d, k) for d in DEP] for m in METH}
nm = f'p5_pooled_all_junctions_vs_depth_min{k}_native'
fig2, ax2 = plt.subplots(figsize=(2.8, 2.3))
for m in METH:
    y = np.array(EA[m])
    ok = ~np.isnan(y)
    ax2.plot(
        list(np.array(DEP)[ok]) + [POOL_ALL[m][1]],
        list(y[ok]) + [expected(*POOL_ALL[m], POOL_ALL[m][1], k)],
        color=COL[m],
        lw=1.3,
        marker='o',
        ms=2.2,
    )
ax2.set_xscale('log')
ax2.set_xticks([1e5, 1e6, 1e7, 1e8])
ax2.set_xticklabels(['100k', '1M', '10M', '100M'])
ax2.set_xlabel('deduplicated molecules\n(all samples of a method pooled, subsampled)')
ax2.set_ylabel(f'unique junctions, annotated or not\n(>= {k} molecules per junction, pooled)')
ax2.set_ylim(0, None)
ax2.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v/1e3:g}k' if v >= 1e3 else f'{v:g}'))
ax2.legend(
    handles=[Line2D([], [], color=COL[m], lw=1.3, marker='o', ms=2.2) for m in METH],
    labels=[f'{m} (n={info[m]["samples"]})' for m in METH],
    loc='upper left',
    frameon=False,
    fontsize=6,
    handlelength=1.8,
    borderaxespad=0.3,
    labelspacing=0.3,
)
save(fig2, nm)
MAIN = f'{P}/preprint_figures_native/main_figure_panels'
import shutil

[shutil.copy2(f'{OUT}/{nm}.{ext}', f'{MAIN}/{nm}.{ext}') for ext in ('svg', 'png')]
json.dump(
    {
        'set': 'native',
        'min_molecules': k,
        'pooled': True,
        'junctions': 'all (annotated or not)',
        'per_method': {
            m: {
                'samples': info[m]['samples'],
                'pooled_molecules': info[m]['pooled_molecules'],
                'pooled_all_junctions_ge1': info[m]['pooled_all_junctions_ge1'],
                'at_full_pool': expected(*POOL_ALL[m], POOL_ALL[m][1], k),
            }
            for m in METH
        },
    },
    open(f'{MAIN}/{nm}.json', 'w'),
    indent=1,
)
with open(f'{MAIN}/values/{nm}.tsv', 'w') as fv:
    fv.write('method\tsamples_pooled\tpooled_molecules\tdepth\tall_junctions_expected\n')
    [
        fv.write(f'{m}\t{info[m]["samples"]}\t{POOL_ALL[m][1]}\t{d}\t{EA[m][i]:.1f}\n')
        for m in METH
        for i, d in enumerate(DEP)
        if not np.isnan(EA[m][i])
    ]
    [
        fv.write(
            f'{m}\t{info[m]["samples"]}\t{POOL_ALL[m][1]}\t{POOL_ALL[m][1]}\t{expected(*POOL_ALL[m], POOL_ALL[m][1], k):.1f}\n'
        )
        for m in METH
    ]
L += [
    '',
    '## All junctions (annotated or not), pooled, min 3 molecules',
    '',
    '| pooled depth | ' + ' | '.join(METH) + ' | BOBseq / DRUG-seq | BOBseq / prime-seq |',
    '|---|' + '---|' * (len(METH) + 2),
]
for i, d in enumerate(DEP):
    v = {m: EA[m][i] for m in METH}
    if np.isnan(v['BOBseq']):
        continue
    L.append(
        f'| {DL[i]} | '
        + ' | '.join('na' if np.isnan(v[m]) else f'{v[m]:,.0f}' for m in METH)
        + ' | '
        + ' | '.join('na' if np.isnan(v[o]) else f"{v['BOBseq'] / v[o]:.1f}x" for o in ('DRUG-seq', 'prime-seq'))
        + ' |'
    )
L.append(
    '| full pool of each | '
    + ' | '.join(f'{expected(*POOL_ALL[m], POOL_ALL[m][1], k):,.0f}' for m in METH)
    + ' | '
    + ' | '.join(
        f"{expected(*POOL_ALL['BOBseq'], POOL_ALL['BOBseq'][1], k) / expected(*POOL_ALL[o], POOL_ALL[o][1], k):.1f}x"
        for o in ('DRUG-seq', 'prime-seq')
    )
    + ' |'
)
open(f'{OUT}/SUMMARY_pooled_junctions.md', 'w').write('\n'.join(L) + '\n')
json.dump(info, open(f'{OUT}/pooled_numbers.json', 'w'), indent=1)
print('DONE')

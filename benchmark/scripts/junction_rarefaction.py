#!/usr/bin/env python3
"""Cumulative unique junctions detected vs depth, at several minimum-support thresholds, from the per-sample
junction BEDs of extract_junctions.sh: counts are UMI-deduplicated molecules with
MAPQ 255, a pair spanning a junction counts once. Rarefaction is exact in expectation without the BAMs: subsampling a sample's
deduplicated molecules to depth d keeps each junction's n_j molecules binomially with p = d / N (N = the sample's dedup MAPQ-255
molecules, records_q255 in the sample sheet), so E[junctions with >= k molecules at d] = sum_j P(Binom(n_j, p) >= k).
Sets: native (bobseq_pe_native, drugseq_native, primeseq_native) and matched (bobseq_50nt, drugseq, primeseq); the 18 benchmark BOBseq
wells (Ris 500 nM excluded, Hek control 2 included, as in every main panel). Output: junctions/all_junctions/ (one SVG+PNG per
set and threshold in the p1 style, an overview sheet per set, values TSV)."""

from settings import WORK, PROJECT
import os, re, sys, csv, gzip, collections, numpy as np, matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats as st
from matplotlib.lines import Line2D

H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, H)
from palette import COL, METH

PR = PROJECT
D = f'{PR}/data/junctions'
OUT = f'{WORK}/junctions/all_junctions'
os.makedirs(OUT, exist_ok=True)
SETS = {
    'native': {'DRUG-seq': 'drugseq_native', 'prime-seq': 'primeseq_native', 'BOBseq': 'bobseq_pe_native'},
    'matched': {'DRUG-seq': 'drugseq', 'prime-seq': 'primeseq', 'BOBseq': 'bobseq_50nt'},
}
K = [1, 2, 3, 5, 10]
KMAIN = 3  # the >= 3 molecules panel is also copied to the main figure sets, the others stay in OUT
P = WORK
MAIN = {
    'native': f'{P}/preprint_figures_native/main_figure_panels',
    'matched': f'{P}/preprint_figures_50nt/main_figure_panels',
}
import json

DEPTHS = [25000, 50000, 100000, 250000, 500000, 1000000, 2000000]
DL = ['25k', '50k', '100k', '250k', '500k', '1M', '2M']
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


def save(fig, name, out=None):
    OUT_ = out or OUT
    fig.tight_layout(pad=0.4)
    fig.savefig(f'{OUT_}/{name}.svg', bbox_inches='tight', pad_inches=0.03)
    fig.savefig(f'{OUT_}/{name}.png', dpi=200, bbox_inches='tight', pad_inches=0.03)
    plt.close(fig)
    s = open(f'{OUT_}/{name}.svg').read()
    s = re.sub(r"font-family:\s*[^;\"]+", 'font-family: Arial', s)
    open(f'{OUT_}/{name}.svg', 'w').write(s)


def band(ax, x, y, col, lw=1.3, alpha=0.3):
    ds = [d for d in sorted(set(x)) if (x == d).sum() >= 3]
    mu = np.array([y[x == d].mean() for d in ds])
    n = np.array([(x == d).sum() for d in ds])
    ci = np.array([st.t.ppf(0.975, k - 1) * y[x == d].std(ddof=1) / np.sqrt(k) for d, k in zip(ds, n)])
    ax.fill_between(ds, mu - ci, mu + ci, color=col, alpha=alpha, lw=0, zorder=2)
    ax.plot(ds, mu, color=col, lw=lw, zorder=3)


def counts(acc):
    c = collections.Counter()
    for l in gzip.open(f'{D}/{acc}_junctions.bed.gz', 'rt'):
        f = l.rstrip('\n').split('\t')
        c[(f[0], f[1], f[2])] += int(f[4])  # unique junction = intron coordinates; support = molecules
    return np.array(sorted(c.values()), dtype=np.int64)


def expected(nj, N, d, k):
    if d > N:
        return np.nan
    p = d / N
    return float(st.binom.sf(k - 1, nj, p).sum())


rows = []
VAL = open(f'{OUT}/values_junction_rarefaction.tsv', 'w')
VAL.write(
    'set\tmethod\tsample\tdedup_molecules_total\tunique_junctions_total\tdepth\tmin_molecules\tjunctions_expected\n'
)
for tag, arms in SETS.items():
    R = {}
    for m, s in arms.items():
        accs = [
            a
            for a, r in meta.items()
            if r['set'] == s and not (m == 'BOBseq' and r['sample_name'].startswith('Ris_dose_500mM'))
        ]
        for a in accs:
            nj = counts(a)
            N = int(meta[a]['records_q255'])
            r = {'n_total': N, 'junctions_total': int(len(nj))}
            for k in K:
                for d in DEPTHS:
                    e = expected(nj, N, d, k)
                    r[(k, d)] = e
                    if not np.isnan(e):
                        VAL.write(f'{tag}\t{m}\t{meta[a]["sample_name"]}\t{N}\t{len(nj)}\t{d}\t{k}\t{e:.1f}\n')
            R.setdefault(m, {})[a] = r
        print(
            tag,
            m,
            len(accs),
            'samples; dedup MAPQ255 molecules median',
            int(np.median([R[m][a]['n_total'] for a in R[m]])),
            '; junctions (>=1) median',
            int(np.median([R[m][a]['junctions_total'] for a in R[m]])),
            flush=True,
        )
    ymax = max(v[(1, d)] for m in R for v in R[m].values() for d in DEPTHS if not np.isnan(v[(1, d)])) * 1.05
    fig, axs = plt.subplots(1, len(K), figsize=(2.4 * len(K), 2.4), sharey=True)
    for ki, k in enumerate(K):
        for ax in (axs[ki],):
            for m in METH:
                x = np.array([d for a in R[m] for d in DEPTHS if not np.isnan(R[m][a][(k, d)])])
                y = np.array([R[m][a][(k, d)] for a in R[m] for d in DEPTHS if not np.isnan(R[m][a][(k, d)])])
                band(ax, x, y, COL[m])
            ax.set_xscale('log')
            ax.set_xticks(DEPTHS)
            ax.set_xticklabels(DL, fontsize=6)
            ax.set_xlabel('deduplicated molecules per sample (subsampled)')
            ax.set_ylim(0, ymax)
            ax.set_title(f'junctions with >= {k} molecule{"s" if k > 1 else ""}', fontsize=7)
            if ki == 0:
                ax.set_ylabel('unique junctions detected')
        fig2, ax2 = plt.subplots(figsize=(2.6, 2.2))
        for m in METH:
            x = np.array([d for a in R[m] for d in DEPTHS if not np.isnan(R[m][a][(k, d)])])
            y = np.array([R[m][a][(k, d)] for a in R[m] for d in DEPTHS if not np.isnan(R[m][a][(k, d)])])
            band(ax2, x, y, COL[m])
        ax2.set_xscale('log')
        ax2.set_xticks(DEPTHS)
        ax2.set_xticklabels(DL)
        ax2.set_xlabel('deduplicated molecules per sample\n(subsampled)')
        ax2.set_ylabel(f'unique junctions (>= {k} molecule{"s" if k > 1 else ""})')
        ax2.set_ylim(0, None)
        ax2.legend(
            handles=[Line2D([], [], color=COL[m], lw=1.3) for m in METH],
            labels=[f'{m} (n={len(R[m])})' for m in METH],
            loc='upper left',
            frameon=False,
            fontsize=6.5,
            handlelength=1.6,
            borderaxespad=0.3,
            labelspacing=0.3,
        )
        save(fig2, f'p5_junctions_vs_depth_min{k}_{tag}')
        if k == KMAIN:  # main figure copy + sidecar + values
            fig3, ax3 = plt.subplots(figsize=(2.6, 2.2))
            nums = {}
            for m in METH:
                x = np.array([d for a in R[m] for d in DEPTHS if not np.isnan(R[m][a][(k, d)])])
                y = np.array([R[m][a][(k, d)] for a in R[m] for d in DEPTHS if not np.isnan(R[m][a][(k, d)])])
                band(ax3, x, y, COL[m])
                nums[m] = {
                    'n_samples': len(R[m]),
                    'per_depth': {
                        str(d): {'n': int((x == d).sum()), 'mean': float(y[x == d].mean())}
                        for d in DEPTHS
                        if (x == d).sum()
                    },
                }
            ax3.set_xscale('log')
            ax3.set_xticks(DEPTHS)
            ax3.set_xticklabels(DL)
            ax3.set_xlabel('deduplicated molecules per sample\n(subsampled)')
            ax3.set_ylabel(f'unique junctions (>= {k} molecules)')
            ax3.set_ylim(0, None)
            ax3.legend(
                handles=[Line2D([], [], color=COL[m], lw=1.3) for m in METH],
                labels=[f'{m} (n={len(R[m])})' for m in METH],
                loc='upper left',
                frameon=False,
                fontsize=6.5,
                handlelength=1.6,
                borderaxespad=0.3,
                labelspacing=0.3,
            )
            save(fig3, f'p5_junctions_vs_depth_min{k}_{tag}', out=MAIN[tag])
            json.dump(
                {
                    'set': tag,
                    'min_molecules': k,
                    'source': 'per-sample junction BEDs (extract_junctions.sh), UMI-deduplicated MAPQ-255 molecules, rarefaction by binomial thinning',
                    'per_method': nums,
                },
                open(f'{MAIN[tag]}/p5_junctions_vs_depth_min{k}_{tag}.json', 'w'),
                indent=1,
            )
            os.makedirs(f'{MAIN[tag]}/values', exist_ok=True)
            with open(f'{MAIN[tag]}/values/p5_junctions_vs_depth_min{k}_{tag}.tsv', 'w') as fv:
                fv.write('method\tsample\tdedup_molecules_total\tdepth\tjunctions_expected\n')
                [
                    fv.write(f"{m}\t{meta[a]['sample_name']}\t{R[m][a]['n_total']}\t{d}\t{R[m][a][(k, d)]:.1f}\n")
                    for m in METH
                    for a in R[m]
                    for d in DEPTHS
                    if not np.isnan(R[m][a][(k, d)])
                ]
    fig.suptitle(
        f'unique splice junctions vs depth, {tag} set: thick line = mean over samples, ribbon = 95% t-interval; '
        + ', '.join(f'{m} {COL[m]}' for m in METH)
        + '; junction support = UMI-deduplicated MAPQ-255 molecules, pair spanning a junction = 1',
        fontsize=7.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(f'{OUT}/overview_junctions_{tag}.png', dpi=170, bbox_inches='tight')
    plt.close(fig)
    for k in K:
        for d in (250000, 500000):
            rows.append(
                (tag, k, d)
                + tuple(
                    (
                        f"{np.mean([R[m][a][(k, d)] for a in R[m] if not np.isnan(R[m][a][(k, d)])]):,.0f} (n={sum(1 for a in R[m] if not np.isnan(R[m][a][(k, d)]))})"
                        if any(not np.isnan(R[m][a][(k, d)]) for a in R[m])
                        else 'na'
                    )
                    for m in METH
                )
            )
VAL.close()
with open(f'{OUT}/SUMMARY_junctions.md', 'w') as f:
    f.write(
        '# Unique junctions vs depth, minimum-support sweep\n\nSource: the per-sample junction BEDs of extract_junctions.sh: junction support counted in UMI-deduplicated MAPQ-255 molecules, a pair spanning a junction = 1. Rarefaction is exact in expectation (binomial thinning of each junction\'s count at p = depth / dedup molecules of the sample). Mean over samples (n = samples with at least that many molecules).\n\n| set | min molecules | depth | '
        + ' | '.join(METH)
        + ' |\n|---|---|---|'
        + '---|' * len(METH)
        + '\n'
    )
    for r in rows:
        f.write(f'| {r[0]} | {r[1]} | {r[2]:,} | ' + ' | '.join(r[3:]) + ' |\n')
print('DONE')

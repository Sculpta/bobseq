#!/usr/bin/env python3
"""Annotated-only version of junction_rarefaction.py (strict reference, with the same minimum-support sweep).
A junction counts only when its intron (chrom, start, end) is an intron of a transcript in the Ensembl 113 GTF the reads were aligned to
(all biotypes, all transcripts; BED start = upstream exon end, BED end = downstream exon start - 1). Everything else is as in
junction_rarefaction.py: support = UMI-deduplicated MAPQ-255 molecules (the per-sample junction BEDs of extract_junctions.sh), rarefaction
by binomial thinning with p = depth / the sample's dedup MAPQ-255 molecules (the depth axis is unchanged, so the curves compare directly
with the all-junction ones).
Output: junctions/annotated/ under WORK (panel per set and threshold, overview per set with the all-junction mean dashed, values, summary with fold values).
"""

from settings import WORK, PROJECT, GTF
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
OUT = f'{WORK}/junctions/annotated'
os.makedirs(OUT, exist_ok=True)
SETS = {
    'native': {'DRUG-seq': 'drugseq_native', 'prime-seq': 'primeseq_native', 'BOBseq': 'bobseq_pe_native'},
    'matched': {'DRUG-seq': 'drugseq', 'prime-seq': 'primeseq', 'BOBseq': 'bobseq_50nt'},
}
K = [1, 2, 3, 5, 10]
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
# ---- annotated introns of the human contigs (cached) ----
CACHE = f'{OUT}/annotated_introns_GRCh38.113.tsv.gz'
if not os.path.exists(CACHE):
    ex = collections.defaultdict(list)
    for l in open(GTF):
        if l.startswith('#') or not l.startswith('HUMAN_'):
            continue
        f = l.split('\t', 9)
        if f[2] != 'exon':
            continue
        ex[(f[0][6:], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))].append((int(f[3]), int(f[4])))
    intr = set()
    for (c, t), e in ex.items():
        e.sort()
        for (s1, e1), (s2, e2) in zip(e, e[1:]):
            intr.add((c, e1, s2 - 1))
    with gzip.open(CACHE, 'wt') as fh:
        [fh.write(f'{c}\t{a}\t{b}\n') for c, a, b in sorted(intr)]
ANN = {tuple(l.rstrip('\n').split('\t')) for l in gzip.open(CACHE, 'rt')}
print('annotated introns (human, Ensembl 113, all transcripts):', f'{len(ANN):,}', flush=True)


def save(fig, name):
    fig.tight_layout(pad=0.4)
    fig.savefig(f'{OUT}/{name}.svg', bbox_inches='tight', pad_inches=0.03)
    fig.savefig(f'{OUT}/{name}.png', dpi=200, bbox_inches='tight', pad_inches=0.03)
    plt.close(fig)
    s = open(f'{OUT}/{name}.svg').read()
    s = re.sub(r"font-family:\s*[^;\"]+", 'font-family: Arial', s)
    open(f'{OUT}/{name}.svg', 'w').write(s)


def mean_ci(x, y):
    ds = [d for d in sorted(set(x)) if (x == d).sum() >= 3]
    mu = np.array([y[x == d].mean() for d in ds])
    n = np.array([(x == d).sum() for d in ds])
    ci = np.array([st.t.ppf(0.975, k - 1) * y[x == d].std(ddof=1) / np.sqrt(k) for d, k in zip(ds, n)])
    return ds, mu, ci


def band(ax, x, y, col, lw=1.3, alpha=0.3):
    ds, mu, ci = mean_ci(x, y)
    ax.fill_between(ds, mu - ci, mu + ci, color=col, alpha=alpha, lw=0, zorder=2)
    ax.plot(ds, mu, color=col, lw=lw, zorder=3)


def counts(acc):
    c = collections.Counter()
    for l in gzip.open(f'{D}/{acc}_junctions.bed.gz', 'rt'):
        f = l.rstrip('\n').split('\t')
        c[(f[0], f[1], f[2])] += int(f[4])
    ann = np.array(sorted(v for k, v in c.items() if k in ANN), dtype=np.int64)
    al = np.array(sorted(c.values()), dtype=np.int64)
    return ann, al


def expected(nj, N, d, k):
    return np.nan if d > N else float(st.binom.sf(k - 1, nj, d / N).sum())


# ---- coordinate convention check: the stated offset must be the only one that matches ----
a0 = next(a for a, r in meta.items() if r['set'] == 'bobseq_pe_native' and r['sample_name'].startswith('Hek'))
J = [l.split('\t')[:3] for l in gzip.open(f'{D}/{a0}_junctions.bed.gz', 'rt')]
off = {
    (ds, de): sum((c, str(int(s) + ds), str(int(e) + de)) in ANN for c, s, e in J) / len(J)
    for ds in (-1, 0, 1)
    for de in (-1, 0, 1)
}
best = max(off, key=off.get)
print('offset check on', a0, {k: round(v, 3) for k, v in off.items()}, flush=True)
assert best == (0, 0) and off[best] > 0.5 and sorted(off.values())[-2] < 0.05, 'coordinate convention does not match'
VAL = open(f'{OUT}/values_junction_rarefaction_annotated.tsv', 'w')
VAL.write(
    'set\tmethod\tsample\tdedup_molecules_total\tjunctions_all\tjunctions_annotated\tmolecules_on_annotated_pct\tdepth\tmin_molecules\tannotated_junctions_expected\tall_junctions_expected\n'
)
S = [
    '# Annotated junctions vs depth, minimum-support sweep',
    '',
    f'A junction counts when its intron is in the Ensembl 113 annotation used for alignment ({len(ANN):,} human introns, all transcripts and biotypes). Support = UMI-deduplicated MAPQ-255 molecules (per-sample junction BEDs of extract_junctions.sh, a pair spanning a junction = 1). Rarefaction: binomial thinning at p = depth / dedup molecules of the sample, exact in expectation. Mean over samples with at least that many molecules (n).',
    '',
]
for tag, arms in SETS.items():
    R = {}
    share = {}
    for m, s in arms.items():
        accs = [
            a
            for a, r in meta.items()
            if r['set'] == s and not (m == 'BOBseq' and r['sample_name'].startswith('Ris_dose_500mM'))
        ]
        for a in accs:
            ann, al = counts(a)
            N = int(meta[a]['records_q255'])
            r = {
                'N': N,
                'ja': len(al),
                'jn': len(ann),
                'pct_j': 100 * len(ann) / max(len(al), 1),
                'pct_m': 100 * ann.sum() / max(al.sum(), 1),
            }
            for k in K:
                for d in DEPTHS:
                    r[('n', k, d)] = expected(ann, N, d, k)
                    r[('a', k, d)] = expected(al, N, d, k)
                    if not np.isnan(r[('n', k, d)]):
                        VAL.write(
                            f"{tag}\t{m}\t{meta[a]['sample_name']}\t{N}\t{len(al)}\t{len(ann)}\t{r['pct_m']:.2f}\t{d}\t{k}\t{r[('n', k, d)]:.1f}\t{r[('a', k, d)]:.1f}\n"
                        )
            R.setdefault(m, {})[a] = r
        share[m] = (np.median([v['pct_j'] for v in R[m].values()]), np.median([v['pct_m'] for v in R[m].values()]))
        print(
            tag,
            m,
            len(accs),
            f'samples; annotated share of distinct junctions median {share[m][0]:.1f}%, of junction molecules {share[m][1]:.1f}%',
            flush=True,
        )

    def xy(m, kind, k):
        x = np.array([d for a in R[m] for d in DEPTHS if not np.isnan(R[m][a][(kind, k, d)])])
        y = np.array([R[m][a][(kind, k, d)] for a in R[m] for d in DEPTHS if not np.isnan(R[m][a][(kind, k, d)])])
        return x, y

    ymax = max(v[('a', 1, d)] for m in R for v in R[m].values() for d in DEPTHS if not np.isnan(v[('a', 1, d)])) * 1.05
    fig, axs = plt.subplots(1, len(K), figsize=(2.6 * len(K), 2.5), sharey=False)
    for ki, k in enumerate(K):
        ax = axs[ki]
        for m in METH:
            band(ax, *xy(m, 'n', k), COL[m])
            ds, mu, _ = mean_ci(*xy(m, 'a', k))
            ax.plot(ds, mu, color=COL[m], lw=0.7, ls=(0, (3, 2)), zorder=1)
        ax.set_xscale('log')
        ax.set_xticks(DEPTHS)
        ax.set_xticklabels(DL, fontsize=6)
        ax.set_xlabel('deduplicated molecules per sample (subsampled)')
        ax.set_ylim(0, None)
        ax.set_title(f'>= {k} molecule{"s" if k > 1 else ""} per junction', fontsize=7)
        if ki == 0:
            ax.set_ylabel('unique junctions detected')
        if ki == 0:
            ax.legend(
                handles=[Line2D([], [], color='#555', lw=1.3), Line2D([], [], color='#555', lw=0.7, ls=(0, (3, 2)))],
                labels=['annotated junctions', 'all junctions (current panel)'],
                loc='upper left',
                frameon=False,
                fontsize=6,
            )
        fig2, ax2 = plt.subplots(figsize=(2.6, 2.2))
        for m in METH:
            band(ax2, *xy(m, 'n', k), COL[m])
        ax2.set_xscale('log')
        ax2.set_xticks(DEPTHS)
        ax2.set_xticklabels(DL)
        ax2.set_xlabel('deduplicated molecules per sample\n(subsampled)')
        ax2.set_ylabel(f'annotated junctions (>= {k} molecule{"s" if k > 1 else ""})')
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
        save(fig2, f'p5_annotated_junctions_vs_depth_min{k}_{tag}')
        if (
            k == 3 and tag == 'native'
        ):  # main figure copy: the annotated per-sample panel beside the all-junction one
            MAIN = f'{WORK}/preprint_figures_native/main_figure_panels'
            os.makedirs(f'{MAIN}/values', exist_ok=True)
            fig3, ax3 = plt.subplots(figsize=(2.6, 2.2))
            nums = {}
            for m in METH:
                x, y = xy(m, 'n', k)
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
            ax3.set_ylabel(f'unique annotated junctions\n(>= {k} molecules per junction)')
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
            fig3.tight_layout(pad=0.4)
            fig3.savefig(
                f'{MAIN}/p5_annotated_junctions_vs_depth_min{k}_{tag}.svg', bbox_inches='tight', pad_inches=0.03
            )
            fig3.savefig(
                f'{MAIN}/p5_annotated_junctions_vs_depth_min{k}_{tag}.png',
                dpi=200,
                bbox_inches='tight',
                pad_inches=0.03,
            )
            plt.close(fig3)
            sv = open(f'{MAIN}/p5_annotated_junctions_vs_depth_min{k}_{tag}.svg').read()
            open(f'{MAIN}/p5_annotated_junctions_vs_depth_min{k}_{tag}.svg', 'w').write(
                re.sub(r"font-family:\s*[^;\"]+", 'font-family: Arial', sv)
            )
            json.dump(
                {
                    'set': tag,
                    'min_molecules': k,
                    'annotated_introns': len(ANN),
                    'source': 'per-sample junction BEDs (extract_junctions.sh), UMI-deduplicated MAPQ-255 molecules, annotated introns only (Ensembl 113), rarefaction by binomial thinning',
                    'per_method': nums,
                },
                open(f'{MAIN}/p5_annotated_junctions_vs_depth_min{k}_{tag}.json', 'w'),
                indent=1,
            )
            with open(f'{MAIN}/values/p5_annotated_junctions_vs_depth_min{k}_{tag}.tsv', 'w') as fv:
                fv.write(
                    'method\tsample\tdedup_molecules_total\tannotated_junctions_total\tdepth\tannotated_junctions_expected\tall_junctions_expected\n'
                )
                [
                    fv.write(
                        f"{m}\t{meta[a]['sample_name']}\t{R[m][a]['N']}\t{R[m][a]['jn']}\t{d}\t{R[m][a][('n', k, d)]:.1f}\t{R[m][a][('a', k, d)]:.1f}\n"
                    )
                    for m in METH
                    for a in R[m]
                    for d in DEPTHS
                    if not np.isnan(R[m][a][('n', k, d)])
                ]
    fig.suptitle(
        f'annotated splice junctions vs depth, {tag} set: line = mean over samples, ribbon = 95% t-interval, dashed = all junctions; '
        + ', '.join(f'{m} {COL[m]}' for m in METH),
        fontsize=7.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(f'{OUT}/overview_annotated_junctions_{tag}.png', dpi=170, bbox_inches='tight')
    plt.close(fig)
    S += [
        f'## {tag} set',
        '',
        'Annotated share per sample (median): '
        + '; '.join(
            f'{m} {share[m][0]:.0f}% of distinct junctions, {share[m][1]:.0f}% of junction molecules' for m in METH
        ),
        '',
        '| min molecules | depth | '
        + ' | '.join(METH)
        + ' | BOBseq / DRUG-seq | BOBseq / prime-seq | same ratios, all junctions |',
        '|---|---|' + '---|' * (len(METH) + 3),
    ]
    for k in K:
        for d in (100000, 250000, 500000):
            mu = {}
            mua = {}
            for m in METH:
                v = [R[m][a][('n', k, d)] for a in R[m] if not np.isnan(R[m][a][('n', k, d)])]
                va = [R[m][a][('a', k, d)] for a in R[m] if not np.isnan(R[m][a][('a', k, d)])]
                mu[m] = (np.mean(v), len(v)) if v else (np.nan, 0)
                mua[m] = np.mean(va) if va else np.nan
            fold = lambda q, o: (
                'na'
                if np.isnan(q['BOBseq'] if not isinstance(q['BOBseq'], tuple) else q['BOBseq'][0])
                or np.isnan(q[o] if not isinstance(q[o], tuple) else q[o][0])
                else f"{(q['BOBseq'][0] if isinstance(q['BOBseq'], tuple) else q['BOBseq']) / (q[o][0] if isinstance(q[o], tuple) else q[o]):.1f}x"
            )
            S.append(
                f'| {k} | {d:,} | '
                + ' | '.join('na' if not mu[m][1] else f'{mu[m][0]:,.0f} (n={mu[m][1]})' for m in METH)
                + f" | {fold(mu, 'DRUG-seq')} | {fold(mu, 'prime-seq')} | {fold(mua, 'DRUG-seq')} / {fold(mua, 'prime-seq')} |"
            )
    S.append('')
VAL.close()
open(f'{OUT}/SUMMARY_annotated_junctions.md', 'w').write('\n'.join(S))
print('DONE')

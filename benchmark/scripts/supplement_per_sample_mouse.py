#!/usr/bin/env python3
"""Per-sample supplement, the three mouse wells of the 24-plex (RAW 264.7 controls 1-3). Same panels, rules and axis ranges as the
kept human per-sample panels (supplement_per_sample.py), species = mouse (GRCm39, Ensembl 113), in their own directory: another cell line and genome,
so they are not drawn beside the human wells. No junction panel (not needed). Inputs from mouse3_build_wells.py, mouse3_tables.py and mouse3_native_bams.py.
Panels: S01 read fate, S02 filtered reads after de-duplication, S03 duplicate rate, S06 composition, S07 intronic %, S08 read length per mate, S09 insert
length, S10 aligned bases, S12 5'-3' balance (Picard), S14 coverage vs distance from the poly(A) site."""

from settings import WORK, PROJECT, COVERAGE
import os, re, sys, csv, json, numpy as np, matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, H)
from palette import CLS

P = WORK
PR = PROJECT
A = COVERAGE
ARM = f'{P}/benchmark_uniform/bob57_24plex_pe_native'
OUT = f'{P}/preprint_figures_native/supplement_per_sample_mouse'
V = f'{OUT}/values'
os.makedirs(V, exist_ok=True)
os.makedirs(f'{OUT}/superseded', exist_ok=True)
MET = f'{P}/preprint_figures_native/coverage_architecture/picard_metrics/native_mouse'
TAG = 'native_mouse'
WELLS = ['TAAGACG', 'GTCCATC', 'GATGACAGTAT']
LAB = {'TAAGACG': 'RAW_control_1', 'GTCCATC': 'RAW_control_2', 'GATGACAGTAT': 'RAW_control_3'}
SHORT = {w: LAB[w].replace('RAW_control_', 'RAW ') for w in WELLS}
MC = '#7A4FA3'
RIBO = re.compile(r'^(RP[LS]|RPLP|MRP[LS])[^K]*$', re.I)
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
    fig.tight_layout(pad=0.4)
    fig.savefig(f'{OUT}/{name}.svg', bbox_inches='tight', pad_inches=0.03)
    fig.savefig(f'{OUT}/{name}.png', dpi=200, bbox_inches='tight', pad_inches=0.03)
    plt.close(fig)
    s = open(f'{OUT}/{name}.svg').read()
    s = re.sub(r"font-family:\s*[^;\"]+", 'font-family: Arial', s)
    open(f'{OUT}/{name}.svg', 'w').write(s)
    print('wrote', name, flush=True)


LS = ['-', (0, (4, 2)), (0, (1, 1.5))]


def bars(ax, vals, ylabel, log=False):
    x = np.arange(3)
    ax.bar(x, [vals.get(w, np.nan) for w in WELLS], color=MC, width=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([SHORT[w] for w in WELLS])
    ax.set_ylabel(ylabel)
    ax.set_xlim(-0.7, 2.7)
    if log:
        ax.set_yscale('log')


def wlegend(ax, loc):
    ax.legend(
        handles=[plt.Line2D([], [], color=MC, lw=0.9, ls=LS[i]) for i in range(3)],
        labels=[SHORT[w] + ' (mouse)' for w in WELLS],
        loc=loc,
        frameon=False,
        fontsize=5.5,
        handlelength=2.2,
        borderaxespad=0.3,
        labelspacing=0.25,
    )


BW = (1.7, 2.3)
# ---- core numbers ----
inp = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(f'{ARM}/input_counts.tsv') if l.strip()}
for l in open(f'{ARM}/mouse3/n_umi_dropped_per_well.tsv'):
    w, d = l.rstrip('\n').split('\t')
    inp[w] -= int(d)
pm = {
    r[0]: (int(r[1]), int(r[2])) for r in (l.rstrip('\n').split('\t') for l in open(f'{ARM}/mouse3/perwell_mapped.tsv'))
}
rpw = {
    r[0]: int(r[1])
    for r in (l.rstrip('\n').split('\t') for l in open(f'{ARM}/stats_E_U_mouse3/basic/reads_per_well.tsv'))
}
pb = {r['well']: r for r in csv.DictReader(open(f'{ARM}/stats_E_U_mouse3/basic/perwell_B.tsv'), delimiter='\t')}
core = {
    w: dict(
        reads_in=inp[w],
        mapped=pm[w][0],
        unique=pm[w][1],
        filtered=rpw[w],
        molecules_UCI=float(pb[w]['B_corr']),
        genes=int(pb[w]['genes_B_noRP']),
        duplicate_rate_pct=100 * (1 - float(pb[w]['B_corr']) / rpw[w]),
    )
    for w in WELLS
}
for w in WELLS:
    assert (
        core[w]['reads_in']
        >= core[w]['mapped']
        >= core[w]['unique']
        >= core[w]['filtered']
        >= core[w]['molecules_UCI'] * 0.5
    ), (w, core[w])
fig, ax = plt.subplots(figsize=(2.8, 2.3))
steps = ['reads_in', 'mapped', 'unique', 'filtered']
for i, w in enumerate(WELLS):
    ax.plot(range(4), [100 * core[w][k] / core[w]['reads_in'] for k in steps], color=MC, lw=0.9, ls=LS[i])
ax.set_xticks(range(4))
ax.set_xticklabels(['reads\ninto STAR', 'mapped', 'uniquely\nmapped', 'filtered'], fontsize=6)
ax.set_ylim(0, 100)
ax.set_ylabel('% of reads into STAR (per well)')
wlegend(ax, 'lower left')
save(fig, f'S01_read_fate_{TAG}')
fig, ax = plt.subplots(figsize=BW)
bars(
    ax,
    {w: core[w]['molecules_UCI'] for w in WELLS},
    'filtered reads per sample\nafter de-duplication (molecules, UMI)',
    log=True,
)
ax.set_ylim(7e4, 2e6)
save(fig, f'S02_reads_after_dedup_{TAG}')
fig, ax = plt.subplots(figsize=BW)
bars(ax, {w: core[w]['duplicate_rate_pct'] for w in WELLS}, 'duplicate rate (%)\nas sequenced')
ax.set_ylim(0, 100)
save(fig, f'S03_duplicate_rate_{TAG}')
# ---- composition ----
comp = json.load(open(f'{PR}/results/supplement_21/composition_rules_mouse3_native.json'))
CK = [k for k, _, _ in CLS]
compv = {}
fig, ax = plt.subplots(figsize=BW)
bottom = np.zeros(3)
x = np.arange(3)
for k, l, c in CLS:
    v = np.array([100 * comp[LAB[w]].get(k, 0) / sum(comp[LAB[w]].get(kk, 0) for kk in CK) for w in WELLS])
    compv[k] = v
    ax.bar(x, v, 0.6, bottom=bottom, color=c, lw=0)
    bottom += v
ax.set_xticks(x)
ax.set_xticklabels([SHORT[w] for w in WELLS])
ax.set_ylabel('% of mapped reads')
ax.set_ylim(0, 100)
ax.set_xlim(-0.7, 2.7)
save(fig, f'S06_composition_{TAG}')
fig, ax = plt.subplots(figsize=BW)
bars(ax, {w: compv['intronic'][i] for i, w in enumerate(WELLS)}, 'intronic reads\n(% of mapped)')
ax.set_ylim(0, None)
save(fig, f'S07_intronic_{TAG}')
# ---- read length, aligned bases, coverage vs poly(A) (positions npz) ----
rl = {}
ab = {}
prof = {}
ngen = {}
W3 = 3000
MINR = 20


def sm(y, k=25):
    y = np.asarray(y, float)
    v = np.where(np.isnan(y), 0, y)
    n = (~np.isnan(y)).astype(float)
    ker = np.ones(k)
    return np.convolve(v, ker, 'same') / np.maximum(np.convolve(n, ker, 'same'), 1)


for k_, w in enumerate(WELLS, 1):
    z = np.load(f'{A}/positions_canonical/bobseq-raw-control-{k_}.npz')
    al = z['alen'].astype(np.int64)
    mate = z['mate']
    nf = int((mate == 1).sum())
    rl[w] = {'m1': float(np.median(al[mate == 1])), 'm2': float(np.median(al[mate == 2]))}
    ab[w] = {'total': int(al.sum()), 'per_fragment': al.sum() / max(nf, 1), 'reads': len(al)}
    genes_ = list(z['genes'])
    tl = z['tx_len']
    g = z['gene']
    T_ = z['t']
    o = np.argsort(g, kind='stable')
    gs = g[o]
    st = np.flatnonzero(np.r_[True, gs[1:] != gs[:-1]])
    en = np.r_[st[1:], len(gs)]
    Uw = []
    for a_, b_ in zip(st, en):
        gi = int(gs[a_])
        ii = o[a_:b_]
        L = int(tl[gi])
        if L < 1000 or RIBO.match(genes_[gi]) or len(ii) < MINR:
            continue
        t = T_[ii].astype(int)
        e = np.minimum(t + al[ii].astype(int), L)
        d = np.zeros(L + 1)
        np.add.at(d, t, 1)
        np.add.at(d, e, -1)
        c = np.cumsum(d)[:L]
        zz = c / c.mean()
        row = np.full(W3, np.nan)
        wl = min(W3, L)
        row[:wl] = zz[::-1][:wl]
        Uw.append(row)
    prof[w] = sm(np.nanmean(np.array(Uw), 0))
    ngen[w] = len(Uw)
fig, ax = plt.subplots(figsize=BW)
x = np.arange(3)
ax.bar(x - 0.17, [rl[w]['m1'] for w in WELLS], 0.34, color=MC)
ax.bar(x + 0.17, [rl[w]['m2'] for w in WELLS], 0.34, color=MC, alpha=0.45)
ax.set_xticks(x)
ax.set_xticklabels([SHORT[w] for w in WELLS])
ax.set_ylabel('median aligned length (nt)\ndark read 2, light read 1')
ax.set_xlim(-0.7, 2.7)
ax.set_ylim(0, 160)
save(fig, f'S08_read_length_{TAG}')
fig, ax = plt.subplots(figsize=BW)
bars(ax, {w: ab[w]['total'] / 1e9 for w in WELLS}, 'aligned bases (Gb)\nunique reads, as sequenced', log=True)
save(fig, f'S10_aligned_bases_{TAG}')
fig, ax = plt.subplots(figsize=(2.8, 2.3))
for i, w in enumerate(WELLS):
    ax.plot(np.arange(W3), prof[w], color=MC, lw=0.9, ls=LS[i])
ax.axhline(1, color='#9a9a9a', lw=0.6, ls=(0, (2, 2)))
ax.set_xlim(W3, 0)
ax.set_ylim(0, 8.5)
ax.set_xlabel('distance from the poly(A) site (nt)')
ax.set_ylabel('coverage / gene mean\n(genes >= 1 kb, RP excluded)')
wlegend(ax, 'upper left')
save(fig, f'superseded/S14_mouse_only_coverage_vs_nt_from_polyA_{TAG}')
# ---- insert length ----
ins = {w: json.load(open(f'{PR}/results/library_metrics/insert_per_well/{LAB[w]}.json')) for w in WELLS}
fig, ax = plt.subplots(figsize=(2.8, 2.3))
xs = np.arange(1, 601)
for i, w in enumerate(WELLS):
    h = np.array(ins[w]['hist_1nt_to_3000'], float)
    ax.plot(xs, np.convolve(100 * h / h.sum(), np.ones(5) / 5, 'same')[1:601], color=MC, lw=0.9, ls=LS[i])
ax.set_xlabel('insert length on the mature transcript (nt)')
ax.set_ylabel('% of deduplicated pairs per nt')
ax.set_xlim(0, 600)
ax.set_ylim(0, None)
wlegend(ax, 'upper right')
save(fig, f'S09_insert_length_{TAG}')


# ---- Picard balance ----
def parse(path):
    hist = []
    inh = False
    for l in open(path):
        r = l.rstrip('\n').split('\t')
        if r and r[0].startswith('## HISTOGRAM'):
            inh = True
            continue
        if inh and len(r) >= 2 and r[0].isdigit():
            hist.append(float(r[1]))
    return np.array(hist)


p = np.arange(101) / 100
pic = {}
for w in WELLS:
    h = parse(f'{MET}/{LAB[w]}.RNA_Metrics.txt')
    assert len(h) == 101 and h.sum() > 0, w
    pic[w] = {'profile': h, 'balance': float(2 * (p * h).sum() / h.sum())}
fig, ax = plt.subplots(figsize=BW)
bars(ax, {w: pic[w]['balance'] for w in WELLS}, "5' to 3' balance\n(2 x centroid)")
ax.axhline(1, color='#7f7f7f', lw=0.7, ls=(0, (4, 3)))
ax.set_ylim(0, 2)
save(fig, f'S12_picard_balance_{TAG}')
# ---- values + README ----
with open(f'{V}/per_well_{TAG}.tsv', 'w') as f:
    f.write(
        'bobcode\tsample\tspecies\treads_in\tmapped\tunique\tfiltered_reads\tmolecules_UCI\tgenes_no_ribosomal_protein\tduplicate_rate_pct\t'
        + '\t'.join(f'pct_{k}' for k in CK)
        + '\trRNA_records_on_human_rdna_loci\tmedian_aligned_len_mate1\tmedian_aligned_len_mate2\ttotal_aligned_bases\tbases_per_fragment\tinsert_median\tinsert_p10\tinsert_p90\tpicard_balance\tpolyA_profile_genes\n'
    )
    for i, w in enumerate(WELLS):
        c = core[w]
        f.write(
            '\t'.join(
                [
                    w,
                    LAB[w],
                    'mouse',
                    str(c['reads_in']),
                    str(c['mapped']),
                    str(c['unique']),
                    str(c['filtered']),
                    f"{c['molecules_UCI']:.0f}",
                    str(c['genes']),
                    f"{c['duplicate_rate_pct']:.2f}",
                ]
                + [f'{compv[k][i]:.2f}' for k in CK]
                + [
                    str(comp[LAB[w]]['rRNA'] - comp[LAB[w]]['rRNA_on_mouse_contigs']),
                    f"{rl[w]['m1']:.0f}",
                    f"{rl[w]['m2']:.0f}",
                    str(ab[w]['total']),
                    f"{ab[w]['per_fragment']:.1f}",
                    f"{ins[w]['insert_median']:.0f}",
                    f"{ins[w]['insert_quantiles']['10']:.0f}",
                    f"{ins[w]['insert_quantiles']['90']:.0f}",
                    f"{pic[w]['balance']:.3f}",
                    str(ngen[w]),
                ]
            )
            + '\n'
        )
with open(f'{V}/polyA_profiles_{TAG}.tsv', 'w') as f:
    f.write('sample\tgenes\t' + '\t'.join(f'nt{d}' for d in range(W3)) + '\n')
    [f.write(f'{LAB[w]}\t{ngen[w]}\t' + '\t'.join(f'{v:.4f}' for v in prof[w]) + '\n') for w in WELLS]
open(f'{OUT}/README.md', 'w').write(
    f"# Per-sample supplement, the three mouse wells of the 24-plex ({__import__('datetime').date.today()})\n\nRAW 264.7 controls 1 to 3, the species controls of the run, through the same path and rules as the human wells with the mouse half of the combined reference (GRCm39, Ensembl 113). Kept apart from the human panels: another cell line and genome. Same panels and axis ranges as the kept human per-sample panels, no junction panel.\n\n"
    "Species-specific points: filtered reads = uniquely mapped, exon of a protein-coding gene, on a mouse contig. Read composition = the classification rules on the mouse-contig records, plus the well's records on the human rDNA loci counted as rRNA: GRCm39 carries almost no rDNA, so the aligner places mouse rRNA on the human rDNA copies of the combined reference (column rRNA_records_on_human_rdna_loci); the few other human-contig records are left out, as mouse-contig records are for the human wells. "
    "Coverage vs distance from the poly(A) site: canonical mouse transcripts of at least 1 kb, ribosomal-protein genes excluded, genes with at least 20 reads in the well (the human panel uses the benchmark gene set shared by the three methods, which has no mouse counterpart). 5'-3' balance: Picard CollectRnaSeqMetrics on a mouse refFlat built with the rule of the human one.\n"
)
print('MOUSE SUPPLEMENT DONE')

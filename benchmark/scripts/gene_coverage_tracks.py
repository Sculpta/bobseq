#!/usr/bin/env python3
"""Browser-style coverage figures for example genes: one filled coverage track per method (method colour, label top
left, the peak depth as the only y tick on the left), the canonical transcript below (exons as boxes on the intron line, gene name left, 3' end marked), hg38 positions
under it. The gene is drawn 5' to 3' from left to right, so minus-strand genes run with decreasing coordinates. Two versions per gene: genomic (introns to scale) and
exon-only (introns removed, exons to scale; ticks give the genomic position of exon starts). Coverage = read depth per base (aligned blocks only, spliced gaps are not
covered) on the pooled depth-matched BAMs (ucsc_pooled: exactly 8.0 M uniquely mapped fragments per method), so peak values are comparable between methods.
    gene_coverage_tracks.py RELA MCM7 ...   ->  preprint_figures_native/gene_coverage_tracks/<gene>_(with_introns|exons_only).(svg|png) + values/<gene>_depth.tsv (GCT_OUT overrides the directory for candidate sheets)
"""

from settings import WORK

# GCT_OUT = alternative output directory (candidate sheets)
import os, re, sys, numpy as np, pysam, matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, H)
import coverage_lib as cl
from palette import COL, METH

P = WORK
B = f'{WORK}/benchmark_uniform/ucsc_pooled'
OUT = os.environ.get('GCT_OUT', f'{P}/preprint_figures_native/gene_coverage_tracks')
os.makedirs(f'{OUT}/values', exist_ok=True)
TAG = {'DRUG-seq': 'drugseq', 'prime-seq': 'primeseq', 'BOBseq': 'bobseq'}
plt.rcParams.update(
    {
        'svg.fonttype': 'none',
        'font.family': ['Arial', 'Nimbus Sans', 'DejaVu Sans'],
        'font.size': 8,
        'axes.linewidth': 0.7,
        'text.color': '#231F20',
        'axes.edgecolor': '#231F20',
        'xtick.color': '#231F20',
        'ytick.color': '#231F20',
    }
)
GENECOL = '#2F6DB5'
TX, _ = cl.build()
PAD = 0.03


def depth(tag, chrom, s0, e0):
    d = np.zeros(e0 - s0 + 1, np.int64)
    for a in pysam.AlignmentFile(f'{B}/{tag}_pooled_matched.bam', 'rb').fetch(chrom, s0, e0):
        for bs, be in a.get_blocks():
            d[min(max(bs - s0, 0), e0 - s0)] += 1
            d[min(max(be - s0, 0), e0 - s0)] -= 1
    return np.cumsum(d)[: e0 - s0]


def save(fig, name):
    for ext, kw in (('svg', {}), ('png', {'dpi': 220})):
        fig.savefig(f'{OUT}/{name}.{ext}', bbox_inches='tight', pad_inches=0.04, **kw)
    s = open(f'{OUT}/{name}.svg').read()
    open(f'{OUT}/{name}.svg', 'w').write(re.sub(r"font-family:\s*[^;\"]+", 'font-family: Arial', s))
    plt.close(fig)
    print('wrote', name, flush=True)


def draw(gene, x, tracks, exon_x, ticks, ticklabels, xlim, name, note, joined=False):
    fig, axs = plt.subplots(
        len(METH) + 1,
        1,
        figsize=(6.8, 0.95 * len(METH) + 0.75),
        gridspec_kw={'height_ratios': [1] * len(METH) + [0.62], 'hspace': 0.12},
        sharex=True,
    )
    for ax, m in zip(axs, METH):
        y = tracks[m]
        pk = int(y.max())
        ax.fill_between(x, 0, y, step='mid', color=COL[m], lw=0)
        ax.set_ylim(0, max(pk, 1) * 1.3)
        ax.set_yticks([pk])
        ax.set_yticklabels([f'{pk:,}'], fontsize=9)
        for sp in ('top', 'right', 'bottom'):
            ax.spines[sp].set_visible(False)
        ax.axhline(0, color=COL[m], lw=0.4, alpha=0.5)
        ax.tick_params(axis='x', bottom=False, labelbottom=False)
        ax.tick_params(axis='y', length=4)
        ax.text(
            0.015,
            0.97,
            m,
            color=COL[m],
            fontsize=9.5,
            transform=ax.transAxes,
            va='top',
            ha='left',
            zorder=6,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, pad=1.2),
        )
    ax = axs[-1]
    ax.set_ylim(-1.35, 1)
    ax.plot([exon_x[0][0], exon_x[-1][1]], [0, 0], color=GENECOL, lw=0.9, zorder=1)
    for a, b in exon_x:
        ax.add_patch(plt.Rectangle((min(a, b), -0.42), abs(b - a), 0.84, color=GENECOL, lw=0, zorder=2))
    if joined:
        [ax.plot([a, a], [-0.42, 0.42], color='white', lw=1.1, zorder=3, solid_capstyle='butt') for a, _ in exon_x[1:]]
    for sp in ('top', 'right', 'left'):
        ax.spines[sp].set_visible(False)
    ax.spines['bottom'].set_position(('data', -1.0))
    ax.set_yticks([])
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticklabels, fontsize=7)
    ax.tick_params(axis='x', length=3)
    ax.text(-0.012, 0.69, gene, color=GENECOL, fontsize=9.5, transform=ax.transAxes, ha='right', va='center')
    ax.text(1.012, 0.69, '3p', color=GENECOL, fontsize=9.5, transform=ax.transAxes, ha='left', va='center')
    ax.text(-0.012, 0.13, 'hg38', fontsize=8.5, transform=ax.transAxes, ha='right', va='center')
    ax.text(1.0, -0.5, note, fontsize=6.5, color='#666', transform=ax.transAxes, ha='right', va='top')
    axs[0].set_xlim(*xlim)
    save(fig, name)


for gene in sys.argv[1:]:
    L, tid, c, strand, starts, ends = TX[gene]
    chrom = 'chr' + c
    s0, e0 = starts[0], ends[-1]
    D = {m: depth(TAG[m], chrom, s0, e0) for m in METH}
    rev = strand == '-'
    span = e0 - s0
    # ---- genomic version: introns to scale, 5' left ----
    gx = np.arange(s0, e0) + 1
    step = [v for v in (500, 1000, 2000, 5000, 10000, 20000) if span / v <= 6][0]
    tk = [t for t in range((s0 // step + 1) * step, e0, step)]
    xlim = (e0 + PAD * span, s0 - PAD * span) if rev else (s0 - PAD * span, e0 + PAD * span)
    draw(
        gene,
        gx,
        D,
        [(a + 1, b) for a, b in zip(starts, ends)],
        tk,
        [f'{t:,}' for t in tk],
        xlim,
        f'{gene}_with_introns',
        f'{chrom}, {tid}, {len(starts)} exons, {strand} strand; pooled benchmark samples, 8.0 M uniquely mapped fragments per method',
    )
    # ---- exon-only version: exons concatenated 5' to 3' ----
    ex = list(zip(starts, ends))
    ex = ex[::-1] if rev else ex
    T = {
        m: np.concatenate([(D[m][a - s0 : b - s0][::-1] if rev else D[m][a - s0 : b - s0]) for a, b in ex])
        for m in METH
    }
    off = np.r_[0, np.cumsum([b - a for a, b in ex])]
    exon_x = [(off[i], off[i + 1]) for i in range(len(ex))]
    sel = [0]
    [sel.append(i) for i in range(1, len(ex)) if off[i] - off[sel[-1]] >= L / 7]
    draw(
        gene,
        np.arange(L) + 0.5,
        T,
        exon_x,
        [off[i] for i in sel],
        [f'{(ex[i][1] if rev else ex[i][0] + 1):,}' for i in sel],
        (-PAD * L, L * (1 + PAD)),
        f'{gene}_exons_only',
        f'{chrom}, {tid}, introns removed ({len(ex)} exons, {L:,} nt, to scale); ticks = hg38 position of the exon start; 8.0 M fragments per method',
        joined=True,
    )
    with open(f'{OUT}/values/{gene}_depth.tsv', 'w') as f:
        f.write('chrom\tpos_hg38\tin_exon\t' + '\t'.join(METH) + '\n')
        inex = np.zeros(span, bool)
        for a, b in zip(starts, ends):
            inex[a - s0 : b - s0] = True
        for i in range(span):
            f.write(f'{chrom}\t{s0 + i + 1}\t{int(inex[i])}\t' + '\t'.join(str(int(D[m][i])) for m in METH) + '\n')
    gm = f'{OUT}/values/gene_models.tsv'
    rows_ = (
        [l for l in open(gm)]
        if os.path.exists(gm)
        else ['gene\tchrom\tstrand\ttranscript\texon_starts_1based\texon_ends\n']
    )
    rows_ = [l for l in rows_ if not l.startswith(gene + '\t')] + [
        f"{gene}\t{chrom}\t{strand}\t{tid}\t{','.join(str(a + 1) for a in starts)}\t{','.join(str(b) for b in ends)}\n"
    ]
    open(gm, 'w').write(''.join(rows_))  # gene model for replotting from the values alone
    print(
        gene,
        chrom,
        f'{s0 + 1:,}-{e0:,}',
        strand,
        tid,
        '| peak depth',
        {m: int(D[m].max()) for m in METH},
        '| exonic peak',
        {m: int(T[m].max()) for m in METH},
        flush=True,
    )

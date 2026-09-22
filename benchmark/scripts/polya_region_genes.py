#!/usr/bin/env python3
"""Per-gene view in nucleotides from the poly(A) site: read 5' ends of every sample along the mature transcript of ACTB, GAPDH (and two
more), plotted as distance from the annotated 3' end. Shows where each method's reads sit relative to the poly(A) site."""
import sys, os, numpy as np, textwrap, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
def _overlap_check(fig, name):
    """no legend / caption / axes / title boxes may intersect (figure-fraction bboxes after layout)"""
    fig.canvas.draw(); r = fig.canvas.get_renderer(); inv = fig.transFigure.inverted()
    boxes = [('axes:' + (a.get_title() or str(i))[:30], a.get_tightbbox(r).transformed(inv)) for i, a in enumerate(fig.axes)]
    for t in fig.texts: boxes.append(('text:' + t.get_text()[:30], t.get_window_extent(r).transformed(inv)))
    for l in fig.legends: boxes.append(('legend', l.get_window_extent(r).transformed(inv)))
    bad = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[j][1]
            if a.x0 < b.x1 - 0.003 and b.x0 < a.x1 - 0.003 and a.y0 < b.y1 - 0.003 and b.y0 < a.y1 - 0.003: bad.append((boxes[i][0], boxes[j][0]))
    if bad: print(f'LAYOUT WARNING {name}: overlapping elements: {bad[:4]}')
    else: print(f'layout ok {name}')
plt.rcParams.update({'svg.fonttype': 'none', 'font.size': 9})
from settings import WORK, COVERAGE
A = COVERAGE; P = WORK
SET = os.environ.get('BM_COVSET', 'native'); SFX = '' if SET == 'native' else f'_{SET}'; FIG = 'preprint_figures_native' if SET == 'native' else f'preprint_figures_{SET}'
BOBNOTE = 'BOBseq: both mates' if SET == 'native' else 'every method one 50-nt single read per fragment'; READNOTE = ('BOBseq both mates of the 2 x 150 pairs; competitors one native-length read per fragment' if SET == 'native' else 'read-length matched set: every method one 50-nt single read per fragment (BOBseq R2 truncated to 50 nt)')
POS = f'{A}/positions_canonical{SFX}'; OUT = f'{P}/{FIG}/coverage_architecture'; os.makedirs(OUT, exist_ok=True)
from palette import COL, METH
S = [l.rstrip('\n').split('\t') for l in open(f'{A}/samples{SFX}.tsv') if l.strip()]
GENES = ['ACTB', 'GAPDH', 'EEF1A1', 'RPLP0']; BIN = 50; XMAX = 2000
z0 = np.load(f'{POS}/{S[0][0]}.npz'); genes = list(z0['genes']); tx_len = z0['tx_len']; gid = {g: genes.index(g) for g in GENES}
fig, ax = plt.subplots(2, 2, figsize=(14, 9)); rows = ['| gene | sample | method | reads | % within 300 nt of the poly(A) site | % within 1000 nt |', '|---|---|---|---|---|---|']
for a, g in zip(ax.flat, GENES):
    L = int(tx_len[gid[g]]); edges = np.arange(0, min(XMAX, L) + BIN, BIN)
    for lab, bam, m in S:
        z = np.load(f'{POS}/{lab}.npz'); sel = z['gene'] == gid[g]
        if m == 'BOBseq': pass  # both mates
        d3 = (z['L'][sel] - 1 - z['t'][sel]); n = sel.sum()
        if n < 50: continue
        h, _ = np.histogram(d3, bins=edges); a.plot(edges[:-1] + BIN / 2, 100 * h / n, color=COL[m], lw=1.6 if m != 'BOBseq' else 1.0, alpha=.9 if m != 'BOBseq' else .7)
        rows.append(f'| {g} | {lab} | {m} | {n:,} | {100*np.mean(d3 <= 300):.0f}% | {100*np.mean(d3 <= 1000):.0f}% |')
    a.set_title(f'{g}: mature transcript {L:,} nt, read 5\' ends by distance from the poly(A) site ({BIN}-nt bins)', fontsize=9); a.set_xlabel("distance from the annotated 3' end / poly(A) site (nt)"); a.set_ylabel(f"% of the gene's read starts per {BIN}-nt bin"); a.grid(ls=':', lw=.5); a.set_xlim(0, min(XMAX, L))
    a.invert_xaxis()
fig.suptitle("Poly(A)-proximal bias, per gene: where do the reads of each sample sit relative to the poly(A) site?", fontsize=11)
cap = ("TESTED: is the poly(A)-region bias visible in the browser (ACTB, GAPDH) captured quantitatively? For four housekeeping genes, the 5' end of every unique read of every sample is placed on the mature transcript "
       "(exons of the Ensembl canonical protein-coding transcript) and binned by its distance from the annotated 3' end (0 = poly(A) site, x axis reversed so the transcript reads left to right as in the browser for a minus-strand gene). "
       "Each line is one sample (" + (f"{sum(1 for s in S if s[2] == 'BOBseq')} BOBseq " + ('both mates' if SET == 'native' else '50-nt R2')) + ", 3 per competitor), normalized to that sample's read count on the gene, so heights compare as fractions. A 3'-end method piles up within a few hundred nt of the poly(A) "
       "site; BOBseq shows what fraction of its reads sit there versus along the body. The table in the summary lists the share within 300 and 1,000 nt per sample.")
W, Hh = fig.get_size_inches(); wrapped = textwrap.fill(cap, int(W * 13.5)); cap_h = ((wrapped.count('\n') + 1) * 0.135 + 0.12) / Hh
fig.text(0.01, 0.006, wrapped, fontsize=7.6, va='bottom', ha='left', color='#333')
fig.legend(handles=[Line2D([], [], color=COL[m], lw=2.5) for m in METH], labels=METH, loc='lower center', ncol=4, frameon=False, fontsize=9, bbox_to_anchor=(0.5, cap_h + 0.01))
fig.tight_layout(rect=(0, cap_h + 0.07, 1, 0.95)); _overlap_check(fig, 'J_polyA_region_genes'); fig.savefig(f'{OUT}/J_polyA_region_genes.png', dpi=140); fig.savefig(f'{OUT}/J_polyA_region_genes.svg')
open(f'{OUT}/J_polyA_region_genes.md', 'w').write('## J. Poly(A)-proximal share per sample for four housekeeping genes\n' + '\n'.join(rows) + '\n'); print('\n'.join(rows[:2] + [r for r in rows if '| ACTB |' in r or '| GAPDH |' in r]))

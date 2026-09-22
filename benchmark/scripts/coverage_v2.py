#!/usr/bin/env python3
"""Base-coverage view of the transcript (what the browser shows), per gene at matched depth. Sections:
 V1 per-gene coverage heatmaps (aggregated browser view), V2 fraction of the transcript covered vs reads per gene, V3 peakiness (share of a
 gene's reads in its densest 300 nt), V4 aligned bases per read / fragment, V5 fraction covered by transcript length, V6 filters and exclusion
 rules (does any rule change the ranking?). Inputs: positions_canonical/*.npz (read 5' end t, aligned length alen on the canonical transcript)."""
import os, sys, json, textwrap, collections, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({'svg.fonttype': 'none', 'font.size': 9})
from settings import WORK, COVERAGE
A = COVERAGE; P = WORK; H = os.path.dirname(os.path.abspath(__file__))   # A: the position tables, H: this directory (sibling scripts)
SET = os.environ.get('BM_COVSET', 'native'); SFX = '' if SET == 'native' else f'_{SET}'; FIG = 'preprint_figures_native' if SET == 'native' else f'preprint_figures_{SET}'
BOBNOTE = 'BOBseq: both mates' if SET == 'native' else 'every method one 50-nt single read per fragment'; READNOTE = ('BOBseq both mates of the 2 x 150 pairs; competitors one native-length read per fragment' if SET == 'native' else 'read-length matched set: every method one 50-nt single read per fragment (BOBseq R2 truncated to 50 nt)')
POS = f'{A}/positions_canonical{SFX}'; OUT = f'{P}/{FIG}/coverage_architecture'; os.makedirs(OUT, exist_ok=True)
from palette import COL, METH
S = [l.rstrip('\n').split('\t') for l in open(f'{A}/samples{SFX}.tsv') if l.strip()]; rng = np.random.default_rng(11); MD = []
def md(s=''): MD.append(s)
exec(open(f'{H}/coverage_architecture.py').read().split('def _overlap_check')[1].split('\ndef save(')[0].join(['def _overlap_check', '']))
def save(fig, name, caption, handles=None, labels=None, legend=True):
    W, Hh = fig.get_size_inches(); wrapped = textwrap.fill(caption, int(W * 13.5)); n = wrapped.count('\n') + 1; cap_h = (n * 0.135 + 0.12) / Hh
    fig.text(0.01, 0.006, wrapped, fontsize=7.6, va='bottom', ha='left', color='#333')
    if legend: fig.legend(handles=handles or [Line2D([], [], color=COL[m], lw=2.5) for m in METH], labels=labels or METH, loc='lower center', ncol=len(labels or METH), frameon=False, fontsize=9, bbox_to_anchor=(0.5, cap_h + 0.01))
    fig.tight_layout(rect=(0, cap_h + (0.07 if legend else 0.01), 1, 0.95)); _overlap_check(fig, name); fig.savefig(f'{OUT}/{name}.png', dpi=140); fig.savefig(f'{OUT}/{name}.svg'); plt.close(fig); print('wrote', name)
# ---------- load, pooled per method (+ BOBseq mate-1 only as a control) ----------
POOL = {}
for lab, bam, m in S:
    z = np.load(f'{POS}/{lab}.npz')
    for k in ('gene', 't', 'L', 'alen', 'mate'): POOL.setdefault(m, {}).setdefault(k, []).append(z[k])
    if 'tx_len' not in globals(): tx_len = z['tx_len']; genes = z['genes']
for m in METH: POOL[m] = {k: np.concatenate(v) for k, v in POOL[m].items()}
sel = POOL['BOBseq']['mate'] == 1
if sel.any() and (~sel).any(): POOL['BOBseq mate 1 only'] = {k: v[sel] for k, v in POOL['BOBseq'].items()}   # paired native set only
def by_gene(d):
    g = d['gene']; o = np.argsort(g, kind='stable'); gs = g[o]; st = np.flatnonzero(np.r_[True, gs[1:] != gs[:-1]]); en = np.r_[st[1:], len(gs)]; return {int(gs[s]): o[s:e] for s, e in zip(st, en)}
IDX = {m: by_gene(POOL[m]) for m in POOL}
def gene_cov(d, ii, L):
    """per-base read coverage (depth) and start density along the mature transcript for the reads ii of one gene"""
    t = d['t'][ii].astype(int); e = np.minimum(t + d['alen'][ii].astype(int), L); diff = np.zeros(L + 1); np.add.at(diff, t, 1); np.add.at(diff, e, -1)
    return np.cumsum(diff)[:L], t
def common_genes(d0, methods=METH, lo=0, hi=10**9): return [g for g in IDX[methods[0]] if all(len(IDX[m].get(g, [])) >= d0 for m in methods) and lo <= tx_len[g] < hi]
def metrics(m, g, d0, exclude3p=0):
    idx = IDX[m][g]; ii = rng.choice(idx, min(d0, len(idx)), replace=False); L = int(tx_len[g]); cov, t = gene_cov(POOL[m], ii, L)
    Lx = L - exclude3p; covx = cov[:Lx]; frac1 = np.mean(covx > 0); frac3 = np.mean(covx >= 3)
    cs = np.cumsum(np.bincount(t, minlength=L)); w = min(300, L); peak = (cs[w - 1:] - np.r_[0, cs[:-w]]).max() / len(ii)
    return frac1, frac3, peak, cov
# ---------- V1 heatmaps ----------
def V1(d0=200, nb=50):
    common = common_genes(d0); common = sorted(common, key=lambda g: tx_len[g]); H = {m: np.zeros((len(common), nb)) for m in METH}
    for m in METH:
        for r, g in enumerate(common):
            _, _, _, cov = metrics(m, g, d0); L = len(cov); idx = np.minimum((np.arange(L) * nb) // L, nb - 1); b = np.bincount(idx, weights=cov, minlength=nb) / np.bincount(idx, minlength=nb); H[m][r] = b / b.max() if b.max() > 0 else b
    fig, ax = plt.subplots(1, len(METH), figsize=(3.75 * len(METH), 6.2), sharey=True)   # one panel per method
    for a, m in zip(ax, METH):
        a.imshow(H[m], aspect='auto', cmap='Greys', vmin=0, vmax=1, interpolation='nearest'); a.set_title(m, fontsize=10, color=COL[m]); a.set_xticks([0, nb / 2, nb - 1]); a.set_xticklabels(["5'", '50%', "3'"]); a.set_xlabel('position along the mature transcript')
    ax[0].set_ylabel(f'{len(common):,} genes, sorted by transcript length (short at top)'); yt = [0, len(common) // 2, len(common) - 1]; ax[0].set_yticks(yt); ax[0].set_yticklabels([f'{int(tx_len[common[i]]):,} nt' for i in yt])
    fig.suptitle(f'V1. Per-gene base coverage, every gene at {d0} unique reads: one row per gene, each row scaled to its own maximum (black = the gene\'s peak)', fontsize=11)
    save(fig, 'V1_coverage_heatmaps', legend=False, caption=('TESTED: what the browser shows, for all genes at once. For every gene with >= 200 unique reads in every method (same genes in every panel), exactly 200 reads are drawn and the per-base '
        'read coverage along the Ensembl canonical mature transcript is computed from each read\'s 5\' position and aligned length (' + BOBNOTE + '). Rows are genes sorted by transcript length, columns 50 bins from 5\' to 3\', '
        'each row divided by its own maximum. A 3\'-end method shows a black stripe at the right edge whose width shrinks with transcript length; end-to-end coverage fills the row. No expression weighting, no isoform averaging.'))
    md(f'## V1. Heatmaps: {len(common):,} genes at {d0} reads each, rows scaled to their own maximum; mean over genes of the fraction of bins above 20% of the gene\'s peak: ' + ', '.join(f'{m} {100*np.mean(H[m] > 0.2):.0f}%' for m in METH))
    return common
# ---------- V2 fraction of transcript covered vs reads per gene ----------
def V2():
    depths = [25, 50, 100, 200, 400, 800]; res = {m: {} for m in POOL}
    for d0 in depths:
        common = common_genes(d0, methods=list(POOL))
        for m in POOL:
            f1, f3, pk = zip(*[metrics(m, g, d0)[:3] for g in common]); res[m][d0] = dict(genes=len(common), frac1=float(np.mean(f1)), frac3=float(np.mean(f3)), peak=float(np.median(pk)))
    fig, ax = plt.subplots(1, 3, figsize=(15, 5.6)); HL = [(Line2D([], [], color=COL[m], lw=2.5), m) for m in METH] + ([(Line2D([], [], color=COL['BOBseq'], lw=2, ls='--'), 'BOBseq mate 1 only')] if 'BOBseq mate 1 only' in POOL else [])
    for m in POOL:
        st = dict(color=COL.get(m, COL['BOBseq']), lw=2, ls='-' if m in COL else '--')
        ax[0].plot(depths, [100 * res[m][d]['frac1'] for d in depths], 'o', **st); ax[1].plot(depths, [100 * res[m][d]['frac3'] for d in depths], 'o', **st); ax[2].plot(depths, [100 * res[m][d]['peak'] for d in depths], 'o', **st)
    for a, t, yl in zip(ax, ['Fraction of the transcript covered by >= 1 read, mean over genes', 'Fraction covered by >= 3 reads', "Peakiness: share of a gene's reads in its densest 300 nt (median)"], ['% of transcript bases', '% of transcript bases', '% of reads']):
        a.set_xscale('log'); a.set_xticks(depths); a.set_xticklabels(depths); a.set_xlabel('unique reads per gene (every gene subsampled to this number)'); a.set_title(t, fontsize=9); a.set_ylabel(yl); a.set_ylim(0, 100); a.grid(ls=':', lw=.5)
    fig.suptitle('V2. Base coverage per gene at matched depth (genes with >= depth reads in every method: ' + ', '.join(f'{d}: {res["BOBseq"][d]["genes"]:,}' for d in depths) + ')', fontsize=11)
    save(fig, 'V2_fraction_covered_vs_depth', handles=[h for h, _ in HL], labels=[l for _, l in HL], caption=('TESTED: does BOBseq cover the transcript, base for base, rather than only start reads in more places? Same subsampling as before (every gene to the same read count), but now each read paints its aligned length '
        'along the canonical mature transcript. Left: % of transcript bases under at least one read; middle: under at least three; right: the share of the gene\'s reads that fall in the densest 300-nt window of that gene (a 3\'-end method '
        'puts most of its reads in one window, the browser spike). ' + ('Dashed: BOBseq with mate 1 (R2, ~134 nt) only, to separate the effect of the second mate from the chemistry.' if 'BOBseq mate 1 only' in POOL else READNOTE + '.') + ''))
    md('## V2. Fraction of the canonical transcript covered and peakiness, at matched reads per gene'); md('| method | ' + ' | '.join(f'd={d} (n={res["BOBseq"][d]["genes"]:,})' for d in depths) + ' |'); md('|---|' + '---|' * len(depths))
    for m in POOL: md(f'| {m} | ' + ' | '.join(f">=1x {100*res[m][d]['frac1']:.0f}%, >=3x {100*res[m][d]['frac3']:.0f}%, densest-300 {100*res[m][d]['peak']:.0f}%" for d in depths) + ' |')
    json.dump(res, open(f'{OUT}/V2_fraction_covered_vs_depth.json', 'w'), indent=1); return res
# ---------- V3 peakiness distributions + V4 aligned bases ----------
def V3V4(d0=200):
    common = common_genes(d0); fig, ax = plt.subplots(1, 3, figsize=(15, 5.6)); dat_pk, dat_f1 = [], []
    md(f'## V3. Per-gene peakiness and coverage at {d0} reads ({len(common):,} genes)'); md("| method | share of reads in the densest 300 nt, median (IQR) | % of transcript covered >= 1x, median (IQR) | aligned nt per read, mean | aligned nt per fragment, mean |"); md('|---|---|---|---|---|')
    for m in METH:
        f1, f3, pk = zip(*[metrics(m, g, d0)[:3] for g in common]); dat_pk.append(100 * np.array(pk)); dat_f1.append(100 * np.array(f1))
        al = POOL[m]['alen']; per_read = al.mean(); per_frag = per_read * (2 if m == 'BOBseq' else 1)
        q = np.percentile(pk, [25, 50, 75]); q2 = np.percentile(f1, [25, 50, 75]); md(f'| {m} | {100*q[1]:.0f}% ({100*q[0]:.0f}-{100*q[2]:.0f}) | {100*q2[1]:.0f}% ({100*q2[0]:.0f}-{100*q2[2]:.0f}) | {per_read:.0f} | {per_frag:.0f} |')
        ax[2].bar(m, per_frag, color=COL[m]); ax[2].text(m, per_frag + 4, f'{per_frag:.0f}', ha='center', fontsize=8)
    for a, dat, t, yl in zip(ax[:2], [dat_pk, dat_f1], ["share of the gene's reads in its densest 300 nt", 'fraction of the transcript covered by >= 1 read'], ['% of reads', '% of transcript bases']):
        v = a.violinplot(dat, showmedians=True); [b.set_facecolor(COL[m]) for b, m in zip(v['bodies'], METH)]; a.set_xticks(range(1, len(METH) + 1)); a.set_xticklabels(METH); a.set_title(t + f', {len(common):,} genes at {d0} reads', fontsize=9); a.set_ylabel(yl); a.set_ylim(0, 100); a.grid(axis='y', ls=':', lw=.5)
    ax[2].set_title('aligned bases per fragment (mean; ' + BOBNOTE + ')', fontsize=9); ax[2].set_ylabel('nt'); ax[2].grid(axis='y', ls=':', lw=.5)
    fig.suptitle('V3/V4. Per-gene peakiness and covered fraction at matched depth; aligned bases per fragment', fontsize=11)
    save(fig, 'V3_peakiness_and_covered_fraction', legend=False, caption=('TESTED: is the browser spike a per-gene property (not an average)? Left: for each gene at 200 reads, the largest share of its reads that any 300-nt window holds (100% = all reads in one window). Middle: '
        '% of the gene\'s canonical transcript under at least one read. Violins over genes, bar = median. Right: mean aligned nucleotides per fragment (one read per fragment for the 3\' methods; ' + ('the sum of both mates for' if SET == 'native' else 'one 50-nt read for') + ' '
        'BOBseq), the raw reason the browser tracks differ: a BOBseq fragment paints about five times more transcript than a DRUG-seq read.'))
# ---------- V5 by transcript length ----------
def V5(d0=200):
    classes = [('0.5-1 kb', 500, 1000), ('1-2 kb', 1000, 2000), ('2-4 kb', 2000, 4000), ('4-8 kb', 4000, 8000), ('> 8 kb', 8000, 10**9)]; fig, ax = plt.subplots(1, 2, figsize=(13, 5.6)); res = {}
    md(f'## V5. Covered fraction (>= 1 read) and peakiness by transcript length, at {d0} reads per gene'); md('| method | ' + ' | '.join(c for c, _, _ in classes) + ' |'); md('|---|' + '---|' * len(classes))
    for m in METH:
        f1s, pks, ns = [], [], []
        for cname, lo, hi in classes:
            gl = common_genes(d0, lo=lo, hi=hi); ns.append(len(gl))
            if not gl: f1s.append(np.nan); pks.append(np.nan); continue
            f1, f3, pk = zip(*[metrics(m, g, d0)[:3] for g in gl]); f1s.append(100 * np.mean(f1)); pks.append(100 * np.median(pk))
        ax[0].plot(range(5), f1s, 'o-', color=COL[m], lw=2); ax[1].plot(range(5), pks, 'o-', color=COL[m], lw=2); md(f'| {m} | ' + ' | '.join(f'covered {a:.0f}%, densest-300 {b:.0f}% (n={n:,})' for a, b, n in zip(f1s, pks, ns)) + ' |')
    for a, t, yl in zip(ax, ['fraction of the transcript covered by >= 1 read (mean over genes)', "share of the gene's reads in its densest 300 nt (median)"], ['% of transcript bases', '% of reads']):
        a.set_xticks(range(5)); a.set_xticklabels([c for c, _, _ in classes]); a.set_xlabel('mature-transcript length'); a.set_title(t, fontsize=9); a.set_ylabel(yl); a.set_ylim(0, 100); a.grid(ls=':', lw=.5)
    fig.suptitle(f'V5. Short versus long transcripts in base-coverage terms, every gene at {d0} reads', fontsize=11)
    save(fig, 'V5_covered_fraction_by_length', caption=('TESTED: "short v long", now in base coverage. Genes with >= 200 reads in every method, split by canonical mature-transcript length, each gene at exactly 200 reads. Left: mean fraction of the transcript under at least '
        'one read. Right: median share of reads in the densest 300-nt window. A fixed 3\' window covers most of a 700-nt transcript and almost nothing of a 10-kb one; full-length coverage keeps its covered fraction with length.'))
# ---------- V6 filters and exclusion rules ----------
def V6(d0=200):
    md(f'## V6. Filters and exclusion rules, {d0} reads per gene: covered fraction (>= 1x) / share in the densest 300 nt, mean / median over genes'); md('| rule | genes | ' + ' | '.join(POOL) + ' |'); md('|---|---|' + '---|' * len(POOL))   # header follows the pools actually present (METH + mate-1 pool on the paired set)
    rules = [('all genes with >= 200 reads in every method', dict()), ('transcripts >= 1 kb only', dict(lo=1000)), ('transcripts >= 2 kb only', dict(lo=2000)), ('top 500 genes by BOBseq reads', dict(top=500)), ('lower half of genes by BOBseq reads', dict(bottom=True)),
             ("exclude the 3'-most 300 nt of every transcript", dict(ex=300)), ("exclude the 3'-most 1,000 nt (transcripts >= 2 kb)", dict(ex=1000, lo=2000)), ('50 reads per gene instead of 200', dict(d=50)), ('800 reads per gene', dict(d=800))]
    for name, o in rules:
        d = o.get('d', d0); gl = common_genes(d, methods=list(POOL), lo=o.get('lo', 0))
        if 'top' in o: gl = sorted(gl, key=lambda g: -len(IDX['BOBseq'][g]))[:o['top']]
        if o.get('bottom'): gl = sorted(gl, key=lambda g: len(IDX['BOBseq'][g]))[:len(gl) // 2]
        cells = []
        for m in [m for m in POOL]:
            f1, f3, pk = zip(*[metrics(m, g, d, exclude3p=o.get('ex', 0))[:3] for g in gl]); cells.append(f'{100*np.mean(f1):.0f}% / {100*np.median(pk):.0f}%')
        md(f'| {name} | {len(gl):,} | ' + ' | '.join(cells) + ' |')
for name, fn in (('V1', V1), ('V2', V2), ('V3V4', V3V4), ('V5', V5), ('V6', V6)):
    try: print('section', name); fn()
    except Exception as e:
        import traceback; traceback.print_exc(); md(f'\n**{name} failed: {e}**\n')
open(f'{OUT}/COVERAGE_V2.md', 'w').write('# Base-coverage view of the transcript (canonical transcripts, unique reads, ' + BOBNOTE + ' unless stated)\n\n' + '\n'.join(MD) + '\n'); print('V2 DONE')

#!/usr/bin/env python3
"""Coverage architecture along the mature transcript, per method, from read 5' positions (build_positions.py output).
Replaces the Picard histogram with per-gene, depth-matched, unscaled and stratified views. Every section is independent; failures are
logged and the rest continues. Usage: coverage_architecture.py <samples.tsv> <npz_dir> <out_dir>"""
import sys, os, csv, json, traceback, collections, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'svg.fonttype': 'none', 'font.size': 9})
SAMPLES, NPZ, OUT = sys.argv[1:4]; os.makedirs(OUT, exist_ok=True)
SET = os.environ.get('BM_COVSET', 'native'); SFX = '' if SET == 'native' else f'_{SET}'; FIG = 'preprint_figures_native' if SET == 'native' else f'preprint_figures_{SET}'
BOBNOTE = 'BOBseq: both mates' if SET == 'native' else 'every method one 50-nt single read per fragment'; READNOTE = ('BOBseq both mates of the 2 x 150 pairs; competitors one native-length read per fragment' if SET == 'native' else 'read-length matched set: every method one 50-nt single read per fragment (BOBseq R2 truncated to 50 nt)')

from palette import COL, METH
AB = {'DRUG-seq': 'DS', 'BRB-seq': 'BRB', 'prime-seq': 'PS', 'BOBseq': 'BOB'}   # short labels for the legend-side gene counts (BRB-seq is not part of the preprint)
S = [l.rstrip('\n').split('\t') for l in open(SAMPLES) if l.strip()]
rng = np.random.default_rng(7); LOG = open(f'{OUT}/analysis.log', 'a'); MD = []
def log(*a): print(*a); print(*a, file=LOG); LOG.flush()
def md(s=''): MD.append(s)
import textwrap
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
def save(fig, name, caption=None, legend=True, handles=None, labels=None):
    """shared legend under the panels, wrapped caption under that; nothing drawn inside the axes. Space is computed from the caption length."""
    for a in fig.axes:
        if a.get_legend(): a.get_legend().remove()
    W, Hh = fig.get_size_inches(); cap_h = 0.0
    if caption:
        wrapped = textwrap.fill(caption, int(W * 13.5)); n = wrapped.count('\n') + 1; cap_h = (n * 0.135 + 0.12) / Hh
        fig.text(0.01, 0.006, wrapped, fontsize=7.6, va='bottom', ha='left', color='#333')
    if legend or handles:
        fig.legend(handles=handles or [Line2D([], [], color=COL[m], lw=2.5) for m in METHS], labels=labels or METHS, loc='lower center', ncol=len(labels or METHS), frameon=False, fontsize=9, bbox_to_anchor=(0.5, cap_h + 0.01))
    fig.tight_layout(rect=(0, cap_h + (0.07 if (legend or handles) else 0.01), 1, 0.96))
    _overlap_check(fig, name)
    fig.savefig(f'{OUT}/{name}.png', dpi=140); fig.savefig(f'{OUT}/{name}.svg'); plt.close(fig); log('wrote', name)
# ---------- load ----------
D = {}
for lab, bam, m in S:
    p = f'{NPZ}/{lab}.npz'
    if not os.path.exists(p): log('missing', p); continue
    z = np.load(p); D[lab] = {k: z[k] for k in ('gene', 't', 'L', 'mate', 'sense')}; D[lab]['method'] = m; D[lab]['n_unique'] = int(z['n_unique'])
    if 'genes' not in globals(): genes = z['genes']; tx_len = z['tx_len']
meth_of = {lab: D[lab]['method'] for lab in D}
POOL = {}
for m in METH:
    labs = [l for l in D if meth_of[l] == m]
    if labs: POOL[m] = {k: np.concatenate([D[l][k] for l in labs]) for k in ('gene', 't', 'L', 'mate', 'sense')}; POOL[m]['samples'] = labs
METHS = [m for m in METH if m in POOL]
log('samples:', {m: POOL[m]['samples'] for m in METHS})
def by_gene(d):
    g = d['gene']; o = np.argsort(g, kind='stable'); gs = g[o]; st = np.flatnonzero(np.r_[True, gs[1:] != gs[:-1]]); en = np.r_[st[1:], len(gs)]
    return {int(gs[s]): o[s:e] for s, e in zip(st, en)}
IDX = {m: by_gene(POOL[m]) for m in METHS}
def frac(d, ii): return d['t'][ii] / np.maximum(d['L'][ii], 1)
def gini(c):
    cs = np.sort(c); n = len(cs); return (2 * np.sum(np.arange(1, n + 1) * cs) / (n * cs.sum())) - (n + 1) / n if cs.sum() else 0
def profile_equal(d, gene_list, nb=20, sub=None):
    """mean scaled read-start density over genes, equal gene weight; optional subsample per gene"""
    P = np.zeros(nb); n = 0
    for g in gene_list:
        ii = IDX_cur[g] if sub is None else rng.choice(IDX_cur[g], sub, replace=False)
        b = np.minimum((frac(d, ii) * nb).astype(int), nb - 1); P += np.bincount(b, minlength=nb) / len(ii); n += 1
    return P / max(n, 1)
SUMMARY = {m: {} for m in METHS}
# ---------- A: breadth vs per-gene depth ----------
def section_A():
    global IDX_cur
    depths = [25, 50, 100, 200, 400, 800]; res = {m: {} for m in METHS}
    for d in depths:
        common = [g for g in IDX[METHS[0]] if all(len(IDX[m].get(g, [])) >= d for m in METHS)]
        for m in METHS:
            occ, hit5, hit3, gi = [], [], [], []
            for g in common:
                ii = rng.choice(IDX[m][g], d, replace=False); b = np.minimum((frac(POOL[m], ii) * 20).astype(int), 19); c = np.bincount(b, minlength=20)
                occ.append((c > 0).sum()); hit5.append(c[:2].sum() > 0); hit3.append(c[18:].sum() > 0); gi.append(gini(c))
            res[m][d] = dict(genes=len(common), occ=float(np.median(occ)), occ_mean=float(np.mean(occ)), hit5=float(np.mean(hit5)), hit3=float(np.mean(hit3)), gini=float(np.median(gi)))
    fig, ax = plt.subplots(1, 3, figsize=(15, 5.6))
    for m in METHS:
        ax[0].plot(depths, [res[m][d]['occ_mean'] for d in depths], 'o-', color=COL[m]); ax[1].plot(depths, [100 * res[m][d]['hit5'] for d in depths], 'o-', color=COL[m]); ax[2].plot(depths, [res[m][d]['gini'] for d in depths], 'o-', color=COL[m])
    for a, t, yl in zip(ax, ['Breadth: transcript bins (of 20) with a read start, mean over genes', "5' end reached: % of genes with a read in the 5'-most 10%", 'Unevenness: median Gini of the 20 bin counts (0 = uniform)'], ['bins of 20', '% of genes', 'Gini']):
        a.set_xscale('log'); a.set_xticks(depths); a.set_xticklabels(depths); a.set_xlabel('unique reads per gene (every gene subsampled to this number)'); a.set_title(t, fontsize=9); a.set_ylabel(yl); a.grid(ls=':', lw=.5)
    ax[0].axhline(20, color='grey', lw=.5, ls=':'); ax[0].set_ylim(0, 20.5)
    fig.suptitle('A. Per-gene coverage breadth at matched per-gene depth', fontsize=11)
    save(fig, 'A_breadth_vs_depth', caption=('TESTED: "random coverage can look like full coverage". Every gene is subsampled to exactly the same number of unique reads (x axis), so depth and expression cannot differ between methods; the mature '
        'transcript (exons of the Ensembl canonical protein-coding transcript) is cut into 20 equal bins and we count how many bins hold at least one read 5\' end (left), whether the 5\'-most 10% holds one (middle), and how uneven the 20 bin '
        'counts are (right). Gene sets per depth: ' + ', '.join(f'{d} reads: {res[METHS[0]][d]["genes"]:,} genes' for d in depths) + ' (genes with at least that many reads in every method). A method whose average profile is flat only because different genes '
        'are covered in different random patches would score low here; per-gene full-length coverage scores high and saturates. Unique reads (MAPQ 255); ' + READNOTE + '.'))
    md('## A. Per-gene breadth at matched depth (subsample exactly d unique reads per gene, 20 bins along the mature transcript)'); md('| method | ' + ' | '.join(f'd={d} (n={res[METHS[0]][d]["genes"]:,})' for d in depths) + ' |'); md('|---|' + '---|' * len(depths))
    for m in METHS: md(f'| {m} | ' + ' | '.join(f"{res[m][d]['occ_mean']:.1f} bins, 5' hit {100*res[m][d]['hit5']:.0f}%, Gini {res[m][d]['gini']:.2f}" for d in depths) + ' |')
    for m in METHS: SUMMARY[m].update(breadth200=res[m][200]['occ_mean'], hit5_200=res[m][200]['hit5'], gini200=res[m][200]['gini'])
    json.dump(res, open(f'{OUT}/A_breadth_vs_depth.json', 'w'), indent=1)
# ---------- B: per-gene position distributions ----------
def section_B():
    d0 = 200; common = [g for g in IDX[METHS[0]] if all(len(IDX[m].get(g, [])) >= d0 for m in METHS)]
    fig, ax = plt.subplots(1, 3, figsize=(15, 5.6)); data_med, data_3p, data_5p = [], [], []
    md(f'## B. Per-gene distributions at 200 reads per gene ({len(common):,} genes): where does the median read sit, and what share of reads is in the last 300 nt / first 200 nt'); md('| method | median of per-gene median position (0 = 5\', 1 = 3\') | IQR | per-gene share of reads in the 3\'-most 300 nt, median | per-gene share in the 5\'-most 200 nt, median |'); md('|---|---|---|---|---|')
    for m in METHS:
        medp, s3, s5 = [], [], []
        for g in common:
            ii = rng.choice(IDX[m][g], d0, replace=False); t = POOL[m]['t'][ii]; L = POOL[m]['L'][ii][0]
            medp.append(np.median(t) / L); s3.append(np.mean(t >= L - 300)); s5.append(np.mean(t < 200))
        data_med.append(medp); data_3p.append(s3); data_5p.append(s5); q = np.percentile(medp, [25, 50, 75])
        md(f'| {m} | {q[1]:.2f} | {q[0]:.2f}-{q[2]:.2f} | {100*np.median(s3):.0f}% | {100*np.median(s5):.0f}% |'); SUMMARY[m].update(median_pos=float(q[1]), share3p300_gene=float(np.median(s3)), share5p200_gene=float(np.median(s5)))
    for a, dat, t in zip(ax, [data_med, data_3p, data_5p], ["position of the median read (0 = 5' end, 1 = 3' end)", "share of the gene's reads in the 3'-most 300 nt", "share of the gene's reads in the 5'-most 200 nt"]):
        v = a.violinplot(dat, showmedians=True); [b.set_facecolor(COL[m]) for b, m in zip(v['bodies'], METHS)]; a.set_xticks(range(1, len(METHS) + 1)); a.set_xticklabels(METHS); a.set_title(t, fontsize=9); a.set_ylabel('per gene'); a.grid(axis='y', ls=':', lw=.5)
    fig.suptitle(f'B. Per-gene distributions at 200 unique reads per gene ({len(common):,} genes)', fontsize=11)
    save(fig, 'B_per_gene_distributions', legend=False, caption=('TESTED: "random coverage can look like full coverage", per gene. Same genes and subsampling as Figure A at 200 reads per gene. For every gene: where its median read sits along the mature transcript '
        '(0 = 5\' end, 1 = 3\' end; left), what share of its reads start in the 3\'-most 300 nt (middle) and in the 5\'-most 200 nt (right). Violins = distribution over genes, bar = median. If a flat average came from random patches, the median '
        'positions would scatter over 0-1 and the 3\' share would be bimodal; end-to-end coverage puts the median near the middle with a narrow spread. Each violin is one method; the x label names it.'))
# ---------- C/D: unscaled distance from the 3' and 5' end ----------
def dist_curves(end):
    global IDX_cur
    classes = [('0.5-1 kb', 500, 1000), ('1-2 kb', 1000, 2000), ('2-4 kb', 2000, 4000), ('4-8 kb', 4000, 8000), ('> 8 kb', 8000, 10**9)]
    xs = np.arange(0, 3001, 25); key = 'C' if end == 3 else 'D'; cuts = [300, 500, 1000] if end == 3 else [100, 200, 500]
    fig, ax = plt.subplots(1, 6, figsize=(22, 5.4), sharey=True); tab = {m: {} for m in METHS}; ngen = {}
    for k, (cname, lo, hi) in enumerate([('all', 500, 10**9)] + classes):
        for m in METHS:
            d = POOL[m]; glist = [g for g, ii in IDX[m].items() if len(ii) >= 100 and lo <= tx_len[g] < hi]
            if not glist: continue
            cdf = np.zeros(len(xs)); nn = 0
            for g in glist:
                ii = IDX[m][g]; dist = (d['L'][ii] - 1 - d['t'][ii]) if end == 3 else d['t'][ii]
                cdf += np.searchsorted(np.sort(dist), xs, side='right') / len(ii); nn += 1
            cdf /= nn; ax[k].plot(xs, 100 * cdf, color=COL[m], lw=2); ngen.setdefault(cname, []).append(f'{AB[m]} {nn:,}')
            tab[m][cname] = {c: float(np.interp(c, xs, cdf)) for c in cuts}
        ax[k].set_title((f'transcripts {cname}' if k else 'all transcripts >= 500 nt') + '\ngenes: ' + ' | '.join(ngen.get(cname, [])), fontsize=8); ax[k].set_xlabel(f"distance from the {end}' end (nt)"); ax[k].grid(ls=':', lw=.5); ax[k].set_ylim(0, 102)
    ax[0].set_ylabel('cumulative % of read starts (equal gene weight)')
    fig.suptitle(f"{key}. Unscaled: cumulative share of read starts within a given distance of the {end}' end, by mature-transcript length", fontsize=11)
    where = "from the 3' end" if end == 3 else "only on short transcripts"
    save(fig, f'{key}_distance_from_{end}p_end', caption=(f'TESTED: "short v long", without any scaling. For each gene with >= 100 unique reads, the distance (in nucleotides along the mature transcript) of every read 5\' end from the {end}\' end is '
        f'turned into a cumulative curve; curves are averaged over genes with equal weight (so abundant genes do not dominate). Panels split the genes by mature-transcript length. A 3\'-end method shows a steep rise within the first few '
        f'hundred nt {where}; the same method looks "even" in a scaled plot on short transcripts (its window is a large fraction) and "noisy-even" on long ones (its window is a small '
        'fraction and the rest is sparse). Reading in nucleotides removes that ambiguity. Panel titles give the number of genes per method.'))
    md(f"## {key}. Share of read starts within a distance of the {end}' end (equal gene weight, genes >= 100 reads)"); md('| method | class | ' + ' | '.join(f'<= {c} nt' for c in cuts) + ' |'); md('|---|---|' + '---|' * len(cuts))
    for m in METHS:
        for cname in tab[m]: md(f'| {m} | {cname} | ' + ' | '.join(f'{100*tab[m][cname][c]:.0f}%' for c in cuts) + ' |')
        if 'all' in tab[m]: SUMMARY[m].update({f'share{end}p_{cuts[0]}': tab[m]['all'][cuts[0]]})
    json.dump(tab, open(f'{OUT}/{key}_distance_from_{end}p_end.json', 'w'), indent=1)
# ---------- E/F: stratified scaled profiles ----------
def section_EF():
    global IDX_cur
    NG = {}
    fig, ax = plt.subplots(2, 4, figsize=(18, 9.6)); md('## E. Scaled read-start profile by expression quartile (genes >= 50 reads; equal gene weight; thirds 5\'/mid/3\')'); md('| method | Q1 (lowest) | Q2 | Q3 | Q4 (highest) |'); md('|---|---|---|---|---|')
    x = np.arange(20) * 5 + 2.5
    for m in METHS:
        IDX_cur = IDX[m]; d = POOL[m]; gl = [g for g, ii in IDX_cur.items() if len(ii) >= 50]; cnt = np.array([len(IDX_cur[g]) for g in gl]); qs = np.percentile(cnt, [25, 50, 75]); row = []
        for qi, (lo, hi) in enumerate([(0, qs[0]), (qs[0], qs[1]), (qs[1], qs[2]), (qs[2], np.inf)]):
            sel = [g for g, c in zip(gl, cnt) if lo <= c < hi]; P = profile_equal(d, sel); ax[0, qi].plot(x, 100 * P, color=COL[m], lw=2); th = [P[:7].sum(), P[7:13].sum(), P[13:].sum()]; row.append(f'{th[0]:.2f}/{th[1]:.2f}/{th[2]:.2f}'); NG.setdefault(('E', qi), []).append(f'{AB[m]} {len(sel):,}')
        md(f'| {m} | ' + ' | '.join(row) + ' |')
    for qi, t in enumerate(['Q1 lowest expressed', 'Q2', 'Q3', 'Q4 highest expressed']): ax[0, qi].set_title(f'E. {t} quartile of genes\ngenes: ' + ', '.join(NG.get(('E', qi), [])), fontsize=8); ax[0, qi].set_xlabel("position along the mature transcript, 5' to 3' (%)"); ax[0, qi].set_ylabel('% of read starts per 5% bin'); ax[0, qi].grid(ls=':', lw=.5); ax[0, qi].axhline(5, color='grey', lw=.5, ls=':')
    classes = [('0.5-1 kb', 500, 1000), ('1-2 kb', 1000, 2000), ('2-5 kb', 2000, 5000), ('> 5 kb', 5000, 10**9)]
    md('## F. Scaled read-start profile by mature-transcript length (genes >= 50 reads; thirds 5\'/mid/3\')'); md('| method | ' + ' | '.join(c for c, _, _ in classes) + ' |'); md('|---|' + '---|' * len(classes))
    for m in METHS:
        IDX_cur = IDX[m]; d = POOL[m]; row = []
        for ci, (cname, lo, hi) in enumerate(classes):
            sel = [g for g, ii in IDX_cur.items() if len(ii) >= 50 and lo <= tx_len[g] < hi]; P = profile_equal(d, sel); ax[1, ci].plot(x, 100 * P, color=COL[m], lw=2); th = [P[:7].sum(), P[7:13].sum(), P[13:].sum()]; row.append(f'{th[0]:.2f}/{th[1]:.2f}/{th[2]:.2f}'); NG.setdefault(('F', ci), []).append(f'{AB[m]} {len(sel):,}')
        md(f'| {m} | ' + ' | '.join(row) + ' |')
    for ci, (cname, _, _) in enumerate(classes): ax[1, ci].set_title(f'F. transcripts {cname}\ngenes: ' + ', '.join(NG.get(('F', ci), [])), fontsize=8); ax[1, ci].set_xlabel("position along the mature transcript, 5' to 3' (%)"); ax[1, ci].set_ylabel('% of read starts per 5% bin'); ax[1, ci].grid(ls=':', lw=.5); ax[1, ci].axhline(5, color='grey', lw=.5, ls=':')
    fig.suptitle('E/F. Scaled read-start profiles stratified by expression (top row) and by mature-transcript length (bottom row)', fontsize=11)
    save(fig, 'EF_stratified_profiles', caption=('TESTED: "removing low expressed genes" (top row) and "short v long" in the scaled view (bottom row). Each gene with >= 50 unique reads contributes its own profile (read 5\' ends in 20 bins '
        'along the mature transcript, normalized to the gene\'s read count), profiles are averaged with equal gene weight; one transcript per gene (Ensembl canonical, MANE Select). Top: genes split into expression quartiles per method '
        '(Q1 = lowest read count). Bottom: genes split by mature-transcript length. The dotted line at 5% is a perfectly uniform profile. If low-expressed genes flattened the average, Q4 would look much more 3\'-biased than Q1.'))
# ---------- G: BOBseq mates ----------
def section_G():
    global IDX_cur
    if 'BOBseq' not in POOL: return
    d = POOL['BOBseq']
    if not (d['mate'] > 0).any(): md('## G. BOBseq mates: skipped, this set has no mate information (single-end 50-nt reads)'); md(''); return
    x = np.arange(20) * 5 + 2.5; fig, ax = plt.subplots(1, 3, figsize=(16, 5.6)); HL = []
    md("## G. BOBseq: the two mates separately (mate 1 = R2, 5' end = template-switch point where the RT stopped; mate 2 = R1, 5' end = random-priming site)")
    md('| read end | % within 50 nt of the transcript 5\' end (cap reached) | % within 200 nt of 5\' end | % within 300 nt of the 3\' end | thirds 5\'/mid/3\' |'); md('|---|---|---|---|---|')
    for mate, lab, c in ((1, "R2 5' end = RT stop / template switch", '#0b6e63'), (2, "R1 5' end = random priming site", '#7fb8b0')):
        sel = d['mate'] == mate; sub = {k: d[k][sel] for k in ('gene', 't', 'L')}; IDX_cur = by_gene(sub); gl = [g for g, ii in IDX_cur.items() if len(ii) >= 50]
        P = profile_equal(sub, gl); ax[0].plot(x, 100 * P, color=c, lw=2); HL.append((Line2D([], [], color=c, lw=2.5), lab))
        for k, end in ((1, 5), (2, 3)):
            xs = np.arange(0, 2001, 25); cdf = np.zeros(len(xs))
            for g in gl:
                ii = IDX_cur[g]; dist = sub['t'][ii] if end == 5 else (sub['L'][ii] - 1 - sub['t'][ii]); cdf += np.searchsorted(np.sort(dist), xs, side='right') / len(ii)
            cdf /= len(gl); ax[k].plot(xs, 100 * cdf, color=c, lw=2)
            if end == 5: c50, c200 = np.interp(50, xs, cdf), np.interp(200, xs, cdf)
            else: c300 = np.interp(300, xs, cdf)
        md(f"| {lab} | {100*c50:.1f}% | {100*c200:.1f}% | {100*c300:.1f}% | {P[:7].sum():.2f}/{P[7:13].sum():.2f}/{P[13:].sum():.2f} |")
    ax[0].set_title('scaled profile of read 5\' ends, genes >= 50 reads', fontsize=9); ax[0].set_xlabel("position along the mature transcript (%)"); ax[1].set_title("cumulative distance from the transcript 5' end (cap)", fontsize=9); ax[1].set_xlabel('nt'); ax[2].set_title("cumulative distance from the 3' end", fontsize=9); ax[2].set_xlabel('nt')
    for a in ax: a.grid(ls=':', lw=.5)
    ax[0].axhline(5, color='grey', lw=.5, ls=':'); ax[0].set_ylabel('% of read 5\' ends per 5% bin'); ax[1].set_ylabel('cumulative % of read 5\' ends'); ax[2].set_ylabel('cumulative % of read 5\' ends')
    fig.suptitle('G. BOBseq mechanism: template-switch points (R2 5\' ends) versus random-priming sites (R1 5\' ends), 10 samples pooled', fontsize=11)
    save(fig, 'G_bobseq_mates', legend=False, handles=[h for h, _ in HL], labels=[l for _, l in HL], caption=('TESTED: is the BOBseq "double peak" (5\' shoulder plus 3\' rise) an artefact or the chemistry? In a BOBseq pair, R2 starts at the template-switch point (where the reverse transcriptase '
        'stopped: at the cap if it went all the way) and R1 starts at the random-priming site. The two 5\' ends are plotted separately: scaled profile over genes with >= 50 reads (left, equal gene weight; dotted = uniform), and '
        'the cumulative distance of each end from the transcript 5\' end (middle) and from the 3\' end (right). A pile-up of R2 ends within the first bins means the RT reached the cap; a pile-up of R1 ends in the last bins means '
        'priming is enriched near the poly(A) end. Both mates are used everywhere else in this report; this figure separates them.'))
# ---------- H: replicate fidelity ----------
def section_H():
    pairs = []
    for m in METHS:
        labs = POOL[m]['samples']
        for i in range(len(labs)):
            for j in range(i + 1, len(labs)):
                if m != 'BOBseq' or (labs[i].split('-rep')[0] == labs[j].split('-rep')[0]): pairs.append((m, labs[i], labs[j]))
    if 'bobseq-hek' in D and 'primeseq-hek' in D: pairs.append(('cross-method control', 'bobseq-hek', 'primeseq-hek'))
    md('## H. Replicate fidelity of the per-gene profile (genes with >= 100 unique reads in both samples, 20 bins, Pearson r per gene)'); md('| pair | genes | median r | % genes r > 0.8 |'); md('|---|---|---|---|')
    res = collections.defaultdict(list); labels = []
    for m, a, b in pairs:
        ia, ib = by_gene(D[a]), by_gene(D[b]); common = [g for g in ia if g in ib and len(ia[g]) >= 100 and len(ib[g]) >= 100]; rs = []
        for g in common:
            pa = np.bincount(np.minimum((frac(D[a], ia[g]) * 20).astype(int), 19), minlength=20); pb = np.bincount(np.minimum((frac(D[b], ib[g]) * 20).astype(int), 19), minlength=20)
            if pa.std() > 0 and pb.std() > 0: rs.append(np.corrcoef(pa, pb)[0, 1])
        rs = np.array(rs); res[m].append(rs); labels.append(f'{a} vs {b}'); md(f'| {a} vs {b} | {len(rs):,} | {np.median(rs):.2f} | {100*np.mean(rs > .8):.0f}% |')
        if m in SUMMARY: SUMMARY[m].setdefault('rep_r', []).append(float(np.median(rs)))
    fig, ax = plt.subplots(figsize=(max(9, 0.95 * len(labels)), 6.2)); allr = [r for m in res for r in res[m]]; cols = [COL.get(m, 'grey') for m in res for _ in res[m]]
    bp = ax.boxplot(allr, showfliers=False, patch_artist=True, medianprops=dict(color='k')); [b.set_facecolor(c) for b, c in zip(bp['boxes'], cols)]; ax.set_xticks(range(1, len(labels) + 1)); ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8); ax.set_ylabel('Pearson r of the 20-bin profile, per gene'); ax.grid(axis='y', ls=':', lw=.5)
    ax.set_title('H. Replicate fidelity of the per-gene coverage profile (genes with >= 100 unique reads in both samples of a pair)', fontsize=11)
    save(fig, 'H_replicate_fidelity', legend=False, caption=('TESTED: is the coverage pattern of a gene reproducible between replicates ("peaks and troughs identical across samples")? For every pair of replicate samples and every gene with '
        '>= 100 unique reads in both, the 20-bin read-start profile of sample 1 is correlated (Pearson) with that of sample 2; boxes show the distribution over genes (median line, quartiles, whiskers 1.5 IQR, outliers hidden). '
        'The grey box is a cross-method pair (BOBseq HEK vs prime-seq HEK) as the negative control. CAVEAT: a profile that is a single 3\' spike correlates trivially, so a high r for a 3\'-end method reflects the spike, not '
        'reproducibility of coverage; BOBseq\'s r is over 20 informative bins. Colours: ' + ', '.join(f'{m} {COL[m]}' for m in METH) + '.'))
# ---------- I: strand ----------
def section_I():
    md('## I. Strand of the read start relative to the gene (unique exonic reads; BOBseq: mate 1 = R2 and mate 2 = R1 separately)'); md('| sample | method | reads positioned | % sense (mate 1 / single) | % sense mate 2 |'); md('|---|---|---|---|---|')
    for lab in D:
        d = D[lab]; m1 = d['mate'] != 2; m2 = d['mate'] == 2
        md(f"| {lab} | {meth_of[lab]} | {len(d['t']):,} | {100*d['sense'][m1].mean():.1f}% | {('%.1f%%' % (100*d['sense'][m2].mean())) if m2.any() else '-'} |")
# ---------- K: summary ----------
def section_K():
    md('## K. One-table summary (candidates to replace the Picard histogram in the paper)'); md("| method | breadth at 200 reads/gene (bins of 20) | 5' end reached at 200 reads | median per-gene median position | read starts within 300 nt of 3' end | within 100 nt of 5' end | replicate profile r (median) |"); md('|---|---|---|---|---|---|---|')
    for m in METHS:
        s = SUMMARY[m]; rr = np.median(s['rep_r']) if s.get('rep_r') else float('nan')
        md(f"| {m} | {s.get('breadth200', float('nan')):.1f} | {100*s.get('hit5_200', float('nan')):.0f}% | {s.get('median_pos', float('nan')):.2f} | {100*s.get('share3p_300', float('nan')):.0f}% | {100*s.get('share5p_100', float('nan')):.0f}% | {rr:.2f} |")
for name, fn in (('A', section_A), ('B', section_B), ('C', lambda: dist_curves(3)), ('D', lambda: dist_curves(5)), ('EF', section_EF), ('G', section_G), ('H', section_H), ('I', section_I), ('K', section_K)):
    try: log('section', name); fn()
    except Exception: log('SECTION', name, 'FAILED'); log(traceback.format_exc()); md(f'\n**Section {name} failed, see analysis.log**\n')
open(f'{OUT}/COVERAGE_ARCHITECTURE.md', 'w').write('# Coverage architecture along the mature transcript (' + ', '.join(f'{m}: {len(POOL[m]["samples"])} samples' for m in METHS) + ')\n\nRead 5\' ends of unique (MAPQ 255) alignments placed on the Ensembl canonical protein-coding transcript of each gene (exons only, Ensembl 113). ' + ('BOBseq = both mates unless stated.' if SET == 'native' else READNOTE + '.') + ' "Equal gene weight" = every gene contributes its own normalized profile, like Picard, but with one transcript per gene and a minimum read count.\n\n' + '\n'.join(MD) + '\n')
log('ANALYSIS DONE')

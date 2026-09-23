#!/usr/bin/env python3
"""Stage figures: browser-style coverage tracks with junction arcs for every locus of config/loci.tsv.

The benchmark's gene_coverage_tracks.py (one filled coverage track per BAM, the peak depth as its only y tick, the canonical
transcript below, 5' to 3' left to right) extended with:
  * four views per locus: gene_with_introns (the canonical transcript span, introns to scale), gene_exons_only (exons
    concatenated, ticks = hg38 position of the exon starts), zoom and zoom_arcs (every exon that carries a site of the LSV's
    junctions, or of the positive control's junctions, plus the event exon, widened to three exons with the nearest canonical
    neighbours, terminal exons cut to EXON_CAP nt, +- ZOOM_PAD; without and with arcs); zoom_site_arcs where config/loci.tsv
    gives a site_window;
  * junction arcs on zoom_arcs: one per junction of the catalogue (ATLAS_JUNCTIONS, shipped for these loci) with >= MIN_ARC
    fragments and both ends in view and >= ARC_MIN_FRAC of the strongest catalogue junction in that track's view, no exception
    for the LSV's own junctions; apex height by intron length, line width by fragment count (log scale), the count at the apex;
  * track sets on the depth-matched BAMs of prepare_bams.py (config/loci.tsv column `sets`), tracks coloured by condition;
  * the event exon in red on the transcript row when it is not part of the canonical transcript.
Coverage = depth per base from the aligned blocks of every alignment record (spliced gaps are not covered); a junction = every N
in a CIGAR with >= MIN_ANCHOR nt aligned on both sides on that read, counted once per read name. Output:
figures/<locus>_<tag>/<locus>_<view>_<set>.svg (text stays text), PNG previews under .cache/png/,
results/values/<locus>_{depth,junctions}_<set>.tsv, results/gene_models.tsv, results/windows.tsv.
    gene_coverage_tracks.py [LOCUS ...] [--sets replicates pseudobulk ...] [--views gene_with_introns gene_exons_only zoom zoom_arcs zoom_site_arcs]"""
import argparse, csv, os, re, sys
import numpy as np, pysam, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coverage_lib as cl
from settings import CONFIG, RESULTS, FIGURES, CACHE, BAMS, TRACK_SETS, MIN_ANCHOR, MIN_ARC, ARC_MIN_FRAC, ATLAS_JUNCTIONS, EXON_CAP, ZOOM_PAD
from palette import COL, ARC, GENECOL, EVENTCOL

# --- fonts: Arial named in the SVG, Liberation Sans (metric-identical, installed) doing the layout
from matplotlib import font_manager as _fm
import dataclasses as _dc
for _p in (q for q in _fm.findSystemFonts() if 'LiberationSans-' in q):
    _fm.fontManager.ttflist.append(_dc.replace(_fm.ttfFontProperty(_fm.get_font(_p)), name='Arial'))
plt.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype': 42, 'font.family': 'Arial', 'font.size': 8, 'axes.linewidth': 0.7,
                     'text.color': '#231F20', 'axes.edgecolor': '#231F20', 'xtick.color': '#231F20', 'ytick.color': '#231F20'})
PAD = 0.03
VALUES = os.path.join(RESULTS, 'values'); PNG = os.path.join(CACHE, 'png')
os.makedirs(VALUES, exist_ok=True); os.makedirs(PNG, exist_ok=True)


# --- data -----------------------------------------------------------------------------------------------------------
def depth(bam, chrom, s0, e0):
    d = np.zeros(e0 - s0 + 1, np.int64)
    for a in pysam.AlignmentFile(bam, 'rb').fetch(chrom, s0, e0):
        for bs, be in a.get_blocks():
            d[min(max(bs - s0, 0), e0 - s0)] += 1
            d[min(max(be - s0, 0), e0 - s0)] -= 1
    return np.cumsum(d)[: e0 - s0]


def junctions(bam, chrom, s0, e0):
    """{(intron_start, intron_end): fragments} for every splice junction of a read overlapping [s0, e0), with >= MIN_ANCHOR nt
    aligned on both sides of the splice on that read; a fragment (read name) counts once per junction."""
    seen = {}
    for a in pysam.AlignmentFile(bam, 'rb').fetch(chrom, s0, e0):
        if a.cigartuples is None or not any(op == 3 for op, _ in a.cigartuples):
            continue
        pos, runs = a.reference_start, []                      # runs of (kind, ref_start, ref_end): 'M' aligned, 'N' intron
        for op, n in a.cigartuples:
            if op in (0, 7, 8):
                if runs and runs[-1][0] == 'M' and runs[-1][2] == pos:
                    runs[-1] = ('M', runs[-1][1], pos + n)
                else:
                    runs.append(('M', pos, pos + n))
                pos += n
            elif op == 3:
                runs.append(('N', pos, pos + n)); pos += n
            elif op == 2:
                pos += n                                        # deletion: reference advances, no aligned bases
        for i, (k, s, e) in enumerate(runs):
            if k != 'N' or i == 0 or i == len(runs) - 1:
                continue
            left, right = runs[i - 1], runs[i + 1]
            if left[0] != 'M' or right[0] != 'M' or left[2] - left[1] < MIN_ANCHOR or right[2] - right[1] < MIN_ANCHOR:
                continue
            if e <= s0 or s >= e0:
                continue
            seen.setdefault((s, e), set()).add(a.query_name)
    return {k: len(v) for k, v in seen.items()}


# --- windows ---------------------------------------------------------------------------------------------------------
def exon_at(exons_all, site, side):
    """The annotated exon whose boundary is `site` (side 'end' = exon ends at site, 'start' = exon starts at site); the shortest
    one, cut to EXON_CAP nt away from the site. None if no exon has that boundary."""
    cands = [(s, e) for (s, e) in exons_all if (e == site if side == 'end' else s == site)]
    if not cands:
        return None
    s, e = min(cands, key=lambda x: x[1] - x[0])
    if e - s > EXON_CAP:
        s, e = (e - EXON_CAP, e) if side == 'end' else (s, s + EXON_CAP)
    return (s, e)


def zoom_window(loc, tx, exons_all):
    """Exons touched by the LSV sites (+ the event exon), widened to >= 3 exons with the nearest canonical neighbours."""
    L, tid, chrom, strand, starts, ends = tx
    canon = list(zip(starts, ends))
    touched = set()
    for j in loc['members'].split(';'):
        c, span, st = j.split(':'); s, e = map(int, span.split('-'))     # intron [s, e): donor exon ends at s, acceptor exon starts at e
        for site, side in ((s, 'end'), (e, 'start')):
            ex = exon_at(exons_all, site, side)
            touched.add(ex if ex else ((site - 50, site) if side == 'end' else (site, site + 50)))
    if loc['event_exon']:
        c, span = loc['event_exon'].split(':'); s, e = map(int, span.split('-')); touched.add((s, e))
    touched = sorted(touched)
    lo, hi = touched[0][0], touched[-1][1]
    n_exons = len({(s, e) for (s, e) in canon if s < hi and e > lo} | set(touched))
    left = [x for x in canon if x[1] <= lo]; right = [x for x in canon if x[0] >= hi]
    while n_exons < 3 and (left or right):
        if left and (not right or (lo - left[-1][1]) <= (right[0][0] - hi)):
            x = left.pop(); lo = min(lo, max(x[0], x[1] - EXON_CAP))
        else:
            x = right.pop(0); hi = max(hi, min(x[1], x[0] + EXON_CAP))
        n_exons += 1
    pad = int(ZOOM_PAD * (hi - lo))
    return lo - pad, hi + pad, touched


# --- drawing ---------------------------------------------------------------------------------------------------------
def arc_patch(x1, y1, x2, y2, h, **kw):
    verts = [(x1, y1), ((x1 + x2) / 2, h), (x2, y2)]
    return PathPatch(Path(verts, [Path.MOVETO, Path.CURVE3, Path.CURVE3]), fill=False, **kw)


def binned(x, y, n=2000):
    """The coverage for drawing: block maxima over <= n bins (the per-base values are in results/values/). Keeps the SVG
    at ~2,000 path points per track instead of one per base -- Illustrator-friendly, and nothing is lost at 6.8 in."""
    x, y = np.asarray(x), np.asarray(y)
    if len(x) <= n:
        return x, y
    k = int(np.ceil(len(x) / n)); m = len(x) // k * k
    xb = x[:m].reshape(-1, k).mean(1); yb = y[:m].reshape(-1, k).max(1)
    if m < len(x):
        xb, yb = np.r_[xb, x[m:].mean()], np.r_[yb, y[m:].max()]
    return xb, yb


def draw(loc, x, tracks, arcs, exon_x, event_x, ticks, ticklabels, xlim, name, note, header, joined=False):
    """tracks: [(track id, label, condition, depth array)]; arcs: None (no arcs: the original's layout, 0.95 in per track,
    y up to 1.3 x peak) or {track id: [(x1, x2, count)]} in the x units of the view (1.15 in per track, arcs above the coverage)."""
    n = len(tracks); with_arcs = arcs is not None
    fig, axs = plt.subplots(n + 1, 1, figsize=(6.8, (1.15 if with_arcs else 0.95) * n + 0.75), gridspec_kw={'height_ratios': [1] * n + [0.62], 'hspace': 0.12}, sharex=True)
    for ax, (tid, label, cond, y) in zip(axs, tracks):
        pk = int(y.max())
        xb, yb = binned(x, y)
        ax.fill_between(xb, 0, yb, step='mid', color=COL[cond], lw=0, alpha=0.7 if with_arcs else 1.0)
        ax.set_ylim(0, max(pk, 1) * (1.8 if with_arcs else 1.3))
        ax.set_yticks([pk]); ax.set_yticklabels([f'{pk:,}'], fontsize=9)
        for sp in ('top', 'right', 'bottom'):
            ax.spines[sp].set_visible(False)
        ax.axhline(0, color=COL[cond], lw=0.4, alpha=0.5)
        ax.tick_params(axis='x', bottom=False, labelbottom=False); ax.tick_params(axis='y', length=4)
        ax.text(0.015, 0.97, label, color=COL[cond], fontsize=9.5, transform=ax.transAxes, va='top', ha='left', zorder=8,
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.95, pad=1.2))
        # junction arcs: apex by intron length rank, width by count, count at the apex
        lo_x, hi_x = min(xlim), max(xlim)
        A = sorted((t for t in (arcs or {}).get(tid, []) if lo_x <= t[0] <= hi_x and lo_x <= t[1] <= hi_x), key=lambda t: abs(t[1] - t[0]))   # both ends in view
        if A:
            nmax = max(c for _, _, c in A); labelled = set(range(len(A)))
            xs = np.asarray(x); labels = []
            for i, (x1, x2, c) in enumerate(A):
                y1 = float(np.interp(x1, xs, y)); y2 = float(np.interp(x2, xs, y))
                h = max(pk, 1) * (1.05 + 0.55 * (i / max(len(A) - 1, 1)))
                lw = 0.4 + 1.3 * (np.log10(c) / np.log10(nmax) if nmax > 1 else 1)
                y1, y2 = min(y1, pk), min(y2, pk)
                ax.add_patch(arc_patch(x1, y1, x2, y2, 2 * h - (y1 + y2) / 2, edgecolor=ARC[cond], lw=lw, alpha=0.9, zorder=4))   # control point so the apex sits at h
                if i in labelled:
                    labels.append(((x1 + x2) / 2, h * 1.02, c))
            # count labels: apex of the arc; labels whose apexes are within 2 % of the view of each other are stacked
            placed = []
            for xm, yb, c in sorted(labels, key=lambda t: (t[1], t[0])):        # lowest arcs first: a stack's order follows the arcs'
                near = [py for px, py in placed if abs(px - xm) < 0.02 * (hi_x - lo_x)]
                yl = max([yb] + [py + 0.2 * max(pk, 1) for py in near]); placed.append((xm, yl))
                ax.text(xm, yl, f'{c:,}', fontsize=6, color=ARC[cond], ha='center', va='bottom', zorder=7,
                        bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=0.6))
            if placed:
                ax.set_ylim(0, max(max(pk, 1) * 1.8, max(py for _, py in placed) + 0.25 * max(pk, 1)))   # room for stacked labels
    ax = axs[-1]
    ax.set_ylim(-1.35, 1)
    ax.plot([exon_x[0][0], exon_x[-1][1]], [0, 0], color=GENECOL, lw=0.9, zorder=1)
    for a, b in exon_x:
        ax.add_patch(plt.Rectangle((min(a, b), -0.42), abs(b - a), 0.84, color=GENECOL, lw=0, zorder=2))
    for a, b, lab in event_x:
        ax.add_patch(plt.Rectangle((min(a, b), -0.42), abs(b - a), 0.84, facecolor=EVENTCOL, edgecolor='none', zorder=3))
        ax.text((a + b) / 2, 0.55, lab, fontsize=6.5, color=EVENTCOL, ha='center', va='bottom', zorder=4)
    if joined:
        [ax.plot([a, a], [-0.42, 0.42], color='white', lw=1.1, zorder=3.5, solid_capstyle='butt') for a, _ in exon_x[1:]]
    for sp in ('top', 'right', 'left'):
        ax.spines[sp].set_visible(False)
    ax.spines['bottom'].set_position(('data', -1.0))
    ax.set_yticks([]); ax.set_xticks(ticks); ax.set_xticklabels(ticklabels, fontsize=7); ax.tick_params(axis='x', length=3)
    ax.text(-0.012, 0.69, loc['gene'], color=GENECOL, fontsize=9.5, transform=ax.transAxes, ha='right', va='center')
    ax.text(1.012, 0.69, '3p', color=GENECOL, fontsize=9.5, transform=ax.transAxes, ha='left', va='center')
    ax.text(-0.012, 0.13, 'hg38', fontsize=8.5, transform=ax.transAxes, ha='right', va='center')
    ax.text(1.0, -0.5, note, fontsize=6.5, color='#666', transform=ax.transAxes, ha='right', va='top')
    axs[0].set_xlim(*xlim)
    axs[0].text(0.0, 1.06, header, fontsize=7, color='#444', transform=axs[0].transAxes, ha='left', va='bottom')
    d = os.path.join(FIGURES, loc['dir']); os.makedirs(d, exist_ok=True)          # <gene>_<source tag>, config/loci.tsv
    fig.savefig(os.path.join(d, name + '.svg'), bbox_inches='tight', pad_inches=0.04)
    fig.savefig(os.path.join(PNG, name + '.png'), bbox_inches='tight', pad_inches=0.04, dpi=160)
    plt.close(fig)
    print('wrote', name, flush=True)


# --- per locus ---------------------------------------------------------------------------------------------------------
def header_of(loc):
    if loc['source'] == 'positive_controls':
        kind = 'risdiplam' if loc['contrast'] == 'ris_low' else 'CHX / NMD'
        tx = f" ({loc['event_transcript']})" if loc['event_transcript'] else ''
        return f"{loc['gene']} {loc['event_label']}{tx} — {kind} positive control"
    ns = '' if int(loc['hit']) else ', not significant (padded top-10)'
    return (f"{loc['gene']} — screen_splicing {loc['contrast'].replace(';', ' + ')} rank {loc['rank']}, ΔPSI {float(loc['dpsi']):+.2f}, "
            f"padj {float(loc['padj']):.2g}{ns}; LSV {loc['lsv'].split(':', 1)[1]} ({loc['members'].count(';') + 1} junctions), top junction {loc['junction']}")


def run_locus(loc, tx, exons_all, sets, views, tracks_by_set):
    L, tid, chrom, strand, starts, ends = tx
    rev = strand == '-'
    s0, e0 = starts[0], ends[-1]
    ev = None
    if loc['event_exon']:
        c, span = loc['event_exon'].split(':'); a, b = map(int, span.split('-'))
        if (a, b) not in set(zip(starts, ends)):
            ev = (a, b, loc['event_label']); s0, e0 = min(s0, a), max(e0, b)
    zlo, zhi, touched = zoom_window(loc, tx, exons_all.get(loc['gene'], {}))
    zlo, zhi = max(zlo, 0), zhi
    lo, hi = min(s0, zlo), max(e0, zhi)
    chrom_bam = chrom
    window_row = f"{loc['locus']}\t{chrom}\t{strand}\t{tid}\t{s0 + 1}-{e0}\t{zlo + 1}-{zhi}\t" + ';'.join(f'{a + 1}-{b}' for a, b in touched) + '\n'
    for setname in sets:
        tr = tracks_by_set[setname]
        D = {t: depth(bam, chrom_bam, lo, hi) for t, _, _, bam in tr}
        J = {t: junctions(bam, chrom_bam, lo, hi) for t, _, _, bam in tr}
        alljct = sorted({k for t in J for k in J[t]}, key=lambda k: (k[0], k[1]))
        with open(os.path.join(VALUES, f"{loc['locus']}_junctions_{setname}.tsv"), 'w') as f:
            f.write('chrom\tintron_start_1based\tintron_end\tlength\t' + '\t'.join(t for t, _, _, _ in tr) + '\n')
            for s, e in alljct:
                f.write(f'{chrom}\t{s + 1}\t{e}\t{e - s}\t' + '\t'.join(str(J[t].get((s, e), 0)) for t, _, _, _ in tr) + '\n')
        with open(os.path.join(VALUES, f"{loc['locus']}_depth_{setname}.tsv"), 'w') as f:
            f.write('chrom\tpos_hg38\tin_canonical_exon\t' + '\t'.join(t for t, _, _, _ in tr) + '\n')
            inex = np.zeros(hi - lo, bool)
            for a, b in zip(starts, ends):
                inex[max(a - lo, 0): max(b - lo, 0)] = True
            for i in range(hi - lo):
                f.write(f'{chrom}\t{lo + i + 1}\t{int(inex[i])}\t' + '\t'.join(str(int(D[t][i])) for t, _, _, _ in tr) + '\n')
        depthnote = f"MAPQ-255 fragments, {setname} set depth-matched (results/track_depths.tsv)"
        # ---- genomic views (whole gene; zoom without / with arcs; the site window): introns to scale, 5' left ----
        for view in [v for v in views if v != 'gene_exons_only']:
            if view == 'zoom_site_arcs':
                if not loc.get('site_window'):
                    continue
                w0, w1 = map(int, loc['site_window'].split(':')[1].split('-')); a0, b0 = w0 - 1, w1
            else:
                a0, b0 = (s0, e0) if view == 'gene_with_introns' else (zlo, zhi)
            gx = np.arange(a0, b0) + 1
            span = b0 - a0
            step = [v for v in (20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000) if span / v <= 6][0]
            tk = [t for t in range((a0 // step + 1) * step, b0, step)]
            xlim = (b0 + PAD * span, a0 - PAD * span) if rev else (a0 - PAD * span, b0 + PAD * span)
            tracks = [(t, lab, cond, D[t][a0 - lo: b0 - lo]) for t, lab, cond, _ in tr]
            arcs = None
            if view in ('zoom_arcs', 'zoom_site_arcs'):
                arcs = {}
                for t, _, _, _ in tr:
                    inview = {k: c for k, c in J[t].items() if c >= MIN_ARC and k[0] >= a0 and k[1] <= b0 and (chrom, k[0], k[1]) in CATALOGUE}
                    nmax = max(inview.values(), default=0)
                    arcs[t] = [(s + 0.5, e + 0.5, c) for (s, e), c in inview.items() if c >= ARC_MIN_FRAC * nmax]
            exon_x = [(a + 1, b) for a, b in zip(starts, ends) if b > a0 and a < b0]
            event_x = [(ev[0] + 1, ev[1], ev[2])] if ev and ev[1] > a0 and ev[0] < b0 else []
            note = (f"chr{chrom}, {tid}, {len(starts)} exons, {strand} strand; {depthnote}" if view == 'gene_with_introns'
                    else f"chr{chrom}:{a0 + 1:,}-{b0:,} ({span:,} nt), {tid}, {strand} strand; {depthnote}")
            if arcs is not None:
                note += f"\narcs: catalogue junctions with >= {MIN_ARC} fragments and >= {ARC_MIN_FRAC:.0%} of the strongest catalogue junction in view"
            draw(loc, gx, tracks, arcs, exon_x, event_x, tk, [f'{t:,}' for t in tk], xlim, f"{loc['locus']}_{view}_{setname}", note, header_of(loc))
        # ---- exon-only view: exons (canonical + the event exon) concatenated 5' to 3' ----
        if 'gene_exons_only' in views:
            ex = sorted(set(zip(starts, ends)) | ({(ev[0], ev[1])} if ev else set()))
            ex = ex[::-1] if rev else ex
            off = np.r_[0, np.cumsum([b - a for a, b in ex])]; Lx = int(off[-1])
            T = {t: np.concatenate([(D[t][a - lo: b - lo][::-1] if rev else D[t][a - lo: b - lo]) for a, b in ex]) for t, _, _, _ in tr}
            exon_x = [(off[i], off[i + 1]) for i in range(len(ex)) if not (ev and ex[i] == (ev[0], ev[1]))]
            event_x = [(off[i], off[i + 1], ev[2]) for i in range(len(ex)) if ev and ex[i] == (ev[0], ev[1])]
            sel = [0]
            [sel.append(i) for i in range(1, len(ex)) if off[i] - off[sel[-1]] >= Lx / 7]
            tracks = [(t, lab, cond, T[t]) for t, lab, cond, _ in tr]
            draw(loc, np.arange(Lx) + 0.5, tracks, None, exon_x, event_x, [off[i] for i in sel], [f'{(ex[i][1] if rev else ex[i][0] + 1):,}' for i in sel],
                 (-PAD * Lx, Lx * (1 + PAD)), f"{loc['locus']}_gene_exons_only_{setname}",
                 f"chr{chrom}, {tid}, introns removed ({len(ex)} exons, {Lx:,} nt, to scale); ticks = hg38 position of the exon start; {depthnote}",
                 header_of(loc), joined=True)
    return tid, chrom, strand, starts, ends, window_row


def load_catalogue(chroms):
    """{(chrom, intron_start, intron_end)} of the junction catalogue on the loci's chromosomes (ATLAS_JUNCTIONS: chrom,
    0-based intron start, intron end = last intron base 1-based, as the junction tables; header line skipped)."""
    cat = set()
    with open(ATLAS_JUNCTIONS) as fh:
        for line in fh:
            c, s, e = line.split('\t', 3)[:3]
            if c in chroms and s.isdigit():
                cat.add((c, int(s), int(e)))
    return cat


CATALOGUE = set()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('loci', nargs='*')
    ap.add_argument('--sets', nargs='+', default=list(TRACK_SETS))
    ap.add_argument('--views', nargs='+', default=['gene_with_introns', 'gene_exons_only', 'zoom', 'zoom_arcs', 'zoom_site_arcs'])
    a = ap.parse_args()
    loci = list(csv.DictReader(open(os.path.join(CONFIG, 'loci.tsv')), delimiter='\t'))
    if a.loci:
        loci = [l for l in loci if l['locus'] in a.loci]
    tracks_by_set = {}
    for r in csv.DictReader(open(os.path.join(RESULTS, 'track_depths.tsv')), delimiter='\t'):
        tracks_by_set.setdefault(r['set'], []).append((r['track'], r['label'], r['condition'], os.path.join(os.path.dirname(RESULTS), r['bam'])))
    tx, exons_all = cl.build([l['gene'] for l in loci])
    CATALOGUE.update(load_catalogue({tx[l['gene']][2] for l in loci}))
    print(f'junction catalogue: {len(CATALOGUE):,} junctions on the loci chromosomes', flush=True)
    gm, wf = os.path.join(RESULTS, 'gene_models.tsv'), os.path.join(RESULTS, 'windows.tsv')
    models = {l.split('\t')[0]: l for l in open(gm)} if os.path.exists(gm) else {}
    windows = {l.split('\t')[0]: l for l in open(wf)} if os.path.exists(wf) else {}
    for loc in loci:
        sets = [st for st in loc['sets'].split(';') if st in a.sets and st in tracks_by_set]      # the locus's own track sets
        tid, chrom, strand, starts, ends, wrow = run_locus(loc, tx[loc['gene']], exons_all, sets, a.views, tracks_by_set)
        models[loc['locus']] = f"{loc['locus']}\t{chrom}\t{strand}\t{tid}\t{','.join(str(s + 1) for s in starts)}\t{','.join(str(e) for e in ends)}\n"
        windows[loc['locus']] = wrow
    order = [l['locus'] for l in csv.DictReader(open(os.path.join(CONFIG, 'loci.tsv')), delimiter='\t')]
    with open(gm, 'w') as f:
        f.write('locus\tchrom\tstrand\ttranscript\texon_starts_1based\texon_ends\n')
        f.writelines(models[k] for k in order if k in models)
    with open(wf, 'w') as f:
        f.write('locus\tchrom\tstrand\ttranscript\tgene_span_1based\tzoom_window_1based\tlsv_exons_1based\n')
        f.writelines(windows[k] for k in order if k in windows)


if __name__ == '__main__':
    main()

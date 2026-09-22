#!/usr/bin/env python3
"""How far reverse transcription got: template-switch site to transcript 3' end.

For a poly-dT primed BOBseq library reverse transcription starts at the poly(A) tail and runs
towards the 5' end of the transcript until it stops or the template switch happens. The
distance from the template-switch site to the transcript 3' end is therefore the length of
the RT product (the cDNA insert). If labelled sites stop RT the distribution shifts left with
labelling density; if RT reads through it reflects processivity and transcript length.

Per read (primary alignment, uniquely mapped = NH:i:1, inside a named protein-coding gene,
not rRNA/rDNA and not mitochondrial):
  anchor   = the aligned end nearest the transcript 5' end, i.e. the alignment start on a
             '+' gene and the alignment end on a '-' gene (the template-switch site; a
             soft-clipped read end cannot move it)
  distance = anchor to the nearest annotated transcript 3' end of that gene downstream of the
             anchor (strand-aware); reads with no downstream end in their gene are skipped

Writes results['rt_extent'] per barcode: n, per-25-nt histogram 0..3000 (last bin >= 3000),
per-10-nt histogram 0..400, quartiles, % under 150 nt, % under 500 nt, % over 2 kb.
The report figure (section 4) reads only this dict, so PDF re-renders do not rescan the BAM.
"""
import os, re, sys, bisect, subprocess
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import star_taxonomy as st
except ImportError:                       # Illumina-suffixed module name
    import star_taxonomy_illumina as st

BINW, MAXD, FBINW, FMAX = 25, 3000, 10, 400
_NAME_RE = re.compile(r'gene_name "([^"]+)"')


def _tx_ends_impl(path):
    """{(contig, gene_name): (strand, sorted transcript 3' ends)} from the GTF."""
    ends = defaultdict(set); strand = {}
    with open(path) as f:
        for ln in f:
            if ln[0] == '#' or '\ttranscript\t' not in ln:
                continue
            p = ln.split('\t')
            m = _NAME_RE.search(p[8])
            if not m:
                continue
            k = (p[0], m.group(1)); strand[k] = p[6]
            ends[k].add(int(p[4]) if p[6] == '+' else int(p[3]))
    return {k: (strand[k], sorted(v)) for k, v in ends.items()}


def load_tx_ends(path):
    return st._disk_cache(path, 'txends', lambda: _tx_ends_impl(path))


def compute_rt_extent(bam, gtf, rdna, txends):
    genes, gst, _exons, _est = gtf
    ridna, rdst = rdna

    def in_rdna(ch, s, e):
        v = ridna.get(ch)
        if not v:
            return False
        i = bisect.bisect_right(rdst[ch], e)
        for j in range(max(0, i - 1), -1, -1):
            lo, hi = v[j]
            if hi < s:
                break
            if lo < e and hi > s:
                return True
        return False

    nb = MAXD // BINW + 1; nf = FMAX // FBINW
    hist = [0] * nb; fhist = [0] * nf; dists = []
    n_prim = n_unique = n_gene = 0
    pv = subprocess.Popen(['samtools', 'view', '-F', '0x904', bam], stdout=subprocess.PIPE, text=True)
    for line in pv.stdout:
        n_prim += 1
        if 'NH:i:1\t' not in line and not line.rstrip('\n').endswith('NH:i:1'):
            continue
        n_unique += 1
        f = line.split('\t', 6)
        ch, pos, cig = f[2], int(f[3]), f[5]
        if ch.endswith('_MT'):
            continue
        s, e = pos - 1, st._refend(pos, cig)
        if in_rdna(ch, s, e):
            continue
        gs = st._overlap(genes, gst, ch, s, e, collect=True)
        if not gs:
            continue
        best = max(gs, key=lambda g: min(g[1], e) - max(g[0], s))
        if not (best[2] == 'protein_coding' and best[3]):
            continue
        te = txends.get((ch, best[3]))
        if te is None:
            continue
        n_gene += 1
        strand, ends = te
        if strand == '+':
            anchor = pos
            j = bisect.bisect_left(ends, anchor); d = (ends[j] - anchor) if j < len(ends) else -1
        else:
            anchor = e
            j = bisect.bisect_right(ends, anchor) - 1; d = (anchor - ends[j]) if j >= 0 else -1
        if d < 0:
            continue
        dists.append(d)
        hist[min(d, MAXD) // BINW] += 1
        if d < FMAX:
            fhist[d // FBINW] += 1
    pv.wait()
    n = len(dists)
    # exact distance counts, so a library-level rt_extent can be merged from per-sample
    # ones with exact quartiles (per-sample mode); ~1 key per distinct distance
    from collections import Counter as _Counter
    dist_counts = {str(k): v for k, v in sorted(_Counter(dists).items())}
    if n:
        dists.sort()
        q = lambda p: dists[min(n - 1, int(p * n))]
        stats = {'q25': q(0.25), 'median': q(0.5), 'q75': q(0.75),
                 'pct_lt150': round(100 * sum(1 for d in dists if d < 150) / n, 1),
                 'pct_lt500': round(100 * sum(1 for d in dists if d < 500) / n, 1),
                 'pct_gt2000': round(100 * sum(1 for d in dists if d > 2000) / n, 1)}
    else:
        stats = {'q25': None, 'median': None, 'q75': None, 'pct_lt150': None, 'pct_lt500': None, 'pct_gt2000': None}
    return {'n': n, 'n_primary': n_prim, 'n_unique': n_unique, 'n_in_gene': n_gene,
            'binw': BINW, 'maxd': MAXD, 'hist': hist, 'fbinw': FBINW, 'fmax': FMAX, 'fhist': fhist,
            'dist_counts': dist_counts,
            'anchor': 'aligned end nearest the transcript 5\' end (template-switch site)',
            'scope': 'uniquely mapped primary reads in named protein-coding genes, rRNA and mitochondrial excluded',
            **stats}


if __name__ == '__main__':
    bam = sys.argv[1]
    gtf = st.parse_gtf(st.GTF); rdna = st.load_rdna(st.RDNA_BED); tx = load_tx_ends(st.GTF)
    r = compute_rt_extent(bam, gtf, rdna, tx)
    print({k: v for k, v in r.items() if k not in ('hist', 'fhist')})

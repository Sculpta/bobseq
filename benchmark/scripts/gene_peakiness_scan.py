#!/usr/bin/env python3
"""Candidate genes for the browser examples, ranked by how peaky the prime-seq coverage is (peaks, not intronic read counts), on the pooled
depth-matched BAMs (8.0 M unique fragments per method). Per protein-coding gene (canonical transcript, span 3-150 kb, >= 1.5 kb transcript): per method the depth
along the span; prime-seq peakiness = max depth / mean depth over the canonical exons, the number of sharp peaks (local maxima >= 5x the exon median depth and
narrower than 150 nt at half height), the share of prime-seq depth in its densest 300 nt, its intronic fragment share; BOBseq completeness = share of the transcript
covered, share of positions >= 1/4 of the mean depth, densest-300-nt share, intronic share; DRUG-seq counts. The position of the highest prime-seq peak is recorded
for the sequence check. Output: gene_candidates/all_genes.tsv, list_A_full_coverage.tsv, list_B_primeseq_peaky_intronic.tsv
"""

from settings import WORK, GTF
import os, re, sys, collections, numpy as np, pysam
from concurrent.futures import ProcessPoolExecutor

H = os.path.dirname(os.path.abspath(__file__))
import coverage_lib as cl

P = WORK
B = f'{WORK}/benchmark_uniform/ucsc_pooled'
OUT = f'{WORK}/gene_candidates'
os.makedirs(OUT, exist_ok=True)
METH = {'DRUG-seq': 'drugseq', 'prime-seq': 'primeseq', 'BOBseq': 'bobseq'}


def load():
    ex = collections.defaultdict(list)
    gn = collections.defaultdict(list)
    for l in open(GTF):
        if not l.startswith('HUMAN_'):
            continue
        f = l.split('\t', 9)
        if f[2] == 'exon':
            ex[f[0][6:]].append((int(f[3]) - 1, int(f[4])))
        elif f[2] == 'gene':
            m = re.search(r'gene_name "([^"]+)"', f[8])
            gn[f[0][6:]].append((int(f[3]) - 1, int(f[4]), m.group(1) if m else ''))
    return ex, gn


def peaks(dep, med):
    """sharp peaks: local maxima >= 5x the exon median depth, width at half height < 150 nt"""
    if med <= 0 or len(dep) < 3:
        return 0
    thr = 5 * med
    idx = np.flatnonzero((dep[1:-1] >= thr) & (dep[1:-1] >= dep[:-2]) & (dep[1:-1] > dep[2:])) + 1
    n = 0
    last = -(10**9)
    for i in idx:
        if i - last < 150:
            continue
        h = dep[i] / 2
        l = i
        r = i
        while l > 0 and dep[l - 1] >= h:
            l -= 1
        while r < len(dep) - 1 and dep[r + 1] >= h:
            r += 1
        if r - l < 150:
            n += 1
            last = i
    return n


def chrom_job(c):
    ex, gn = EX[c], GN[c]
    e = np.array(sorted(set(ex)), dtype=np.int64)
    es, ee = e[:, 0], e[:, 1]
    ml = int((ee - es).max())
    rows = []
    bams = {m: pysam.AlignmentFile(f'{B}/{t}_pooled_matched.bam', 'rb') for m, t in METH.items()}
    u = 'chrM' if c == 'MT' else 'chr' + c

    def exonic(blocks):
        for bs, be in blocks:
            i = np.searchsorted(es, be)
            j = np.searchsorted(es, bs - ml)
            if i > j and (ee[j:i] > bs).any():
                return True
        return False

    for g, (L, tid, chrom, strand, starts, ends) in TX.items():
        if chrom != c:
            continue
        s0, e0 = starts[0], ends[-1]
        if not (3000 <= e0 - s0 <= 150000) or L < 1500:
            continue
        inex = np.zeros(e0 - s0, bool)
        for a, b in zip(starts, ends):
            inex[a - s0 : b - s0] = True
        anyex = np.zeros(
            e0 - s0, bool
        )  # any annotated exon of any transcript or gene in the span: intronic depth is measured outside these
        for a, b in zip(es, ee):
            if b > s0 and a < e0:
                anyex[max(a - s0, 0) : min(b - s0, e0 - s0)] = True
        r = {
            'gene': g,
            'chrom': u,
            'start': s0 + 1,
            'end': e0,
            'strand': strand,
            'span_kb': (e0 - s0) / 1e3,
            'tx_len': L,
            'exons': len(starts),
            'other_genes_in_span': sum(1 for a, b, n in gn if a < e0 and b > s0 and n != g),
        }
        for m, bam in bams.items():
            nex = set()
            nin = set()
            d = np.zeros(e0 - s0 + 1, np.int32)
            for a in bam.fetch(u, s0, e0):
                bl = a.get_blocks()
                (nex if exonic(bl) else nin).add(a.query_name)
                for bs, be in bl:
                    d[min(max(bs - s0, 0), e0 - s0)] += 1
                    d[min(max(be - s0, 0), e0 - s0)] -= 1
            nin -= nex
            dep = np.cumsum(d)[: e0 - s0]
            td = dep[inex]
            tot = len(nex) + len(nin)
            r[f'{m}_frag'] = tot
            r[f'{m}_intronic_pct'] = 100 * len(nin) / tot if tot else np.nan
            r[f'{m}_tx_covered_pct'] = 100 * float((td > 0).mean())
            mean = td.mean() if len(td) else 0
            r[f'{m}_even_pct'] = 100 * float((td >= 0.25 * mean).mean()) if mean > 0 else 0.0
            k = np.convolve(dep, np.ones(300), 'valid')
            r[f'{m}_densest300_pct'] = 100 * float(k.max() / max(dep.sum(), 1)) if len(k) else np.nan
            r[f'{m}_peak_over_mean'] = float(dep.max() / mean) if mean > 0 else np.nan
            r[f'{m}_sharp_peaks'] = peaks(dep, float(np.median(td)) if len(td) else 0)
            if m == 'prime-seq':
                i = int(dep.argmax())
                r['primeseq_top_peak_pos'] = s0 + i + 1
                r['primeseq_top_peak_depth'] = int(dep[i])
                r['primeseq_top_peak_in_exon'] = int(inex[i])
            # intron visibility: intronic depth relative to the gene's exonic depth, and how many introns carry visible signal
            idep = dep[~anyex]
            exmax = float(td.max()) if len(td) else 0.0
            r[f'{m}_intron_mean_over_exon_mean'] = float(idep.mean() / mean) if mean > 0 and len(idep) else np.nan
            r[f'{m}_intron_max_over_exon_max'] = float(idep.max() / exmax) if exmax > 0 and len(idep) else np.nan
            r[f'{m}_intron_max_depth'] = int(idep.max()) if len(idep) else 0
            vis = 0
            for a_, b_ in zip(ends[:-1], starts[1:]):
                seg = dep[a_ - s0 : b_ - s0] * (~anyex[a_ - s0 : b_ - s0])
                if len(seg) and exmax > 0 and seg.max() >= 0.2 * exmax:
                    vis += 1
            r[f'{m}_introns_with_visible_signal'] = vis
            r[f'{m}_introns_with_visible_signal_pct'] = 100 * vis / max(len(starts) - 1, 1)
        rows.append(r)
    return rows


EX, GN = load()
TX, _ = cl.build()
if __name__ == '__main__':
    with ProcessPoolExecutor(12) as ex_:
        rows = [r for rs in ex_.map(chrom_job, [str(i) for i in range(1, 23)] + ['X']) for r in rs]
    cols = list(rows[0])
    fmt = lambda v: f'{v:.2f}' if isinstance(v, float) else str(v)

    def write(name, rs):
        with open(f'{OUT}/{name}.tsv', 'w') as f:
            f.write('\t'.join(cols) + '\n')
            [f.write('\t'.join(fmt(r[c]) for c in cols) + '\n') for r in rs]

    write('all_genes', rows)
    A = [
        r
        for r in rows
        if r['BOBseq_frag'] >= 1500
        and r['BOBseq_tx_covered_pct'] >= 97
        and r['BOBseq_even_pct'] >= 75
        and r['BOBseq_intronic_pct'] <= 15
        and r['prime-seq_frag'] >= 200
        and r['DRUG-seq_frag'] >= 100
        and r['other_genes_in_span'] <= 1
    ]
    A.sort(key=lambda r: (-r['BOBseq_even_pct'], r['BOBseq_densest300_pct']))
    write('list_A_full_coverage', A)
    Bl = [
        r
        for r in rows
        if r['BOBseq_frag'] >= 1000
        and r['BOBseq_tx_covered_pct'] >= 95
        and r['BOBseq_even_pct'] >= 65
        and r['BOBseq_intronic_pct'] <= 15
        and r['prime-seq_frag'] >= 300
        and r['prime-seq_intronic_pct'] >= 35
    ]
    Bl.sort(key=lambda r: (-r['prime-seq_sharp_peaks'], -r['prime-seq_peak_over_mean']))
    write('list_B_primeseq_peaky_intronic', Bl)
    C = {}
    C['C1_strict'] = [
        r
        for r in rows
        if r['BOBseq_frag'] >= 1000
        and r['BOBseq_tx_covered_pct'] >= 95
        and r['BOBseq_even_pct'] >= 70
        and r['BOBseq_intronic_pct'] <= 8
        and r['prime-seq_frag'] >= 300
        and r['prime-seq_intronic_pct'] >= 30
        and r['prime-seq_intron_max_over_exon_max'] >= 0.5
        and r['span_kb'] <= 60
    ]
    C['C2_visible_introns'] = [
        r
        for r in rows
        if r['BOBseq_frag'] >= 800
        and r['BOBseq_tx_covered_pct'] >= 93
        and r['BOBseq_even_pct'] >= 60
        and r['BOBseq_intronic_pct'] <= 10
        and r['prime-seq_frag'] >= 200
        and r['prime-seq_introns_with_visible_signal_pct'] >= 50
        and r['prime-seq_intron_max_over_exon_max'] >= 0.35
        and r['span_kb'] <= 80
    ]
    C['C3_intron_depth_ratio'] = [
        r
        for r in rows
        if r['BOBseq_frag'] >= 800
        and r['BOBseq_tx_covered_pct'] >= 93
        and r['BOBseq_intronic_pct'] <= 10
        and r['prime-seq_frag'] >= 150
        and r['prime-seq_intron_mean_over_exon_mean'] >= 0.25
        and r['prime-seq_intron_mean_over_exon_mean']
        >= 4
        * (
            r['BOBseq_intron_mean_over_exon_mean']
            if r['BOBseq_intron_mean_over_exon_mean'] == r['BOBseq_intron_mean_over_exon_mean']
            else 0
        )
        and r['span_kb'] <= 100
    ]
    C['C4_short_genes'] = [
        r
        for r in rows
        if r['span_kb'] <= 25
        and r['BOBseq_frag'] >= 600
        and r['BOBseq_tx_covered_pct'] >= 95
        and r['BOBseq_intronic_pct'] <= 10
        and r['prime-seq_frag'] >= 150
        and r['prime-seq_intronic_pct'] >= 25
        and r['prime-seq_intron_max_over_exon_max'] >= 0.4
    ]
    C['C5_deep_intron_peaks'] = [
        r
        for r in rows
        if '-' not in r['gene']
        and r['other_genes_in_span'] <= 2
        and r['prime-seq_intron_max_depth'] >= 50
        and r['BOBseq_frag'] >= 800
        and r['BOBseq_tx_covered_pct'] >= 93
        and r['BOBseq_intronic_pct'] <= 12
        and r['BOBseq_even_pct'] >= 55
    ]
    C['C5_deep_intron_peaks'].sort(key=lambda r: -r['prime-seq_intron_max_depth'])
    for nm, rs in C.items():
        rs.sort(
            key=lambda r: (
                -r['prime-seq_intron_max_over_exon_max']
                * min(r['prime-seq_introns_with_visible_signal_pct'], 100)
                / 100,
                r['BOBseq_intronic_pct'],
            )
        )
        write(f'list_{nm}', rs)
        print(nm, len(rs), 'genes:', ' '.join(r['gene'] for r in rs[:30]))
    print(
        len(rows),
        'genes scanned |',
        len(A),
        'in list A (clean full coverage) |',
        len(Bl),
        'in list B (prime-seq peaky and intronic, BOBseq full)',
    )
    for name, rs in (('A', A[:25]), ('B', Bl[:25])):
        print(f'--- list {name}')
        for r in rs:
            print(
                f"{r['gene']:10s} {r['chrom']}:{r['start']:,}-{r['end']:,} {r['span_kb']:5.1f}kb {r['exons']:2d}ex | BOB {r['BOBseq_frag']:6d} fr, even {r['BOBseq_even_pct']:3.0f}%, in {r['BOBseq_intronic_pct']:4.1f}% | prime {r['prime-seq_frag']:5d} fr, in {r['prime-seq_intronic_pct']:4.1f}%, peaks {r['prime-seq_sharp_peaks']:2d}, max/mean {r['prime-seq_peak_over_mean']:5.1f}, top300 {r['prime-seq_densest300_pct']:4.0f}% | DRUG {r['DRUG-seq_frag']:5d} fr, peaks {r['DRUG-seq_sharp_peaks']}"
            )

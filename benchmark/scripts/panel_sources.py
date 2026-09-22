#!/usr/bin/env python3
"""Which data each kept panel uses, per figure set, verified from the data itself so that the native and 50-nt sets cannot be
confused: the arm or BAM directory behind every panel, STAR's average input read length from each arm's Log.final.out, the read length
distribution measured on records of each BAM directory, and the max aligned length in the coverage position files. Writes
<fig set>/main_figure_panels/PANEL_SOURCES.md for both sets."""

from settings import WORK, COVERAGE
import os, re, glob, collections, numpy as np, pysam

P = WORK
U = f'{P}/benchmark_uniform'
A = COVERAGE
PB = f'{U}/per_sample_bams'
import os as _os, sys as _sys

_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from palette import METH

ARMS = {
    'native': {
        'DRUG-seq': 'drugseq_native',
        'BRB-seq': 'brbseq_native',
        'prime-seq': 'primeseq_native',
        'BOBseq': 'bob57_24plex_pe_native',
    },
    'matched': {'DRUG-seq': 'drugseq', 'BRB-seq': 'brbseq', 'prime-seq': 'primeseq', 'BOBseq': 'bob57_24plex_pe_hs'},
}
BAMD = {
    'native': {
        'DRUG-seq': 'drugseq_native',
        'BRB-seq': 'brbseq_native',
        'prime-seq': 'primeseq_native',
        'BOBseq': 'bobseq_pe_native',
    },
    'matched': {'DRUG-seq': 'drugseq', 'BRB-seq': 'brbseq', 'prime-seq': 'primeseq', 'BOBseq': 'bobseq_50nt'},
}


def star_len(arm):
    for f in (f'{U}/{arm}/Log.final.out', f'{U}/{arm}/Log.final.selected.out'):
        if os.path.exists(f):
            for l in open(f):
                if 'Average input read length' in l:
                    return l.split('|')[1].strip()
    return 'n/a'


def bam_len(d, n=20000):
    f = sorted(glob.glob(f'{PB}/{d}/*.bam'))[0]
    ls = []
    paired = 0
    r1 = 0
    for i, a in enumerate(pysam.AlignmentFile(f).fetch(until_eof=True)):
        if a.is_unmapped:
            continue
        ls.append(a.query_length or a.infer_query_length() or 0)
        paired += a.is_paired
        r1 += a.is_read1
        if len(ls) >= n:
            break
    ls = np.array(ls)
    return f'median {int(np.median(ls))}, max {ls.max()}, paired {100*paired/len(ls):.0f}% ({os.path.basename(f)})'


def pos_len(sfx, m):
    S = [l.rstrip('\n').split('\t') for l in open(f'{A}/samples{sfx}.tsv') if l.strip()]
    mx = 0
    k = 0
    for lab, bam, mm in S:
        if mm != m:
            continue
        z = np.load(f'{A}/positions_canonical{sfx}/{lab}.npz')
        mx = max(mx, int(z['alen'].max()))
        k += 1
    return f'{k} samples, max aligned length {mx}'


for SET, TAG, sfx in (('native', 'native', ''), ('50nt', 'matched', '_50nt')):
    D = f'{P}/{"preprint_figures_native" if SET == "native" else "preprint_figures_50nt"}/main_figure_panels'
    L = [
        f'# Data sources of the kept panels, {SET} set',
        '',
        '| panel | source | ' + ' | '.join(METH) + ' |',
        '|---|---|' + '---|' * 4,
    ]
    L.append(
        f'| p0a, p1, p2, p3, p3b | benchmark arm (reads.tsv, rarefied.tsv, composition); STAR average input read length | '
        + ' | '.join(f'{ARMS[TAG][m]}: {star_len(ARMS[TAG][m])}' for m in METH)
        + ' |'
    )
    L.append(
        f'| p4i2 / p4m2 threshold panels | coverage positions samples{sfx}.tsv -> positions_canonical{sfx}/ | '
        + ' | '.join(pos_len(sfx, m) for m in METH)
        + ' |'
    )
    L.append(
        f'| p6 / p6b Picard | per-sample BAMs per_sample_bams/<dir>; read length measured on the first 20k mapped records of one BAM | '
        + ' | '.join(f'{BAMD[TAG][m]}: {bam_len(BAMD[TAG][m])}' for m in METH)
        + ' |'
    )
    L += [
        '',
        'Expected: native = DRUG-seq 52 nt, prime-seq 50 nt (its native read), BOBseq 2 x 150 paired (read 2 trimmed to ~134 nt); 50-nt set = every method 50 nt single-end. The BOBseq 50-nt arm is bob57_24plex_pe_hs (read 2 truncated to 50 nt before alignment).',
    ]
    open(f'{D}/PANEL_SOURCES.md', 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))
    print()

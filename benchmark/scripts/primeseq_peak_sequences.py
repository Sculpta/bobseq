#!/usr/bin/env python3
"""What lies under prime-seq's sharp coverage peaks. Peaks = the highest prime-seq position of every scanned gene with peak/mean >= 8, depth >= 30 and at least 1 kb
away from the transcript 3' end (the poly(A) site peak is expected) (gene_peakiness_scan.py), sequence +-60 nt around it on the transcript strand (combined GRCh38 FASTA); background = one random exonic position
of the same genes. Reports base composition by position, A-run content up- and downstream, and the share of peaks with an A-rich stretch (>= 9 A in 12 nt)
within 60 nt. Output: gene_candidates/primeseq_peak_sequences.{tsv,md}"""

from settings import WORK, FASTA
import os, sys, csv, random, numpy as np, pysam

P = WORK
OUT = f'{WORK}/gene_candidates'
FA = pysam.FastaFile(FASTA)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coverage_lib as cl

TX, _ = cl.build()
rng = random.Random(1)
W = 60
RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def seq(c, pos1, strand):
    s = FA.fetch('HUMAN_' + c, max(pos1 - 1 - W, 0), pos1 + W).upper()
    return s.translate(RC)[::-1] if strand == '-' else s


rows = [
    r
    for r in csv.DictReader(open(f'{OUT}/all_genes.tsv'), delimiter='\t')
    if float(r['prime-seq_peak_over_mean']) >= 8 and int(r['primeseq_top_peak_depth']) >= 30
]
pk = []
bg = []
out = []
for r in rows:
    L, tid, c, strand, starts, ends = TX[r['gene']]
    pos = int(r['primeseq_top_peak_pos'])
    end3 = ends[-1] if strand == '+' else starts[0] + 1
    if abs(pos - end3) < 1000:
        continue  # the genuine 3' end (poly(A) site) is not an anomalous peak
    s = seq(c, pos, strand)
    if len(s) != 2 * W + 1:
        continue
    ex = rng.choice(list(zip(starts, ends)))
    b = seq(c, rng.randrange(ex[0] + 1, ex[1] + 1), strand)
    if len(b) != 2 * W + 1:
        continue
    pk.append(s)
    bg.append(b)
    out.append(
        (
            r['gene'],
            r['chrom'],
            r['primeseq_top_peak_pos'],
            strand,
            int(r['primeseq_top_peak_in_exon']),
            r['primeseq_top_peak_depth'],
            s,
        )
    )


def arun(s, side):
    """max number of A in any 12-nt window on one side of the peak (transcript strand)"""
    t = s[W + 1 :] if side == 'down' else s[:W]
    return max(t[i : i + 12].count('A') for i in range(len(t) - 11))


def comp(ss, lo, hi):
    t = ''.join(x[W + lo : W + hi + 1] for x in ss)
    return {b: 100 * t.count(b) / len(t) for b in 'ACGT'}


M = [
    f'# Sequence under the sharpest prime-seq peaks ({len(pk)} peaks, one per gene; background = {len(bg)} random exonic positions of the same genes)',
    '',
    'Sequence on the transcript strand, peak position = 0, windows in nt; A-rich stretch = at least 9 A in a 12-nt window on that side.',
    '',
    '| window | set | A | C | G | T |',
    '|---|---|---|---|---|---|',
]
for lo, hi, lab in (
    (-60, -21, '-60..-21 (upstream)'),
    (-20, -1, '-20..-1'),
    (1, 20, '+1..+20 (downstream)'),
    (21, 60, '+21..+60'),
):
    for nm, ss in (('peaks', pk), ('background', bg)):
        c_ = comp(ss, lo, hi)
        M.append(f"| {lab} | {nm} | {c_['A']:.1f} | {c_['C']:.1f} | {c_['G']:.1f} | {c_['T']:.1f} |")
M += ['', '| A-rich stretch (>= 9 A / 12 nt) | peaks | background |', '|---|---|---|']
for side in ('up', 'down'):
    M.append(
        f"| {side}stream within 60 nt | {100 * np.mean([arun(s, side) >= 9 for s in pk]):.0f}% | {100 * np.mean([arun(s, side) >= 9 for s in bg]):.0f}% |"
    )
M.append(
    f"| T-rich stretch (>= 9 T / 12 nt) either side | {100 * np.mean([max(s[:W][i:i+12].count('T') for i in range(W-11)) >= 9 or max(s[W+1:][i:i+12].count('T') for i in range(W-11)) >= 9 for s in pk]):.0f}% | {100 * np.mean([max(s[:W][i:i+12].count('T') for i in range(W-11)) >= 9 or max(s[W+1:][i:i+12].count('T') for i in range(W-11)) >= 9 for s in bg]):.0f}% |"
)
M += ['', f"Peaks inside an annotated exon: {100 * np.mean([o[4] for o in out]):.0f}%.", '']
open(f'{OUT}/primeseq_peak_sequences.md', 'w').write('\n'.join(M) + '\n')
with open(f'{OUT}/primeseq_peak_sequences.tsv', 'w') as f:
    f.write(
        'gene\tchrom\tpeak_pos_hg38\tstrand\tpeak_in_exon\tpeak_depth\tsequence_transcript_strand_minus60_to_plus60\n'
    )
    [f.write('\t'.join(map(str, o)) + '\n') for o in out]
print('\n'.join(M))

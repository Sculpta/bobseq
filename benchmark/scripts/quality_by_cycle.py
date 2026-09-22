#!/usr/bin/env python3
"""quality_by_cycle.py <label> <out.json> [--reads N] < FASTQ   (uncompressed on stdin; default 200,000 reads)
Per sequencing cycle: mean Phred quality, fraction of bases with Q >= 30, fraction of bases called N, and the base composition
(A/C/G/T share), so a poly(dT) or poly(A) homopolymer and the quality collapse after it are both visible (both reads of a
PE150 run are usable for BOBseq, whereas the poly(dT) methods lose quality after the poly(dT) region)."""

import sys, json, numpy as np

lab, out = sys.argv[1:3]
N = int(sys.argv[sys.argv.index('--reads') + 1]) if '--reads' in sys.argv else 200000
S = C = Q30 = NN = None
B = None
n = 0
seq = None
for i, line in enumerate(sys.stdin):
    r = i % 4
    if r == 1:
        seq = line.rstrip('\n')
    elif r == 3:
        q = np.frombuffer(line.rstrip('\n').encode(), dtype=np.uint8).astype(np.int16) - 33
        s = np.frombuffer(seq.encode(), dtype=np.uint8)
        if S is None:
            L = len(q)
            S = np.zeros(L)
            C = np.zeros(L)
            Q30 = np.zeros(L)
            NN = np.zeros(L)
            B = {b: np.zeros(L) for b in 'ACGT'}
        k = min(len(q), L)
        S[:k] += q[:k]
        C[:k] += 1
        Q30[:k] += q[:k] >= 30
        NN[:k] += s[:k] == ord('N')
        for b in 'ACGT':
            B[b][:k] += s[:k] == ord(b)
        n += 1
        if n >= N:
            break
C = np.maximum(C, 1)
res = {
    'label': lab,
    'reads': n,
    'cycles': int(len(S)),
    'mean_q': (S / C).round(3).tolist(),
    'frac_q30': (Q30 / C).round(4).tolist(),
    'frac_N': (NN / C).round(4).tolist(),
    'base_share': {b: (B[b] / C).round(4).tolist() for b in 'ACGT'},
}
json.dump(res, open(out, 'w'))
m = S / C
print(
    lab,
    f'{n:,} reads, {len(S)} cycles; mean Q at cycle 1/25/50/75/100/125/150:',
    ' '.join(f'{m[c-1]:.0f}' for c in (1, 25, 50, 75, 100, 125, 150) if c <= len(m)),
)

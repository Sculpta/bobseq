#!/usr/bin/env python3
"""build_positions.py <label> <bam> <out.npz> [--max N]: unique (MAPQ 255) primary alignments -> mature-transcript position of each
read 5' end on the longest protein-coding transcript of its gene. Stores gene id (int), t, L, mate (1/2/0), sense (read strand == gene strand).
"""

import sys, os, numpy as np, pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coverage_lib as cl

lab, bam, out = sys.argv[1:4]
mx = int(sys.argv[sys.argv.index('--max') + 1]) if '--max' in sys.argv else None
tx, idx = cl.build()
genes = sorted(tx)
gid = {g: i for i, g in enumerate(genes)}
G, T, L, M, S, AL = [], [], [], [], [], []
n = n_pos = 0
for a in pysam.AlignmentFile(bam, 'rb').fetch(until_eof=True):
    if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality != 255:
        continue
    n += 1
    if mx and n > mx:
        break
    p5 = a.reference_end - 1 if a.is_reverse else a.reference_start
    r = cl.position(idx, a.reference_name, p5)
    if r is None:
        continue
    g, t, Lg, gs = r
    n_pos += 1
    G.append(gid[g])
    T.append(t)
    L.append(Lg)
    M.append(1 if a.is_read1 else (2 if a.is_read2 else 0))
    S.append((a.is_reverse and gs == '-') or (not a.is_reverse and gs == '+'))
    AL.append(
        sum(ln for op, ln in a.cigartuples if op in (0, 7, 8))
    )  # aligned bases = span on the mature transcript (exon-aware by construction)
np.savez_compressed(
    out,
    gene=np.array(G, np.int32),
    t=np.array(T, np.int32),
    L=np.array(L, np.int32),
    mate=np.array(M, np.int8),
    sense=np.array(S, bool),
    alen=np.array(AL, np.int16),
    n_unique=n,
    genes=np.array(genes),
    tx_len=np.array([tx[g][0] for g in genes], np.int32),
)
print(
    f'{lab}: {n:,} unique reads, {n_pos:,} positioned on canonical protein-coding transcripts ({100*n_pos/max(n,1):.1f}%)'
)

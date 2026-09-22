#!/usr/bin/env python3
"""junctions_per_sample.py <label> <bam> <out.tsv>: junctions per sample at matched depth. Unique (MAPQ 255) primary reads are reservoir-
sampled to the largest depth, shuffled, and for each depth the first d reads are scanned: spliced reads (N in CIGAR), distinct junctions
(contig, start, end) seen >= 1 and >= 2 times. BOBseq: both mates count as reads (a junction spanned by both mates of one pair counts twice
for the 'reads' but once for 'distinct junctions'); seed 1."""

import sys, random, collections, pysam

lab, bam, out = sys.argv[1:4]
DEPTHS = [100_000, 250_000, 500_000, 1_000_000, 2_000_000]
K = DEPTHS[-1]
rng = random.Random(1)
res = []
n = 0
for a in pysam.AlignmentFile(bam, 'rb').fetch(until_eof=True):
    if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality != 255:
        continue
    n += 1
    js = []
    if 'N' in a.cigarstring:
        pos = a.reference_start
        for op, ln in a.cigartuples:
            if op == 3:
                js.append((a.reference_name, pos, pos + ln))
            if op in (0, 2, 3, 7, 8):
                pos += ln
    rec = tuple(js)
    if len(res) < K:
        res.append(rec)
    else:
        j = rng.randrange(n)
        if j < K:
            res[j] = rec
rng.shuffle(res)
with open(out, 'w') as f:
    f.write(
        'sample\tdepth\tunique_reads_total\treads_used\tspliced_reads\tjunctions_ge1\tjunctions_ge2\tjunctions_ge3\tjunctions_ge5\tat_native\n'
    )
    for d in DEPTHS:
        use = res[:d] if n >= d else res
        cnt = collections.Counter(j for r in use for j in r)
        spl = sum(1 for r in use if r)
        f.write(
            f'{lab}\t{d}\t{n}\t{len(use)}\t{spl}\t{len(cnt)}\t{sum(1 for v in cnt.values() if v >= 2)}\t{sum(1 for v in cnt.values() if v >= 3)}\t{sum(1 for v in cnt.values() if v >= 5)}\t{int(n < d)}\n'
        )
        if n < d:
            break
print(f'{lab}: {n:,} unique reads')

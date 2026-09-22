#!/usr/bin/env python3
"""Are prime-seq's sharp peaks PCR duplicates? The pooled tracks are raw (pre-deduplication) alignments and prime-seq read names carry sample + 10-nt UMI.
For the peaks of primeseq_peak_sequences.tsv (sharp, >= 1 kb from the 3' end): reads whose 5' end lies within 25 nt of the peak -> distinct (sample, UMI) pairs
/ reads = molecules per read; the same for reads in a control window 2 kb away in the same gene. Output: gene_candidates/primeseq_peak_umis.md
"""

from settings import WORK
import os, csv, collections, numpy as np, pysam

P = WORK
OUT = f'{WORK}/gene_candidates'
b = pysam.AlignmentFile(f'{WORK}/benchmark_uniform/ucsc_pooled/primeseq_pooled_matched.bam', 'rb')


def stats(chrom, pos):
    reads = 0
    mols = set()
    samples = collections.Counter()
    for a in b.fetch(chrom, max(pos - 200, 0), pos + 200):
        p5 = a.reference_end - 1 if a.is_reverse else a.reference_start
        if abs(p5 - (pos - 1)) > 25:
            continue
        reads += 1
        s, u = a.query_name.split(':', 1)[1].rsplit('_', 2)[1:]
        mols.add((s, u))
        samples[s] += 1
    return reads, len(mols), samples


rows = list(csv.DictReader(open(f'{OUT}/primeseq_peak_sequences.tsv'), delimiter='\t'))
res = []
for r in rows:
    pos = int(r['peak_pos_hg38'])
    pr, pm, ps = stats(r['chrom'], pos)
    if pr < 20:
        continue
    cr, cm, _ = stats(r['chrom'], pos + 2000)
    res.append((r['gene'], pr, pm, pm / pr, ps.most_common(1)[0][1] / pr, cr, cm, cm / cr if cr else np.nan))
mr = np.array([x[3] for x in res])
ms = np.array([x[4] for x in res])
cr = np.array([x[7] for x in res if not np.isnan(x[7]) and x[5] >= 20])
M = [
    f'# UMI diversity under the sharp prime-seq peaks ({len(res)} peaks with >= 20 reads within 25 nt)',
    '',
    f'- molecules per read at the peak (distinct sample + UMI / reads): median {np.median(mr):.2f}, quartiles {np.percentile(mr, 25):.2f} to {np.percentile(mr, 75):.2f}',
    f'- control windows 2 kb away with >= 20 reads (n={len(cr)}): median {np.median(cr):.2f}, quartiles {np.percentile(cr, 25):.2f} to {np.percentile(cr, 75):.2f}',
    f'- share of the peak reads from its single largest sample: median {100 * np.median(ms):.0f}%; peaks where one sample gives >= 80% of the reads: {100 * np.mean(ms >= 0.8):.0f}%',
    f'- peaks that are mostly duplicates (< 0.2 molecules per read): {100 * np.mean(mr < 0.2):.0f}%; peaks that are mostly distinct molecules (> 0.8): {100 * np.mean(mr > 0.8):.0f}%',
    '',
    '| gene | reads at peak | molecules | molecules / read | largest sample share | control reads | control molecules / read |',
    '|---|---|---|---|---|---|---|',
]
for x in sorted(res, key=lambda x: -x[1])[:15]:
    M.append(
        f'| {x[0]} | {x[1]} | {x[2]} | {x[3]:.2f} | {100 * x[4]:.0f}% | {x[5]} | {"na" if np.isnan(x[7]) else f"{x[7]:.2f}"} |'
    )
open(f'{OUT}/primeseq_peak_umis.md', 'w').write('\n'.join(M) + '\n')
print('\n'.join(M[:7]))

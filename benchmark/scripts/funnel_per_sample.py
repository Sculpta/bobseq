#!/usr/bin/env python3
"""funnel_per_sample.py <arm_dir> : per-sample read fate from the pooled STAR BAM (unmapped reads kept, sample tag = 2nd-last '_' field of QNAME).
Counts fragments (read 1 only for paired data; secondary/supplementary skipped): total, mapped, uniquely mapped (MAPQ 255).
Writes <arm_dir>/funnel_per_sample.tsv. The 'usable' step comes from stats_E_U/rarefied.tsv (reads_native) in read_fate_panels.py.
"""

import sys, collections, pysam

arm = sys.argv[1].rstrip('/')
tot = collections.Counter()
mp = collections.Counter()
uq = collections.Counter()
n = 0
for a in pysam.AlignmentFile(f'{arm}/Aligned.out.bam', 'rb').fetch(until_eof=True):
    if a.is_secondary or a.is_supplementary or (a.is_paired and not a.is_read1):
        continue
    s = a.query_name.rsplit('_', 2)[1]
    tot[s] += 1
    n += 1
    if not a.is_unmapped:
        mp[s] += 1
        if a.mapping_quality == 255:
            uq[s] += 1
with open(f'{arm}/funnel_per_sample.tsv', 'w') as f:
    f.write('sample\treads_in\tmapped\tunique\n')
    for s in sorted(tot):
        f.write(f'{s}\t{tot[s]}\t{mp[s]}\t{uq[s]}\n')
print(f'{arm}: {n:,} fragments, {len(tot)} samples')

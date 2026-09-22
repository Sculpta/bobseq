#!/usr/bin/env python3
"""Mate-1 records of one raw (pre-deduplication) pair BAM of the pipeline, in the benchmark's read-name convention.
    native_arm_mate1.py <sample_star_pairs.bam> <out.bam>
Keeps mate 1 (read 2 of the sequencer: the bobcode side) primary mapped records (flag & 0x904 == 0), one record per fragment,
and rewrites the pipeline's read name '<id>__<umi>_<bobcode>_<grun>_<rt>' to the benchmark's '<id>_<bobcode>_<umi>'
(read_table.py: well, umi = rsplit('_', 2)[1:])."""
import sys

import pysam

src, out = sys.argv[1:3]
b = pysam.AlignmentFile(src, 'rb')
o = pysam.AlignmentFile(out, 'wb', template=b)
n = kept = 0
for a in b.fetch(until_eof=True):
    n += 1
    if not a.is_read1 or a.flag & 0x904:
        continue
    rid, rest = a.query_name.split('__', 1)
    umi, bobcode = rest.split('_')[:2]
    a.query_name = f'{rid}_{bobcode}_{umi}'
    o.write(a)
    kept += 1
o.close()
print(f'{src}: {n:,} records, {kept:,} mate-1 primary mapped records written')

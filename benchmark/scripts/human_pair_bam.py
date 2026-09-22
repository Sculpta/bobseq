#!/usr/bin/env python3
"""Per-sample human BAM from one raw pair BAM of the pipeline: primary mapped records of both mates (flag & 0x904 == 0) on the
human contigs, with the HUMAN_ prefix stripped from the contig names (Ensembl names: 1..22, X, Y, MT, scaffolds), coordinate-sorted
and indexed. Read names are left as the pipeline wrote them.
    human_pair_bam.py <sample_star_pairs.bam> <out.bam>
"""
import os
import sys

import pysam

src, out = sys.argv[1:3]
b = pysam.AlignmentFile(src, 'rb')
refs = [(i, n) for i, n in enumerate(b.references) if n.startswith('HUMAN_')]
remap = {i: k for k, (i, n) in enumerate(refs)}
hdr = {'HD': {'VN': '1.4', 'SO': 'unsorted'}, 'SQ': [{'SN': n[6:], 'LN': b.lengths[i]} for i, n in refs]}
tmp = out + '.unsorted.bam'
o = pysam.AlignmentFile(tmp, 'wb', header=hdr)
n = kept = 0
for a in b.fetch(until_eof=True):
    n += 1
    if a.flag & 0x904 or a.reference_id not in remap:
        continue
    d = a.to_dict()
    d['ref_name'] = a.reference_name[6:]
    nr = a.next_reference_name
    d['next_ref_name'] = '=' if a.next_reference_id == a.reference_id else (nr[6:] if nr and nr.startswith('HUMAN_') else '*')
    if d['next_ref_name'] == '*':
        d['next_ref_pos'] = '0'
        d['length'] = '0'
    o.write(pysam.AlignedSegment.from_dict(d, o.header))
    kept += 1
o.close()
pysam.sort('-@', '4', '-o', out, tmp)
pysam.index(out)
os.remove(tmp)
print(f'{os.path.basename(out)}: {n:,} records, {kept:,} primary mapped records on human contigs')

#!/usr/bin/env python3
"""Split one pooled benchmark BAM (50-nt common path, combined GRCh38+GRCm39 reference, sample in the QNAME as
<readid>_<SAMPLE>_<UMI>) into per-sample HUMAN-only BAMs with Ensembl contig names (HUMAN_1 -> 1, HUMAN_MT -> MT,
scaffolds kept with the prefix stripped), primary mapped records only. Streams the BAM once.
    split_pooled.py <bam> <outdir> <method> <samples.tsv: key \t sample_name>
"""

import subprocess, sys, os

bam, out, method, smap = sys.argv[1:5]
names = dict(l.rstrip('\n').split('\t') for l in open(smap) if l.strip())
os.makedirs(out, exist_ok=True)
hdr_lines = subprocess.run(
    ['samtools', 'view', '-H', bam], capture_output=True, text=True, check=True
).stdout.splitlines()
header = []
for h in hdr_lines:
    if h.startswith('@SQ'):
        f = h.split('\t')
        if f[1].startswith('SN:HUMAN_'):
            header.append('\t'.join([f[0], 'SN:' + f[1][len('SN:HUMAN_') :]] + f[2:]))
    elif not h.startswith('@PG'):
        header.append(h)
header.append(
    f'@PG\tID:split_pooled\tPN:split_pooled.py\tDS:{method}: human primary alignments of the benchmark common path (50 nt), one sample'
)
header = '\n'.join(header) + '\n'
writers = {}


def w(sample):
    if sample not in writers:
        p = subprocess.Popen(
            ['samtools', 'view', '-b', '-o', f'{out}/{sample}.unsorted.bam', '-'], stdin=subprocess.PIPE, text=True
        )
        p.stdin.write(header)
        writers[sample] = p
    return writers[sample].stdin


n = kept = 0
per = {}
rd = subprocess.Popen(['samtools', 'view', '-F', '0x904', bam], stdout=subprocess.PIPE, text=True)
for line in rd.stdout:
    n += 1
    f = line.split('\t', 3)
    if not f[2].startswith('HUMAN_'):
        continue
    key = f[0].rsplit('_', 2)[1]
    s = names.get(key)
    if s is None:
        continue
    w(s).write(f[0] + '\t' + f[1] + '\t' + f[2][6:] + '\t' + f[3])
    kept += 1
    per[s] = per.get(s, 0) + 1
rd.wait()
for p in writers.values():
    p.stdin.close()
for p in writers.values():
    p.wait()
for s, c in sorted(per.items()):
    print(f'{method}\t{s}\t{c}')
print(f'# {method}: {n:,} primary mapped records streamed, {kept:,} kept in {len(per)} samples', file=sys.stderr)

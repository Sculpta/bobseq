#!/usr/bin/env python3
"""Well-barcode error correction for the DRUG-seq arm: every 10-nt sequence within one substitution of exactly one whitelist
barcode (config/drugseq_barcodes_384.tsv) maps to that barcode; the barcodes map to themselves. slice.sh applies the map to
the observed barcode in the read name and drops reads whose barcode is not in it. Also writes wells_keep.txt: the 24 DMSO
wells (config/drugseq_dmso_wells.txt), the benchmark samples.
    drugseq_well_map.py <arm_dir>
"""
import os
import sys

from settings import CONFIG

arm = sys.argv[1]
whitelist = [l.split('\t')[1].strip() for l in open(f'{CONFIG}/drugseq_barcodes_384.tsv').read().splitlines()[1:] if l.strip()]
dmso = [l.strip() for l in open(f'{CONFIG}/drugseq_dmso_wells.txt') if l.strip()]
assert len(whitelist) == 384 and len(set(whitelist)) == 384 and len(dmso) == 24 and set(dmso) <= set(whitelist)
wmap = {}
ambiguous = set()
for w in whitelist:
    for i in range(len(w)):
        for b in 'ACGT':
            if b == w[i]:
                continue
            v = w[:i] + b + w[i + 1:]
            if v in wmap and wmap[v] != w:
                ambiguous.add(v)
            wmap[v] = w
for v in ambiguous:
    del wmap[v]
for w in whitelist:
    wmap[w] = w
os.makedirs(arm, exist_ok=True)
with open(f'{arm}/well_map.tsv', 'w') as f:
    for v in sorted(wmap):
        f.write(f'{v}\t{wmap[v]}\n')
with open(f'{arm}/wells_keep.txt', 'w') as f:
    f.write('\n'.join(dmso) + '\n')
print(f'{arm}: well_map.tsv {len(wmap):,} sequences -> 384 barcodes ({len(ambiguous)} ambiguous dropped); wells_keep.txt {len(dmso)} DMSO wells')

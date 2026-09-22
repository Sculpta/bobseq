#!/usr/bin/env python3
"""Annotation index for the strict barcode-accuracy funnel (the preprint's read-classification rules).
From the combined GRCh38 + GRCm39 Ensembl 113 GTF: per contig, sorted exon intervals with gene biotype and name (all transcripts),
and the rDNA loci bed. Written to inputs/annotation_index.npz (+ gene lists) so the per-well script needs no GTF parsing.
"""

from settings import REF_DIR
import re, sys, os, collections, numpy as np

REF = REF_DIR
OUT = sys.argv[1] if len(sys.argv) > 1 else 'inputs/annotation_index.npz'
RIBO = re.compile(r'^(RP[LS]|RPLP|MRP[LS]|Rp[ls]|Rplp|Mrp[ls])[^K]*$', re.I)
RRNA_BT = {'rRNA', 'Mt_rRNA', 'rRNA_pseudogene'}
bt_re = re.compile(r'gene_biotype "([^"]+)"')
nm_re = re.compile(r'gene_name "([^"]+)"')
ex = collections.defaultdict(list)
nrp = set()
nrr = set()
for l in open(f'{REF}/combined_genome.gtf'):
    if l[0] == '#':
        continue
    f = l.split('\t', 9)
    if f[2] != 'exon':
        continue
    bt = bt_re.search(f[8])
    nm = nm_re.search(f[8])
    bt = bt.group(1) if bt else ''
    nm = nm.group(1) if nm else ''
    cls = (
        2
        if bt in RRNA_BT
        else (1 if (bt == 'protein_coding' and RIBO.match(nm)) else (0 if bt == 'protein_coding' else 3))
    )  # 0 mRNA, 1 ribosomal-protein, 2 rRNA biotype, 3 other biotype
    ex[f[0]].append((int(f[3]) - 1, int(f[4]), cls))
    if cls == 1:
        nrp.add(nm)
    if cls == 2:
        nrr.add(nm)
arr = {}
for c, v in ex.items():
    v = sorted(set(v))
    a = np.array(v, dtype=np.int64)
    arr[f'{c}__s'] = a[:, 0]
    arr[f'{c}__e'] = a[:, 1]
    arr[f'{c}__c'] = a[:, 2].astype(np.int8)
    arr[f'{c}__maxlen'] = np.array([int((a[:, 1] - a[:, 0]).max())])
rd = collections.defaultdict(list)
for l in open(f'{REF}/rdna_loci.bed'):
    if l[0] == '#' or l.startswith(('track', 'browser')) or not l.strip():
        continue
    c, s, e = l.split()[:3]
    rd[c].append((int(s), int(e)))
for c, v in rd.items():
    a = np.array(sorted(v), dtype=np.int64)
    arr[f'rdna__{c}__s'] = a[:, 0]
    arr[f'rdna__{c}__e'] = a[:, 1]
np.savez_compressed(OUT, **arr)
print(
    'contigs with exons',
    len(ex),
    '| exon intervals',
    sum(len(v) for v in ex.values()),
    '| ribosomal-protein genes',
    len(nrp),
    '| rRNA-biotype genes',
    len(nrr),
    '| rDNA loci',
    sum(len(v) for v in rd.values()),
)
open(os.path.join(os.path.dirname(OUT) or '.', 'ribosomal_protein_genes.txt'), 'w').write('\n'.join(sorted(nrp)) + '\n')

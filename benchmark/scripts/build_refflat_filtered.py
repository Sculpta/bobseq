#!/usr/bin/env python3
"""Filtered refFlat for Picard CollectRnaSeqMetrics: protein-coding Ensembl 113 transcripts (transcript biotype
protein_coding) of the human genome, with rRNA (not protein-coding anyway), mitochondrial (chromosome MT) and ribosomal-protein genes
(RPL/RPS/MRPL/MRPS/RPLP, RPS6K kinases spared) excluded. Chromosome names are BARE (1, 2, ..., X) because the benchmark's per-sample
BAMs (per_sample_bams) carry Ensembl names; a contig mismatch makes Picard report an empty histogram, not an error."""

from settings import GTF
import re, sys, collections

PREFIX = sys.argv[1] if len(sys.argv) > 1 else 'HUMAN_'  # HUMAN_ or MOUSE_ (both 6 characters)
OUT = sys.argv[2] if len(sys.argv) > 2 else 'refflat_pc_noRP_noMT_bare.txt'
RIBO = re.compile(r'^(RP[LS]|RPLP|MRP[LS])[^K]*$', re.I)
attr = re.compile(r'(\S+) "([^"]*)"')
tx = {}
ex = collections.defaultdict(list)
cds = collections.defaultdict(list)
n_rp = n_mt = 0
for l in open(GTF):
    if l[0] == '#' or not l.startswith(PREFIX):
        continue
    f = l.rstrip('\n').split('\t')
    if f[2] not in ('exon', 'CDS'):
        continue
    a = dict(attr.findall(f[8]))
    if a.get('transcript_biotype') != 'protein_coding':
        continue
    chrom = f[0][6:]
    if chrom == 'MT':
        n_mt += 1
        continue
    g = a.get('gene_name', a.get('gene_id', ''))
    t = a['transcript_id']
    if RIBO.match(g):
        n_rp += 1
        continue
    tx[t] = (g, chrom, f[6])
    (ex if f[2] == 'exon' else cds)[t].append((int(f[3]) - 1, int(f[4])))
out = open(OUT, 'w')
n = 0
genes = set()
for t, (g, chrom, strand) in tx.items():
    e = sorted(ex[t])
    if not e:
        continue
    c = sorted(cds[t])
    cs, ce = (c[0][0], c[-1][1]) if c else (e[-1][1], e[-1][1])
    out.write(
        f'{g}\t{t}\t{chrom}\t{strand}\t{e[0][0]}\t{e[-1][1]}\t{cs}\t{ce}\t{len(e)}\t{",".join(str(s) for s, _ in e)},\t{",".join(str(x) for _, x in e)},\n'
    )
    n += 1
    genes.add(g)
out.close()
print(
    f'refFlat: {n} protein-coding transcripts, {len(genes)} genes; excluded exon/CDS lines: ribosomal-protein {n_rp}, MT {n_mt}'
)

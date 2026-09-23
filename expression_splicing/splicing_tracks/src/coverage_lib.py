"""Transcript models from the Ensembl 113 GTF (adapted from the benchmark coverage_lib.py: the GTF alone is the source --
no refFlat, no species prefix -- and the gene set is the handful of loci asked for, any biotype).
Contig names are Ensembl (1..22, X, Y, MT) because the BAMs use them.

build(genes) -> tx, exons
  tx[gene]    = (L, tid, chrom, strand, starts, ends)   the Ensembl canonical transcript (MANE Select where one exists; the
                longest transcript for genes without the tag), exons 0-based half-open, ascending; L = mature length --
                the tuple layout of the benchmark script, which draws it as the transcript row
  exons[gene] = {(start, end): [transcript ids]}        every annotated exon of the gene, any transcript -- used to place
                the zoom window on whole exons (an LSV site may belong to a non-canonical exon: the NMD poison exons)
"""
from settings import GTF, open_text
import re


def _attr(a, key):
    i = a.find(key + ' "')
    return a[i + len(key) + 2 : a.find('"', i + len(key) + 2)] if i >= 0 else ''


def build(genes):
    genes = set(genes)
    tx_exons, tx_meta, exons = {}, {}, {}
    for l in open_text(GTF):
        if l[0] == '#':
            continue
        f = l.split('\t', 9)
        if f[2] not in ('transcript', 'exon'):
            continue
        a = f[8]
        gname, gid = _attr(a, 'gene_name'), _attr(a, 'gene_id')
        gene = gname if gname in genes else gid if gid in genes else None
        if gene is None:
            continue
        tid = _attr(a, 'transcript_id')
        if f[2] == 'transcript':
            tx_meta[tid] = (gene, f[0], f[6], 'Ensembl_canonical' in a, _attr(a, 'transcript_biotype'))
        else:
            s, e = int(f[3]) - 1, int(f[4])
            tx_exons.setdefault(tid, []).append((s, e))
            exons.setdefault(gene, {}).setdefault((s, e), []).append(tid)
    tx = {}
    rank = {}
    for tid, (gene, chrom, strand, canonical, biotype) in tx_meta.items():
        ex = sorted(tx_exons.get(tid, []))
        if not ex:
            continue
        L = sum(e - s for s, e in ex)
        r = (1 if canonical else 0, L)
        if gene not in tx or r > rank[gene]:
            tx[gene] = (L, tid, chrom, strand, [s for s, _ in ex], [e for _, e in ex])
            rank[gene] = r
    missing = genes - set(tx)
    if missing:
        raise SystemExit(f'not in the GTF: {sorted(missing)}')
    return tx, exons


def transcript_biotypes():
    """transcript id -> transcript_biotype, for labelling the event exon's transcript (e.g. nonsense_mediated_decay)."""
    out = {}
    for l in open_text(GTF):
        if l[0] == '#':
            continue
        f = l.split('\t', 9)
        if f[2] == 'transcript':
            out[_attr(f[8], 'transcript_id')] = (_attr(f[8], 'transcript_name'), _attr(f[8], 'transcript_biotype'))
    return out

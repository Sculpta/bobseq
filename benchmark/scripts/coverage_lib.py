"""Transcript models (longest transcript per human gene, Ensembl 113) and the mature-transcript position of a read 5' end.
Contig names are Ensembl (1..22, X, Y, MT) because the per-sample BAMs use them.
"""

from settings import GTF, REFFLAT
import bisect, collections, os

# species switch for the mouse wells of the 24-plex: BM_SPECIES_PREFIX=MOUSE_ ; default = human. Both prefixes are 6 characters.
PREFIX = os.environ.get('BM_SPECIES_PREFIX', 'HUMAN_')
assert PREFIX in ('HUMAN_', 'MOUSE_')
REF = REFFLAT


def protein_coding_genes():
    pc = set()
    for l in open(GTF):
        if l[0] == '#' or not l.startswith(PREFIX):
            continue
        f = l.split('\t', 9)
        if f[2] != 'gene' or 'gene_biotype "protein_coding"' not in f[8]:
            continue
        i = f[8].find('gene_name "')
        pc.add(f[8][i + 11 : f[8].find('"', i + 11)])
    return pc


def canonical_transcripts():
    """Ensembl canonical transcript per human gene (MANE Select where one exists), from the GTF tags."""
    can = {}
    for l in open(GTF):
        if l[0] == '#' or not l.startswith(PREFIX):
            continue
        f = l.split('\t', 9)
        if f[2] != 'transcript' or 'Ensembl_canonical' not in f[8]:
            continue
        i = f[8].find('gene_name "')
        g = f[8][i + 11 : f[8].find('"', i + 11)]
        j = f[8].find('transcript_id "')
        can[g] = f[8][j + 15 : f[8].find('"', j + 15)]
    return can


def build(protein_coding_only=True, model='canonical'):
    """model = 'canonical' (Ensembl canonical / MANE Select transcript per gene; the longest is used only for genes without one) or 'longest'.
    'longest' systematically picks isoforms with extended 3' UTRs (ACTB: 3' end 738 nt beyond the used poly(A) site), which
    shifts every distance-from-3'-end measure for every method; canonical avoids that."""
    pc = protein_coding_genes() if protein_coding_only else None
    can = canonical_transcripts() if model == 'canonical' else {}
    tx = {}
    rank = {}
    for l in open(REF):
        f = l.rstrip('\n').split('\t')
        gene, tid, chrom, strand = f[0], f[1], f[2], f[3]
        if not chrom.startswith(PREFIX) or (pc is not None and gene not in pc):
            continue
        starts = [int(x) for x in f[9].rstrip(',').split(',')]
        ends = [int(x) for x in f[10].rstrip(',').split(',')]
        L = sum(e - s for s, e in zip(starts, ends))
        r = (1 if can.get(gene) == tid else 0, L)
        if gene not in tx or r > rank[gene]:
            tx[gene] = (L, tid, chrom[6:], strand, starts, ends)
            rank[gene] = r
    by_chr = collections.defaultdict(list)
    for gene, (L, tid, chrom, strand, starts, ends) in tx.items():
        off = 0
        for s, e in zip(starts, ends):
            by_chr[chrom].append((s, e, gene, off, strand, L))
            off += e - s
    idx = {}
    for chrom, iv in by_chr.items():
        iv.sort()
        idx[chrom] = ([x[0] for x in iv], iv)
    return tx, idx


def position(idx, chrom, pos5):
    """(gene, t, L, gene_strand): t = distance from the transcript 5' end along the mature transcript (exons only)."""
    if chrom not in idx:
        return None
    starts, iv = idx[chrom]
    i = bisect.bisect_right(starts, pos5) - 1
    for j in range(i, max(i - 6, -1), -1):
        s, e, gene, off, strand, L = iv[j]
        if s <= pos5 < e:
            t = off + (pos5 - s) if strand == '+' else L - 1 - (off + (pos5 - s))
            return gene, t, L, strand
    return None

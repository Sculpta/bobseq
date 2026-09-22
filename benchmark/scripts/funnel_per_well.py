#!/usr/bin/env python3
"""funnel_per_well.py <sample> <bam> <annotation_index.npz> <exact_code_keys.npz> <out.json> [--reads]
Strict barcode-accuracy funnel for one well of the 24-plex Illumina run, following the preprint's nanopore rules step by step
(Illumina-adapted: the bobcode is positional, read at prep from the first bases of read 2; only mate 1 = trimmed read 2 is scored,
which is what the pipeline's deduplicated analysis BAM holds; a second barcode unit is not observable in a 150-nt mate).
Levels (cumulative), accuracy after each:
  L0  mapped molecules carrying a bobcode (species from the primary contig)
  L1  + species called from ALL alignments (human or mouse only; alignments on both genomes = ambiguous, not scored)
  L2  + not rRNA (any alignment on an rDNA locus, or the primary alignment on an exon of an rRNA-biotype gene) and not mitochondrial
  L3  + mRNA: the primary alignment overlaps an exon of a protein-coding gene that is not a ribosomal-protein gene
      (ribosomal-protein, non-coding, intronic and intergenic = other, not scored)
  L4  + exact bobcode (zero mismatches; from exact_code_flags.py, the pipeline accepts one mismatch at prep)
No low-complexity (entropy) level is applied.
Correct = the aligned species equals the species of the well's bobcode. Wilson 95% intervals. --reads: score every read pair
(read-1 records) of a raw pair BAM instead of deduplicated molecules."""

from settings import RUN_JSON
import sys, json, math, collections, numpy as np, pysam

sample, bam, annot, exact, out = sys.argv[1:6]
READS = '--reads' in sys.argv
A = np.load(annot)
EX = {}
RD = {}
for k in A.files:
    if k.startswith('rdna__'):
        c = k.split('__')[1]
        RD.setdefault(c, {})[k.split('__')[2]] = A[k]
    elif k.endswith('__s'):
        c = k[:-3]
        EX[c] = (A[f'{c}__s'], A[f'{c}__e'], A[f'{c}__c'], int(A[f'{c}__maxlen'][0]))
KEYS = np.load(exact)['keys']
RUN = json.load(open(RUN_JSON))
SPM = RUN['tso_species_map']


def key_of(rid):
    p = rid.split(':')
    return (int(p[3]) << 56) | (int(p[4]) << 40) | (int(p[5]) << 20) | int(p[6])


def is_exact(qname):
    k = key_of(qname.split('__')[0])
    i = np.searchsorted(KEYS, np.uint64(k))
    return bool(i < len(KEYS) and KEYS[i] == k)


def exon_classes(c, blocks):
    """set of exon classes (0 mRNA, 1 ribosomal-protein, 2 rRNA biotype, 3 other biotype) overlapping any aligned block"""
    if c not in EX:
        return set()
    s, e, cl, ml = EX[c]
    out = set()
    for bs, be in blocks:
        i = np.searchsorted(s, be)
        j = np.searchsorted(s, bs - ml)  # candidates: start < be and start >= bs - maxlen
        if i > j:
            sel = e[j:i] > bs
            out.update(cl[j:i][sel].tolist())
    return out


def in_rdna(c, s, e):
    if c not in RD:
        return False
    st, en = RD[c]['s'], RD[c]['e']
    i = np.searchsorted(st, e)
    return bool(i > 0 and (en[:i] > s).any())


def wilson(k, n, z=1.96):
    if not n:
        return (None, None)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (100 * (c - h), 100 * (c + h))


mol = {}  # qname -> [primary(contig, blocks, aseq) or None, set(genomes), rdna_any]
f = pysam.AlignmentFile(bam, 'rb')
for a in f.fetch(until_eof=True):
    if a.is_unmapped or a.is_supplementary:
        continue
    if READS and not a.is_read1:
        continue
    q = a.query_name
    m = mol.get(q)
    if m is None:
        m = mol[q] = [None, set(), False]
    c = a.reference_name
    m[1].add(c.split('_')[0])
    if in_rdna(c, a.reference_start, a.reference_end):
        m[2] = True
    if not a.is_secondary:
        m[0] = (c, a.get_blocks(), a.query_alignment_sequence or '')
LEV = ['L0 mapped, bobcode', 'L1 single species', 'L2 not rRNA/MT', 'L3 mRNA', 'L4 exact bobcode']
n = [0] * 5
co = [0] * 5
cat = collections.Counter()
comp = collections.Counter()
amb = 0
cross = collections.Counter()
wellsp = None
for q, (prim, gen, rdna) in mol.items():
    if prim is None:
        continue
    bc = q.split('__')[1].split('_')[1]
    exp = SPM[bc]
    wellsp = wellsp or exp
    c, blocks, aseq = prim
    sp0 = 'human' if c.startswith('HUMAN_') else 'mouse'
    lvl = 0
    ok0 = sp0 == exp
    sp = 'human' if gen == {'HUMAN'} else 'mouse' if gen == {'MOUSE'} else 'ambiguous'
    ok = sp == exp
    if rdna:
        k = 'rRNA'
    elif c.endswith('_MT'):
        k = 'mitochondrial'
    else:
        cls = exon_classes(c, blocks)
        k = 'rRNA' if 2 in cls else 'ribosomal-protein' if 1 in cls else 'mRNA' if 0 in cls else 'other'
    comp[k] += 1  # composition over ALL mapped molecules (ambiguous included): the rRNA-content number
    if sp != 'ambiguous':
        lvl = 1
        cat[k] += 1
        if k not in ('rRNA', 'mitochondrial'):
            lvl = 2
            if k == 'mRNA':
                lvl = 3
                if is_exact(q):
                    lvl = 4
    else:
        amb += 1
    n[0] += 1
    co[0] += ok0
    for i in range(1, lvl + 1):
        n[i] += 1
        co[i] += ok
    if lvl == 4 and not ok:
        cross[sp] += 1
res = {
    'sample': sample,
    'well_species': wellsp,
    'unit': 'reads' if READS else 'molecules',
    'levels': [],
    'composition_all_mapped': dict(comp),
    'categories_at_L1': dict(cat),
    'n_ambiguous': amb,
    'cross_species_at_final_level_by_aligned_species': dict(cross),
}
for i, lab in enumerate(LEV):
    lo, hi = wilson(co[i], n[i])
    res['levels'].append(
        {
            'level': lab,
            'n': n[i],
            'correct': co[i],
            'accuracy': 100 * co[i] / n[i] if n[i] else None,
            'wilson_lo': lo,
            'wilson_hi': hi,
        }
    )
json.dump(res, open(out, 'w'), indent=1)
print(
    sample,
    wellsp,
    ' | '.join(f"{lab.split()[0]} {100*co[i]/n[i]:.2f}% (n={n[i]:,})" for i, lab in enumerate(LEV) if n[i]),
    flush=True,
)

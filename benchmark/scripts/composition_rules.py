#!/usr/bin/env python3
"""Read composition per benchmark sample with the preprint's classification rules, straight from the per-sample human
BAMs (benchmark_uniform/per_sample_bams; primary mapped records, multimappers included; native BOBseq: mate 1 only = one record per
fragment, as in the pipeline and the read table). Rules, first match wins, on the primary alignment:
  rRNA          span overlaps an rDNA locus (rdna_loci.bed) or an aligned block overlaps an exon of an rRNA / Mt_rRNA / rRNA_pseudogene gene
  mitochondrial contig MT
  ribosomal-protein  aligned block overlaps an exon of a ribosomal-protein gene (RPL/RPS/MRPL/MRPS/RPLP, RPS6K kinases spared)
  mRNA          aligned block overlaps an exon of any other protein-coding gene
  other, split into: exonic other biotype (block overlaps an exon of a non-coding gene), intronic (span inside a gene, no exon overlap), intergenic
"Any alignment" of the rDNA rule cannot be applied here (these BAMs hold primary records only); the pipeline's own composition applies it.
Output: <fig set>/main_figure_panels/values/p3_composition_rules_<tag>_per_sample.tsv (counts) and a per-method summary json.
Usage: BM_COVSET=native|50nt composition_rules.py [--par 8]"""

from settings import CODE, WORK, COVERAGE, REF_DIR, RUN_JSON
import os, re, sys, json, collections, numpy as np, pysam
from concurrent.futures import ProcessPoolExecutor

P = WORK
PB = f'{P}/benchmark_uniform/per_sample_bams'
C = CODE
SET = os.environ.get('BM_COVSET', 'native')
TAG = 'native' if SET == 'native' else 'matched'
FIG = 'preprint_figures_native' if SET == 'native' else f'preprint_figures_{SET}'
OUT = f'{P}/{FIG}/main_figure_panels'
VAL = f'{OUT}/values'
os.makedirs(VAL, exist_ok=True)
IDXF = f'{COVERAGE}/composition_rules_index.npz'
# species switch for the mouse wells of the 24-plex: BM_SPECIES_PREFIX=MOUSE_ uses its own index file; default = human
PREFIX = os.environ.get('BM_SPECIES_PREFIX', 'HUMAN_')
assert PREFIX in ('HUMAN_', 'MOUSE_')
IDXF = IDXF if PREFIX == 'HUMAN_' else IDXF.replace('.npz', '_mouse.npz')
REF = REF_DIR
RIBO = re.compile(r'^(RP[LS]|RPLP|MRP[LS])[^K]*$', re.I)
RRNA_BT = {'rRNA', 'Mt_rRNA', 'rRNA_pseudogene'}
from palette import METH

CLASSES = ['rRNA', 'mitochondrial', 'ribosomal-protein', 'mRNA', 'exonic other biotype', 'intronic', 'intergenic']
if not os.path.exists(
    IDXF
):  # human exons (class 0 mRNA, 1 RP, 2 rRNA biotype, 3 other biotype), human gene spans, rDNA loci; bare contig names
    attr_bt = re.compile(r'gene_biotype "([^"]+)"')
    attr_nm = re.compile(r'gene_name "([^"]+)"')
    ex = collections.defaultdict(set)
    gn = collections.defaultdict(set)
    for l in open(f'{REF}/combined_genome.gtf'):
        if l[0] == '#' or not l.startswith(PREFIX):
            continue
        f = l.split('\t', 9)
        if f[2] not in ('exon', 'gene'):
            continue
        bt = attr_bt.search(f[8])
        nm = attr_nm.search(f[8])
        bt = bt.group(1) if bt else ''
        nm = nm.group(1) if nm else ''
        c = f[0][6:]
        if f[2] == 'gene':
            gn[c].add((int(f[3]) - 1, int(f[4])))
            continue
        cls = (
            2
            if bt in RRNA_BT
            else (1 if (bt == 'protein_coding' and RIBO.match(nm)) else (0 if bt == 'protein_coding' else 3))
        )
        ex[c].add((int(f[3]) - 1, int(f[4]), cls))
    arr = {}
    for c, v in ex.items():
        a = np.array(sorted(v), dtype=np.int64)
        arr[f'ex__{c}__s'] = a[:, 0]
        arr[f'ex__{c}__e'] = a[:, 1]
        arr[f'ex__{c}__c'] = a[:, 2].astype(np.int8)
        arr[f'ex__{c}__ml'] = np.array([int((a[:, 1] - a[:, 0]).max())])
    for c, v in gn.items():
        a = np.array(sorted(v), dtype=np.int64)
        arr[f'gn__{c}__s'] = a[:, 0]
        arr[f'gn__{c}__e'] = a[:, 1]
        arr[f'gn__{c}__ml'] = np.array([int((a[:, 1] - a[:, 0]).max())])
    rd = collections.defaultdict(list)
    for l in open(f'{REF}/rdna_loci.bed'):
        if l[0] == '#' or l.startswith(('track', 'browser')) or not l.strip() or not l.startswith(PREFIX):
            continue
        c, s, e = l.split()[:3]
        rd[c[6:]].append((int(s), int(e)))
    for c, v in rd.items():
        a = np.array(sorted(v), dtype=np.int64)
        arr[f'rd__{c}__s'] = a[:, 0]
        arr[f'rd__{c}__e'] = a[:, 1]
    np.savez_compressed(IDXF, **arr)
    print('index built')
A = np.load(IDXF)
EX, GN, RD = {}, {}, {}
for k in A.files:
    kind, c, part = k.split('__')
    (EX if kind == 'ex' else GN if kind == 'gn' else RD).setdefault(c, {})[part] = A[k]


def ex_classes(c, blocks):
    if c not in EX:
        return set()
    d = EX[c]
    s, e, cl, ml = d['s'], d['e'], d['c'], int(d['ml'][0])
    out = set()
    for bs, be in blocks:
        i = np.searchsorted(s, be)
        j = np.searchsorted(s, bs - ml)
        if i > j:
            out.update(cl[j:i][e[j:i] > bs].tolist())
    return out


def in_gene(c, s0, e0):
    if c not in GN:
        return False
    d = GN[c]
    i = np.searchsorted(d['s'], e0)
    j = np.searchsorted(d['s'], s0 - int(d['ml'][0]))
    return bool(i > j and (d['e'][j:i] > s0).any())


def in_rdna(c, s0, e0):
    if c not in RD:
        return False
    d = RD[c]
    i = np.searchsorted(d['s'], e0)
    return bool(i > 0 and (d['e'][:i] > s0).any())


def classify_bam(args):
    m, n, b, mate1_only = args
    cnt = collections.Counter()
    f = pysam.AlignmentFile(f'{PB}/{b}', 'rb')
    for a in f.fetch(until_eof=True):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or (mate1_only and a.is_paired and not a.is_read1):
            continue
        c, s0, e0 = a.reference_name, a.reference_start, a.reference_end
        if in_rdna(c, s0, e0):
            cnt['rRNA'] += 1
            continue
        if c == 'MT':
            cnt['mitochondrial'] += 1
            continue
        cls = ex_classes(c, a.get_blocks())
        if 2 in cls:
            cnt['rRNA'] += 1
        elif 1 in cls:
            cnt['ribosomal-protein'] += 1
        elif 0 in cls:
            cnt['mRNA'] += 1
        elif 3 in cls:
            cnt['exonic other biotype'] += 1
        elif in_gene(c, s0, e0):
            cnt['intronic'] += 1
        else:
            cnt['intergenic'] += 1
    return m, n, dict(cnt)


if __name__ == '__main__':
    par = int(sys.argv[sys.argv.index('--par') + 1]) if '--par' in sys.argv else 8
    SETNAME = 'native' if SET == 'native' else 'matched_50nt'
    samples = []
    for l in open(f'{PB}/metadata.tsv'):
        f = l.rstrip('\n').split('\t')
        if (
            f[0] in ('drugseq', 'brbseq', 'primeseq')
            and {'drugseq': 'DRUG-seq', 'brbseq': 'BRB-seq', 'primeseq': 'prime-seq'}[f[0]] in METH
            and f[1] == SETNAME
            and f[8] == 'yes'
        ):
            samples.append(
                (
                    {'drugseq': 'DRUG-seq', 'brbseq': 'BRB-seq', 'primeseq': 'prime-seq'}[f[0]],
                    f[2],
                    f[3],
                    False,
                )
            )
    RUN = json.load(open(RUN_JSON))
    ARM = f'{P}/benchmark_uniform/bob57_24plex_pe_native'
    excl = {l.split('\t')[0] for i, l in enumerate(open(f'{ARM}/units_excluded.tsv')) if i > 0}
    bdir = 'bobseq_pe_native' if SET == 'native' else 'bobseq_50nt'
    for w in (l.strip() for l in open(f'{ARM}/wells_keep.txt') if l.strip()):
        if w in excl or RUN['tso_species_map'][w] != 'human':
            continue
        n = RUN['bobcode_labels'][w].replace(' ', '_').replace('/', '_')
        samples.append(('BOBseq', n, f'{bdir}/{n}.bam', SET == 'native'))
    print(TAG, collections.Counter(m for m, *_ in samples), flush=True)
    res = {}
    with ProcessPoolExecutor(par) as ex_:
        for m, n, cnt in ex_.map(classify_bam, samples):
            res[(m, n)] = cnt
            print(' ', m, n, sum(cnt.values()), flush=True)
    with open(f'{VAL}/p3_composition_rules_{TAG}_per_sample.tsv', 'w') as f:
        f.write('method\tsample\tmapped_records\t' + '\t'.join(CLASSES) + '\n')
        for (m, n), cnt in res.items():
            f.write(f'{m}\t{n}\t{sum(cnt.values())}\t' + '\t'.join(str(cnt.get(k, 0)) for k in CLASSES) + '\n')
    summ = {
        m: {
            k: 100
            * sum(res[(mm, n)].get(k, 0) for (mm, n) in res if mm == m)
            / sum(sum(v.values()) for (mm, n), v in res.items() if mm == m)
            for k in CLASSES
        }
        for m in METH
    }
    json.dump(
        {
            'set': SET,
            'rules': 'preprint methods, primary alignment, multimappers included, native BOBseq mate 1 only',
            'per_method_pct_of_mapped': summ,
            'per_sample': {f'{m}|{n}': v for (m, n), v in res.items()},
        },
        open(f'{OUT}/p3_composition_rules_{TAG}.json', 'w'),
        indent=1,
    )
    for m in METH:
        print(m, ' '.join(f'{k}={summ[m][k]:.2f}' for k in CLASSES))
    print('DONE')

#!/usr/bin/env python3
"""Pooled, depth-matched genome-browser tracks per method (the profile matters, not the read count: pool the samples, normalize the depth).
Input: the per-sample human BAMs of the native benchmark set (coverage_architecture/samples.tsv: 24 DRUG-seq, 8 prime-seq, 18 BOBseq; raw pre-dedup alignments).
Unit of depth = a uniquely mapped fragment (MAPQ 255, primary; a BOBseq pair = one fragment, both mates kept). Two variants per method:
  matched   exactly TARGET fragments: a 64-bit hash of sample + read name orders all unique fragments of the method, the TARGET smallest are kept (deterministic,
            uniform, mates together) -> <method>_pooled_matched.bam (+ .bai) and an unscaled bigWig (1-nt bins)
  all       every unique fragment of the method -> bigWig scaled to fragments per 10 M (scale factor = 1e7 / fragments); BAM kept locally only
Contigs: primary chromosomes with UCSC names (chr1..22, chrX, chrY, chrM), scaffolds dropped.   ucsc_pooled_tracks.py <step: hash|write|finish>
"""

from settings import WORK, COVERAGE
import os, sys, json, hashlib, subprocess, numpy as np, pysam
from concurrent.futures import ProcessPoolExecutor

P = WORK
PB = f'{P}/benchmark_uniform/per_sample_bams'
OUT = f'{WORK}/benchmark_uniform/ucsc_pooled'
os.makedirs(f'{OUT}/parts', exist_ok=True)
TARGET = 8_000_000
PRIM = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
UCSC = {c: ('chrM' if c == 'MT' else 'chr' + c) for c in PRIM}
S = [l.rstrip('\n').split('\t') for l in open(f'{COVERAGE}/samples.tsv') if l.strip()]
METH = {'DRUG-seq': 'drugseq', 'prime-seq': 'primeseq', 'BOBseq': 'bobseq'}


def h64(lab, name):
    return int.from_bytes(hashlib.blake2b((lab + '|' + name).encode(), digest_size=8).digest(), 'little')


def keep(a):
    return not (a.flag & 0x904) and a.mapping_quality == 255 and a.reference_name in UCSC


def hash_sample(arg):
    lab, bam, m = arg
    hs = np.fromiter(
        (h64(lab, a.query_name) for a in pysam.AlignmentFile(f'{PB}/{bam}', 'rb').fetch(until_eof=True) if keep(a)),
        dtype=np.uint64,
    )
    n_rec = len(hs)
    hs = np.unique(hs)
    np.save(f'{OUT}/parts/{lab}.hash.npy', hs)
    return lab, m, n_rec, len(hs)


def write_sample(arg):
    lab, bam, m, thr = arg
    src = pysam.AlignmentFile(f'{PB}/{bam}', 'rb')
    hdr = {
        'HD': {'VN': '1.4', 'SO': 'unsorted'},
        'SQ': [{'SN': UCSC[c], 'LN': src.get_reference_length(c)} for c in PRIM],
    }
    rid = {c: i for i, c in enumerate(PRIM)}
    oa = pysam.AlignmentFile(f'{OUT}/parts/{lab}.all.bam', 'wb', header=hdr)
    om = pysam.AlignmentFile(f'{OUT}/parts/{lab}.matched.bam', 'wb', header=hdr)
    na = nm = 0
    for a in src.fetch(until_eof=True):
        if not keep(a):
            continue
        d = a.to_dict()
        d['ref_name'] = UCSC[a.reference_name]
        nr = a.next_reference_name
        d['next_ref_name'] = (
            '=' if (a.is_paired and a.next_reference_id == a.reference_id) else (UCSC[nr] if nr in UCSC else '*')
        )
        if d['next_ref_name'] == '*':
            d['next_ref_pos'] = '0'
            d['length'] = '0'
        d['name'] = lab + ':' + d['name']
        oa.write(pysam.AlignedSegment.from_dict(d, oa.header))
        na += 1
        if h64(lab, a.query_name) <= thr:
            om.write(pysam.AlignedSegment.from_dict(d, om.header))
            nm += 1
    oa.close()
    om.close()
    return lab, m, na, nm


step = sys.argv[1]
if step == 'hash':
    with ProcessPoolExecutor(14) as ex:
        res = list(ex.map(hash_sample, [(r[0], r[1], r[2]) for r in S]))
    info = {}
    for m in METH:
        hs = np.concatenate([np.load(f'{OUT}/parts/{lab}.hash.npy') for lab, mm, _, _ in res if mm == m])
        assert len(np.unique(hs)) == len(hs), f'{m}: 64-bit hash collision between fragments'
        assert len(hs) >= TARGET, f'{m}: only {len(hs):,} unique fragments'
        thr = int(np.partition(hs, TARGET - 1)[TARGET - 1])
        info[m] = {
            'samples': sum(1 for r in res if r[1] == m),
            'unique_records': int(sum(r[2] for r in res if r[1] == m)),
            'unique_fragments': int(len(hs)),
            'hash_threshold': thr,
            'target_fragments': TARGET,
        }
        print(m, info[m], flush=True)
    json.dump(info, open(f'{OUT}/pooling_numbers.json', 'w'), indent=1)
elif step == 'write':
    info = json.load(open(f'{OUT}/pooling_numbers.json'))
    with ProcessPoolExecutor(14) as ex:
        res = list(ex.map(write_sample, [(r[0], r[1], r[2], info[r[2]]['hash_threshold']) for r in S]))
    for m, tag in METH.items():
        for v in ('all', 'matched'):
            parts = [f'{OUT}/parts/{lab}.{v}.bam' for lab, mm, _, _ in res if mm == m]
            cat = f'{OUT}/parts/{tag}.{v}.cat.bam'
            pysam.cat('-o', cat, *parts)
            out = f'{OUT}/{tag}_pooled_{v}.bam'
            pysam.sort('-@', '8', '-m', '2G', '-o', out, cat)
            pysam.index(out)
            os.remove(cat)
            [os.remove(x) for x in parts]
            names = set()
            n = 0
            if v == 'matched':
                for a in pysam.AlignmentFile(out, 'rb').fetch(until_eof=True):
                    names.add(a.query_name)
                    n += 1
                info[m]['matched_fragments_in_bam'] = len(names)
                info[m]['matched_records_in_bam'] = n
                assert len(names) == TARGET, (m, len(names))
            else:
                info[m]['all_records_in_bam'] = int(sum(r[2] for r in res if r[1] == m))
        print(m, 'written', flush=True)
    json.dump(info, open(f'{OUT}/pooling_numbers.json', 'w'), indent=1)
elif step == 'finish':
    info = json.load(open(f'{OUT}/pooling_numbers.json'))
    procs = []
    for m, tag in METH.items():
        for v, sf in (('matched', 1.0), ('all', 1e7 / info[m]['unique_fragments'])):
            name = f'{tag}_pooled_{v}' + ('' if v == 'matched' else '_per10M')
            info[m][f'scale_factor_{v}'] = sf
            cmd = [
                'bamCoverage',
                '-b',
                f'{OUT}/{tag}_pooled_{v}.bam',
                '-o',
                f'{OUT}/{name}.bw',
                '--binSize',
                '1',
                '--scaleFactor',
                f'{sf:.8f}',
                '-p',
                '4',
            ]
            procs.append(
                (
                    name,
                    subprocess.Popen(cmd, stdout=open(f'{OUT}/{name}.bamCoverage.log', 'w'), stderr=subprocess.STDOUT),
                )
            )
    for name, p in procs:
        print(name, 'bigWig exit', p.wait(), flush=True)
    json.dump(info, open(f'{OUT}/pooling_numbers.json', 'w'), indent=1)
print(f'UCSC POOLED {step} DONE')

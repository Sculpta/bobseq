#!/usr/bin/env python3
"""Stage bams: the depth-matched track BAMs -> .cache/bams/<set>_<track>.bam (+ .bai), results/track_depths.tsv.

For each of the twelve wells (control x3, CHX 1 x3, CHX 50 x3, Ris 25 x3): read the per-sample BAM of the benchmark store
(primary alignments, all MAPQ, not deduplicated), keep MAPQ 255 (uniquely mapped), count fragments (read pairs + unpaired
reads). Then, per
track set (settings.TRACK_SETS):
  replicates / ris_replicates  the six wells of control + CHX 50 / control + Ris 25 one by one, each subsampled to the smallest of the six
  pseudobulk / ris_pseudobulk  the conditions, replicate BAMs merged, each subsampled to the smallest condition of the set
Subsampling is exact and by fragment: the read names are sampled without replacement (numpy Generator, SEED) and the
BAM filtered to them (samtools view -N), so a pair is kept or dropped as a unit (the benchmark's "subsampled to exactly N
uniquely mapped fragments"); the smallest member of a set is kept whole. The merged per-condition BAMs are deleted once every
set is built; the MAPQ-255 wells are kept under .cache/bams/q255/ so the matching can be redone without refiltering.
    prepare_bams.py [--force]
"""
import argparse, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from settings import BAMS, RESULTS, RAW_BAMS, SET, WELLS, TRACK_SETS, MAPQ, SEED, THREADS

Q = os.path.join(BAMS, 'q255')
for d in (Q, RESULTS):
    os.makedirs(d, exist_ok=True)


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True, **kw).stdout


def fragments(bam):
    """read pairs (read1 records) + unpaired reads = fragments."""
    return int(sh(f'samtools view -c -f 64 {bam}')) + int(sh(f'samtools view -c -F 1 {bam}'))


def q255(well):
    out = os.path.join(Q, f'{well}.bam')
    if os.path.exists(out + '.done'):
        return out
    raw = os.path.join(RAW_BAMS, f'{well}.bam')
    if not os.path.exists(raw):
        raise SystemExit(f'per-sample BAM not found: {raw} (benchmark/build_per_sample_bams.sh; see locations.py)')
    sh(f'samtools view -b -q {MAPQ} -@ 2 {raw} > {out}.tmp && mv {out}.tmp {out} && samtools index {out}')
    open(out + '.done', 'w').close()
    return out


def subsample(src, n, out, seed):
    """Keep exactly n fragments of src (read names sampled without replacement), sorted + indexed."""
    total = fragments(src)
    if n >= total:                                              # the smallest member of the set: kept whole
        sh(f'cp {src} {out} && samtools index {out}')
        return total
    names = sh(f"samtools view -f 64 {src} | cut -f1") + sh(f"samtools view -F 1 {src} | cut -f1")
    names = np.array(names.split())
    rng = np.random.default_rng(seed)
    keep = rng.choice(names, size=n, replace=False)
    lst = out + '.names'
    open(lst, 'w').write('\n'.join(keep) + '\n')
    sh(f'samtools view -b -N {lst} {src} > {out}.tmp && mv {out}.tmp {out} && samtools index {out}')
    os.remove(lst)
    return len(names)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--force', action='store_true')
    a = ap.parse_args()
    wells = [w for c in WELLS for w in WELLS[c]]
    with ThreadPoolExecutor(3) as ex:
        list(ex.map(q255, wells))
    rows = []
    for setname, tracks in TRACK_SETS.items():
        srcs = {}
        for tid, label, cond, ws in tracks:
            if len(ws) == 1:
                srcs[tid] = os.path.join(Q, f'{ws[0]}.bam')
            else:
                m = os.path.join(Q, f'merged_{cond}.bam')
                if not os.path.exists(m + '.done'):
                    sh(f"samtools merge -f -@ {THREADS} {m}.tmp " + ' '.join(os.path.join(Q, f'{w}.bam') for w in ws) + f" && mv {m}.tmp {m} && samtools index {m}")
                    open(m + '.done', 'w').close()
                srcs[tid] = m
        n = {tid: fragments(srcs[tid]) for tid in srcs}
        target = min(n.values())
        print(f'{setname}: fragments ' + ', '.join(f'{t} {v:,}' for t, v in n.items()) + f' -> matched to {target:,}', flush=True)
        for k, (tid, label, cond, ws) in enumerate(tracks):
            out = os.path.join(BAMS, f'{setname}_{tid}.bam')
            if a.force or not os.path.exists(out + '.done'):
                subsample(srcs[tid], target, out, SEED + k)
                open(out + '.done', 'w').close()
            got = fragments(out)
            rows.append(dict(set=setname, track=tid, label=label, condition=cond, wells=';'.join(ws), fragments_mapq255=n[tid], fragments_matched=got, bam=os.path.relpath(out, os.path.dirname(RESULTS))))
            print(f'  {setname}_{tid}: {got:,} fragments', flush=True)
    for f in os.listdir(Q):                                     # the merged per-condition BAMs (3-4 GB each) are only an intermediate
        if f.startswith('merged_'):
            os.remove(os.path.join(Q, f))
    cols = ['set', 'track', 'label', 'condition', 'wells', 'fragments_mapq255', 'fragments_matched', 'bam']
    with open(os.path.join(RESULTS, 'track_depths.tsv'), 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[c]) for c in cols) + '\n')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Per-well metrics at MATCHED depth, streaming (memory is bounded by wells x max depth).

Every well is subsampled to the same number of usable reads before molecules (definition D)
and genes are counted, so an arm sequenced deeper does not look better for that reason
alone. Rarefaction only goes DOWN: a well with fewer reads than a depth is reported at its
native depth and flagged (at_native = 1).

    rarefy.py <arm_dir> <REG_NH> --depths 50000,100000,250000,500000 [--seed 1]

Reads <arm>/stats_<REG_NH>/wupg.tsv (the same slice every other metric uses) once, keeping a
uniform reservoir of max(depths) rows per well; smaller depths are prefixes of the shuffled
reservoir (a random prefix of a uniform sample is uniform). Writes rarefied.tsv:
    well  depth  reads_native  reads_used  molecules  genes  at_native
Molecules (D) use the same dedup binary and per-arm UMI edit distance as molecule_defs.sh; B_raw/B_kept/B_corr/genes_B
come from umi_b_lib.b_stats, the same estimator as umi_b.py (the common currency).
genes_B_noRP (last column) = genes_B without ribosomal-protein genes (umi_b_lib.RIBO); p2 plots it. Depth, molecules and genes_B are unchanged.
"""

import argparse, os, random, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from umi_b_lib import b_stats

ap = argparse.ArgumentParser()
ap.add_argument('arm')
ap.add_argument('slice')
ap.add_argument('--depths', default='50000,100000,250000,500000,1000000,2000000')
ap.add_argument('--seed', type=int, default=1)
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__))
S = os.path.join(a.arm, f'stats_{a.slice}')
wupg = os.path.join(S, 'wupg.tsv')
if not os.path.isfile(wupg):
    sys.exit(f'rarefy: no {wupg}')
ed = '1'
for ln in open(os.path.join(HERE, 'settings.sh')):
    if ln.startswith('BM_UMI_MAXED='):
        ed = ln.split('=')[1].split()[0]
envf = os.path.join(a.arm, 'arm.env')
if os.path.isfile(envf):
    for ln in open(envf):
        if ln.startswith('UMI_MAXED='):
            ed = ln.split('=')[1].strip()
depths = sorted(int(x) for x in a.depths.split(','))
K = depths[-1]
rng = random.Random(a.seed)
res, seen = {}, {}  # well -> reservoir list, well -> rows seen
for ln in open(wupg):
    w = ln[: ln.index('\t')]
    n = seen.get(w, 0) + 1
    seen[w] = n
    r = res.setdefault(w, [])
    if len(r) < K:
        r.append(ln)
    else:
        j = rng.randrange(n)
        if j < K:
            r[j] = ln
out = open(os.path.join(S, 'rarefied.tsv'), 'w')
out.write(
    'well\tdepth\treads_native\treads_used\tmolecules\tgenes\tat_native\tB_raw\tB_kept\tB_corr\tgenes_B\tgenes_saturated\tgenes_B_noRP\n'
)
tmp = tempfile.mkdtemp(dir=S)
for w in sorted(res):
    r = res[w]
    rng.shuffle(r)
    n = seen[w]
    for d in depths:
        use = r[:d] if n >= d else r
        sub = os.path.join(tmp, 'sub.tsv')
        with open(sub, 'w') as fh:
            fh.writelines(use)
        srt = subprocess.run(
            ['sort', '-k1,1', '-k3,3', '-k6,6', '-k4,4n', sub],
            capture_output=True,
            text=True,
            env={'LC_ALL': 'C', 'PATH': os.environ['PATH']},
        ).stdout
        mol = subprocess.run(
            [os.path.join(HERE, 'dedup'), ed, '0', '/dev/null'], input=srt, capture_output=True, text=True
        ).stdout
        genes = {x.split('\t')[1].strip() for x in mol.splitlines() if '\t' in x}
        bs = b_stats([x.split('\t')[1] for x in use], [x.split('\t')[4] for x in use])
        out.write(
            f'{w}\t{d}\t{n}\t{len(use)}\t{mol.count(chr(10))}\t{len(genes)}\t{int(n < d)}\t{bs["B_raw"]}\t{bs["B_kept"]}\t{bs["B_corr"]}\t{bs["genes"]}\t{bs["saturated"]}\t{bs["genes_noRP"]}\n'
        )
        if n < d:
            break  # smaller-than-depth wells: one native row is enough
out.close()
os.remove(os.path.join(tmp, 'sub.tsv')) if os.path.exists(os.path.join(tmp, 'sub.tsv')) else None
os.rmdir(tmp)
print(f'  {os.path.basename(a.arm)}: {len(res)} wells x {len(depths)} depths -> {S}/rarefied.tsv  (UMI ed<={ed})')

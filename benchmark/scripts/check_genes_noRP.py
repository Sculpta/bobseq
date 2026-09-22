#!/usr/bin/env python3
"""Independent check of rarefied.tsv genes_B_noRP. Rebuilds the rarefy.py reservoir (same seed and call sequence, as genes_thresholds.py),
takes the kept (UMI, gene) pairs per gene from b_stats(per_gene=True) and counts the genes that are not ribosomal-protein genes with a rule written
WITHOUT the regex (prefix test + no K in the remainder). Also checks genes_B is reproduced and lists the ribosomal-protein genes seen in the arm.
    check_genes_noRP.py <arm_dir> [slice]  -> <arm_dir>/stats_<slice>/check_genes_noRP.tsv, exit 1 on any mismatch"""

import os, sys, random

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from umi_b_lib import b_stats

arm = sys.argv[1].rstrip('/')
SL = sys.argv[2] if len(sys.argv) > 2 else 'E_U'
S = f'{arm}/stats_{SL}'


def is_rp(g):
    u = g.upper()
    if u[:4] in ('MRPL', 'MRPS'):
        rest = u[4:]
    elif u[:3] in ('RPL', 'RPS'):
        rest = u[3:]
    else:
        return False
    return 'K' not in rest


depths = [50000, 100000, 250000, 500000, 1000000, 2000000]
K = depths[-1]
rng = random.Random(1)
ref = {}
for ln in open(f'{S}/rarefied.tsv'):
    f = ln.rstrip('\n').split('\t')
    if f[0] != 'well':
        ref[(f[0], int(f[1]))] = (int(f[10]), int(f[12]), f[6])
res, seen = {}, {}
for ln in open(f'{S}/wupg.tsv'):
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
bad = 0
rows = 0
rp_seen = set()
out = open(f'{S}/check_genes_noRP.tsv', 'w')
out.write(
    'well\tdepth\tgenes_B_table\tgenes_B_recount\tgenes_noRP_table\tgenes_noRP_recount\trp_genes_detected\tstatus\n'
)
for w in sorted(res):
    r = res[w]
    rng.shuffle(r)
    n = seen[w]
    for d in depths:
        use = r[:d] if n >= d else r
        gk = b_stats([x.split('\t')[1] for x in use], [x.split('\t')[4] for x in use], per_gene=True).get('gene_k', {})
        rp = {g for g in gk if is_rp(g)}
        rp_seen |= rp
        allg = len(gk)
        norp = allg - len(rp)
        t = ref.get((w, d))
        ok = t is not None and t[0] == allg and t[1] == norp
        bad += not ok
        rows += 1
        out.write(
            f'{w}\t{d}\t{t[0] if t else ""}\t{allg}\t{t[1] if t else ""}\t{norp}\t{len(rp)}\t{"ok" if ok else "MISMATCH"}\n'
        )
        if n < d:
            break
out.close()
open(f'{S}/check_genes_noRP_rp_genes.txt', 'w').write('\n'.join(sorted(rp_seen)) + '\n')
print(
    f'{os.path.basename(arm)} {SL}: {rows} rows checked, {bad} mismatches, {len(rp_seen)} ribosomal-protein genes seen in the arm'
)
sys.exit(1 if bad else 0)

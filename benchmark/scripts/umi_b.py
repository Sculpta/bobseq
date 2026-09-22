#!/usr/bin/env python3
"""Definition B, (well, UMI, gene) = one capture event, with a collision correction that
uses each well's OWN UMI usage instead of the uniform 4^L assumption.

    umi_b.py <arm_dir> <REG> <NH>

Reads  stats_<REG>_<NH>/wupg.tsv   (well  umi  contig  pos  gene  strand), the slice.
Writes stats_<REG>_<NH>/basic/perwell_B.tsv, umi_qc.tsv and prints summary key/values.

Why (measured on a 6-mer UMI library): random hexamers are not used uniformly (effective
tags 1,345 of 4,096 by Simpson), and the over-used tags are homopolymer-like (AAAAAA,
GGGGGG, AGGGGG ...) = failed UMI calls, so the textbook -4^L ln(1 - k/4^L) under-corrects
the top genes by ~20%. Here, per well:
  1. tags flagged as failed calls are excluded: >= L-1 copies of one base (homopolymer
     leakage; 1.9% of uniform hexamers) or used more than BM_UMI_OVERUSE_FOLD times the
     median tag in that well (data-driven);
  2. p_t = that well's usage of the remaining tags; for a gene with k distinct tags the
     molecule count m solves  k = sum_t 1 - (1 - p_t)^m  (bisection on a grid);
  3. a gene whose k sits above the reachable plateau is counted at the plateau and
     reported as saturated (a 6-mer cannot resolve it; D is the fallback there).
For UMIs of >= 9 nt the correction is < 0.1% at any depth we run and m = k is used, so
the 10/14-mer arms are unchanged by construction.
"""

import os, sys, math, collections, subprocess
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
arm, REG, NH = sys.argv[1], (sys.argv[2] if len(sys.argv) > 2 else 'E'), (sys.argv[3] if len(sys.argv) > 3 else 'U')
S = os.path.join(arm, f'stats_{REG}_{NH}')
W = os.path.join(S, 'basic')
os.makedirs(W, exist_ok=True)
OVERUSE_FOLD = float(os.environ.get('BM_UMI_OVERUSE_FOLD', 8))
OVERUSE_MODE = os.environ.get(
    'BM_UMI_OVERUSE_MODE', 'composition'
)  # 'composition' (default) | 'median' (the 6-mer median rule)
CORR_MAX_L = int(os.environ.get('BM_UMI_CORR_MAX_L', 8))  # correct only UMIs this short or shorter
# nominal tag space per UMI length when the UMI has degenerate-but-not-N positions (arm.env:
# UMI_TAG_SPACE="8:36864,11:294912,14:2359296,7:9216" for the 6N+spacer UMI); 4^L otherwise
TAG_SPACE = {}
try:
    for ln in open(os.path.join(arm, 'arm.env')):
        if ln.startswith('UMI_TAG_SPACE='):
            TAG_SPACE = {
                int(k): int(v) for k, v in (x.split(':') for x in ln.split('=', 1)[1].strip().strip('"').split(','))
            }
except FileNotFoundError:
    pass


def homopolymer(u):
    return max(u.count(b) for b in 'ACGT') >= len(u) - 1


# unique (well, umi, gene) rows, grouped by well: sort once, stream once
cmd = f"LC_ALL=C sort -S 8G --parallel=8 -T {W} -u -k1,1 -k2,2 -k5,5 {S}/wupg.tsv | cut -f1,2,5"
proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)

rows_out, qc_out = [], []
tot = collections.Counter()


def flush(well, umis, genes):
    if not umis:
        return
    L = len(umis[0])
    Sfull = TAG_SPACE.get(L, 4**L)  # reported only; the correction uses the measured tag usage
    u = np.array(umis)
    g = np.array(genes)
    tags, tag_idx, tag_n = np.unique(u, return_inverse=True, return_counts=True)
    n_rows = len(u)
    # failed-call flags
    hp = np.array([homopolymer(t) for t in tags])
    if OVERUSE_MODE == 'composition':
        # over-used relative to what the tag's OWN base composition predicts (per-position base frequencies of
        # this well's rows): with a G-rich 6N and 2/3-state spacer positions the legitimate tag frequencies spread
        # over three orders of magnitude, so "8x the median tag" (the 6-mer rule) flags ordinary G-rich tags. A
        # failed call (GGGGGG-like leakage) is still far above its composition expectation.
        M = np.array([list(t) for t in tags])
        freq = np.zeros((L, 4))
        for j in range(L):
            for bi, base in enumerate('ACGT'):
                freq[j, bi] = tag_n[M[:, j] == base].sum() / max(tag_n.sum(), 1)
        expected = np.ones(len(tags))
        for j in range(L):
            for bi, base in enumerate('ACGT'):
                expected[M[:, j] == base] *= max(freq[j, bi], 1e-6)
        expected *= tag_n.sum()
        ov = tag_n > OVERUSE_FOLD * np.maximum(expected, 1) if len(tags) >= 100 else np.zeros(len(tags), bool)
        med = float(np.median(tag_n)) if len(tag_n) else 0
    else:
        med = np.median(tag_n) if len(tag_n) else 0
        ov = tag_n > OVERUSE_FOLD * max(med, 1) if len(tags) >= 100 else np.zeros(len(tags), bool)
    flag = hp | ov
    keep = ~flag[tag_idx]
    B_raw = n_rows
    B_kept = int(keep.sum())
    excluded = n_rows - B_kept
    # usage distribution over kept tags
    kn = tag_n[~flag].astype(float)
    p = kn / kn.sum() if kn.sum() else kn
    simpson = (1.0 / float((p * p).sum())) if len(p) else 0.0
    shannon = float(math.exp(-(p[p > 0] * np.log(p[p > 0])).sum())) if len(p) else 0.0
    # per-gene distinct kept tags
    genes_u, gidx = np.unique(g[keep], return_inverse=True)
    k = np.bincount(gidx, minlength=len(genes_u))
    if L <= CORR_MAX_L and len(p):
        # E[distinct](m) on a log grid, monotone; invert by searchsorted
        mgrid = np.unique(np.round(np.logspace(0, math.log10(max(50 * len(p), 10)), 600)).astype(int))
        E = np.empty(len(mgrid))
        for i, m in enumerate(mgrid):
            E[i] = float((1.0 - np.power(1.0 - p, m)).sum())
        plateau = E[-1]
        sat = k >= plateau * 0.995
        kk = np.minimum(k, plateau * 0.995)
        mhat = np.interp(kk, E, mgrid)  # E is increasing in m
        mhat = np.where(k <= 1, k, mhat)  # a single tag is one molecule
        B_corr = float(mhat.sum())
        n_sat = int(sat.sum())
        top = sorted(zip(k, genes_u, mhat), reverse=True)[:3]
        top_s = '; '.join(f'{gg}:{int(kk_)}->{int(mm)}' for kk_, gg, mm in top)
    else:
        B_corr = float(k.sum())
        n_sat = 0
        top_s = 'no correction (UMI >= %d nt)' % (CORR_MAX_L + 1)
    rows_out.append(
        (
            well,
            n_rows,
            B_kept,
            excluded,
            int(round(B_corr)),
            len(genes_u),
            n_sat,
            len(tags),
            int(hp.sum()),
            int(ov.sum()),
            round(simpson),
            round(shannon),
            Sfull,
            L,
            top_s,
        )
    )
    tot['B_raw'] += B_raw
    tot['B_kept'] += B_kept
    tot['B_corr'] += B_corr
    tot['excluded'] += excluded
    tot['sat_genes'] += n_sat
    tot['wells'] += 1
    tot['hp_tags'] += int(hp.sum())
    tot['ov_tags'] += int(ov.sum())
    qc_out.append(
        (
            well,
            L,
            Sfull,
            len(tags),
            int(hp.sum()),
            int(ov.sum()),
            round(simpson),
            round(shannon),
            round(100.0 * excluded / max(n_rows, 1), 2),
        )
    )


cur, umis, genes = None, [], []
for line in proc.stdout:
    w, u, gname = line.rstrip('\n').split('\t')
    if w != cur:
        flush(cur, umis, genes)
        cur, umis, genes = w, [], []
    umis.append(u)
    genes.append(gname)
flush(cur, umis, genes)
if proc.wait() != 0:
    sys.exit('sort failed')

with open(os.path.join(W, 'perwell_B.tsv'), 'w') as f:
    f.write(
        'well\tB_raw\tB_kept\tumis_excluded\tB_corr\tgenes_B\tgenes_saturated\ttags_seen\ttags_homopolymer\ttags_overused\teff_tags_simpson\teff_tags_shannon\ttag_space\tumi_len\ttop_genes_k_to_m\n'
    )
    for r in rows_out:
        f.write('\t'.join(str(x) for x in r) + '\n')
with open(os.path.join(W, 'umi_qc.tsv'), 'w') as f:
    f.write(
        'well\tumi_len\ttag_space\ttags_seen\ttags_homopolymer\ttags_overused\teff_tags_simpson\teff_tags_shannon\tpct_rows_excluded\n'
    )
    for r in qc_out:
        f.write('\t'.join(str(x) for x in r) + '\n')
med_simpson = int(np.median([r[6] for r in qc_out])) if qc_out else 0
print(f"mol_B_raw\t{tot['B_raw']}")
print(f"mol_B_kept\t{tot['B_kept']}")
print(f"mol_B_corr\t{int(round(tot['B_corr']))}")
print(f"uci\t{int(round(tot['B_corr']))}")  # the common currency by its name
print(f"umi_rows_excluded_pct\t{100.0 * tot['excluded'] / max(tot['B_raw'], 1):.2f}")
print(f"umi_eff_tags_simpson_median\t{med_simpson}")
print(f"genes_saturated_total\t{tot['sat_genes']}")

#!/usr/bin/env python3
"""Per-well definition B with the skew-aware collision correction (see umi_b.py for the why).
b_stats(umis, genes) takes the (umi, gene) pairs of ONE well (duplicates allowed) and returns
dict(B_raw, B_kept, excluded, B_corr, genes, saturated, tags, hp, ov, simpson, shannon, L)."""

import math, os, re
import numpy as np

# ribosomal-protein genes (RPL/RPS/RPLP/MRPL/MRPS; RPS6K kinases spared): the same rule as composition_rules.py, picard_profile.py and
# threshold_panels.py. genes_noRP = detected genes without them (ribosomal-protein genes are always removed from the gene set)
RIBO = re.compile(r'^(RP[LS]|RPLP|MRP[LS])[^K]*$', re.I)
OVERUSE_FOLD = float(os.environ.get('BM_UMI_OVERUSE_FOLD', 8))
OVERUSE_MODE = os.environ.get('BM_UMI_OVERUSE_MODE', 'composition')  # must match umi_b.py
CORR_MAX_L = int(os.environ.get('BM_UMI_CORR_MAX_L', 8))


def homopolymer(u):
    return max(u.count(b) for b in 'ACGT') >= len(u) - 1


def b_stats(umis, genes, per_gene=False):
    """per_gene=True adds gene_k: kept (UMI, gene) pairs per gene, the count the collision correction starts from (used by the threshold variants of p2; rule unchanged)"""
    pairs = np.unique(np.array([u + '\t' + g for u, g in zip(umis, genes)]))
    if not len(pairs):
        return dict(
            B_raw=0,
            B_kept=0,
            excluded=0,
            B_corr=0,
            genes=0,
            genes_noRP=0,
            saturated=0,
            tags=0,
            hp=0,
            ov=0,
            simpson=0,
            shannon=0,
            L=0,
        )
    u = np.array([x.split('\t')[0] for x in pairs])
    g = np.array([x.split('\t')[1] for x in pairs])
    L = len(u[0])
    S = 4**L
    tags, tag_idx, tag_n = np.unique(u, return_inverse=True, return_counts=True)
    hp = np.array([homopolymer(t) for t in tags])
    if OVERUSE_MODE == 'composition':
        # same rule as umi_b.py flush(): over-used relative to the tag's own base-composition expectation
        # (the two must stay identical so that rarefied.tsv agrees with perwell_B.tsv)
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
    else:
        med = np.median(tag_n) if len(tag_n) else 0
        ov = tag_n > OVERUSE_FOLD * max(med, 1) if len(tags) >= 100 else np.zeros(len(tags), bool)
    flag = hp | ov
    keep = ~flag[tag_idx]
    B_raw = len(u)
    B_kept = int(keep.sum())
    kn = tag_n[~flag].astype(float)
    p = kn / kn.sum() if kn.sum() else kn
    simpson = (1.0 / float((p * p).sum())) if len(p) else 0.0
    shannon = float(math.exp(-(p[p > 0] * np.log(p[p > 0])).sum())) if len(p) else 0.0
    genes_u, gidx = np.unique(g[keep], return_inverse=True)
    k = np.bincount(gidx, minlength=len(genes_u))
    n_sat = 0
    if L <= CORR_MAX_L and len(p):
        mgrid = np.unique(np.round(np.logspace(0, math.log10(max(50 * len(p), 10)), 600)).astype(int))
        E = np.array([float((1.0 - np.power(1.0 - p, m)).sum()) for m in mgrid])
        plateau = E[-1]
        sat = k >= plateau * 0.995
        mhat = np.interp(np.minimum(k, plateau * 0.995), E, mgrid)
        mhat = np.where(k <= 1, k, mhat)
        B_corr = float(mhat.sum())
        n_sat = int(sat.sum())
    else:
        B_corr = float(k.sum())
    out = dict(
        B_raw=B_raw,
        B_kept=B_kept,
        excluded=B_raw - B_kept,
        B_corr=int(round(B_corr)),
        genes=len(genes_u),
        genes_noRP=sum(1 for x in genes_u.tolist() if not RIBO.match(x)),
        saturated=n_sat,
        tags=len(tags),
        hp=int(hp.sum()),
        ov=int(ov.sum()),
        simpson=round(simpson),
        shannon=round(shannon),
        L=L,
    )
    if per_gene:
        out['gene_k'] = dict(zip(genes_u.tolist(), k.tolist()))
    return out

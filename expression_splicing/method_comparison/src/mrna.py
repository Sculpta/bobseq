"""Stage mrna (the quant.sf route, after stage `counts`): restrict the gene tables to the mRNA gene set, fold the haplotype
copies, depth-match the counts and complete the gene lengths.
  1. the mRNA gene set of the preprint (S10): every primary-assembly Ensembl 113 protein_coding gene except the 13 on
     chromosome MT and the 169 ribosomal-protein genes (symbol ^(RPL|RPS|MRPL|MRPS), RPS6K* kept), 19,934 genes;
  2. alternative-haplotype / patch copies (Salmon splits reads between a gene and its copies) are added to their primary
     gene when the copy's symbol matches exactly one primary protein-coding gene: counts and TPM summed, the effective length
     combined count-weighted harmonically; copies without a unique match and non-coding genes are dropped;
  3. counts rounded to integers, every sample thinned without replacement (seed 0) to the smallest mRNA library;
  4. gene effective length per sample, complete: the TPM-weighted value where the gene has reads, else the mean over the same
     method's samples with reads, else the sample's unweighted mean transcript length, else 1 (genes without a Salmon row).
Writes results/gene_counts_mrna.csv.gz, gene_counts_mrna_thinned.csv.gz, gene_tpm_mrna.csv.gz, gene_efflen_mrna.csv.gz,
mrna_summary.tsv (per sample: counts kept and dropped by class, thinned depth, genes detected, length fallbacks) and
haplotype_collapse.tsv (copy -> primary gene, with the counts moved)."""
import collections, os, re, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, GTF, MRNA, MRNA_THIN, EFFLEN, table, open_text, downsample, write_tsv
SEED = 0
GENE_SET = table("S10")                                                   # the one mRNA gene set of every expression analysis

def genes_from_gtf():
    """gene_id -> (gene_name, gene_biotype) for the primary assembly (the Ensembl 113 GTF); chromosome-MT protein-coding genes are
    relabelled mt_protein_coding so that only NUCLEAR protein-coding genes pass the protein_coding filter."""
    out = {}
    for line in open_text(GTF):
        if line[0] == "#": continue
        f = line.split("\t", 9)
        if f[2] != "gene": continue
        gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1); name = re.search(r'gene_name "([^"]+)"', f[8]); bt = re.search(r'gene_biotype "([^"]+)"', f[8]).group(1)
        out[gid] = (name.group(1) if name else "", "mt_protein_coding" if (bt == "protein_coding" and f[0] == "MT") else bt)   # the 13 mitochondrially encoded mRNAs are not mRNA for this analysis
    return out

def load(name):
    return pd.read_csv(os.path.join(RESULTS, f"gene_{name}.csv.gz"), index_col=0).drop(columns=["SYMBOL"])

def main():
    ann = genes_from_gtf(); M = pd.read_csv(os.path.join(RESULTS, "gene_counts.csv.gz"), index_col=0)
    sym = M["SYMBOL"].fillna("").str.replace(r" \[ENSG\d+\]$", "", regex=True); C = M.drop(columns=["SYMBOL"]); samples = list(C.columns)
    A, W, U = load("tpm"), load("efflen"), load("efflen_unweighted")                   # Salmon gene TPM, weighted / unweighted effective length
    method = dict(zip(*(lambda t: (t["sample"], t["method"]))(pd.read_csv(os.path.join(RESULTS, "sample_table.tsv"), sep="\t"))))
    # every primary-assembly protein-coding gene of the GTF, not only those in the Salmon matrix: 540 are absent from quant.sf
    # because all their transcripts are exact duplicates of a retained transcript (Salmon's duplicate_clusters.tsv), which in
    # 2,792 of 2,957 cases is the haplotype / patch copy -- their reads sit under the copy's gene id and come back here
    gene_set = pd.read_csv(GENE_SET, sep="\t"); pc = sorted(gene_set.gene_id); assert all(g in ann for g in pc), "gene set / GTF mismatch"
    off = [g for g in C.index if g not in ann]                                                   # alt-haplotype / patch gene ids (not in the primary-assembly GTF)
    by_sym = collections.defaultdict(list)
    for g in pc:
        if ann[g][0]: by_sym[ann[g][0]].append(g)
    target = {g: by_sym[sym[g]][0] for g in off if sym[g] in by_sym and len(by_sym[sym[g]]) == 1}
    P = C.reindex(pc, fill_value=0.0)
    moved = C.loc[list(target)].groupby(pd.Series(target)).sum(); P.loc[moved.index] += moved
    # --- the same collapse for TPM (summed) and the effective lengths (count-weighted harmonic; unweighted: primary, else copies)
    tgt = pd.Series(target); PA = A.reindex(pc, fill_value=0.0); PA.loc[moved.index] += A.loc[list(target)].groupby(tgt).sum()
    inv = (C / W).where(C > 0)                                                       # counts / L per gene, NaN without reads
    PI = inv.reindex(pc); PI.loc[moved.index] = PI.loc[moved.index].fillna(0.0) + inv.loc[list(target)].groupby(tgt).sum(min_count=1).fillna(0.0)
    PW = (P / PI).where(P > 0)                                                        # L = sum(counts) / sum(counts / L_i); NaN without reads
    PU = U.reindex(pc); cu = U.loc[list(target)].groupby(tgt).mean(); PU.loc[cu.index] = PU.loc[cu.index].fillna(cu)
    coll = [{"copy_gene_id": g, "symbol": sym[g], "primary_gene_id": t, **{s: round(float(C.at[g, s]), 3) for s in samples}} for g, t in target.items()]
    write_tsv(os.path.join(RESULTS, "haplotype_collapse.tsv"), coll, ["copy_gene_id", "symbol", "primary_gene_id"] + samples)
    out = P.copy(); out.insert(0, "SYMBOL", [ann[g][0] for g in out.index]); out.index.name = "GENEID"; out.to_csv(MRNA)
    R = P.round().astype(int); depth = int(R.sum().min()); rng = np.random.default_rng(SEED)
    T = pd.DataFrame({s: downsample(R[s].to_numpy(), depth, rng) for s in samples}, index=R.index)
    outT = T.copy(); outT.insert(0, "SYMBOL", [ann[g][0] for g in outT.index]); outT.index.name = "GENEID"; outT.to_csv(MRNA_THIN)
    outA = PA.copy(); outA.insert(0, "SYMBOL", [ann[g][0] for g in outA.index]); outA.index.name = "GENEID"; outA.to_csv(os.path.join(RESULTS, "gene_tpm_mrna.csv.gz"))
    # --- complete the effective-length matrix: weighted where the gene has reads -> same-method mean -> unweighted -> 1 (never used)
    E = PW.copy(); fb_method, fb_unw, fb_one = {}, {}, {}
    for s_ in samples:
        peers = [x for x in samples if method[x] == method[s_]]; m1 = PW[peers].mean(1)
        need = E[s_].isna(); E.loc[need, s_] = m1[need]; fb_method[s_] = int((need & m1.notna()).sum())
        need = E[s_].isna(); E.loc[need, s_] = PU.loc[need, s_]; fb_unw[s_] = int((need & PU[s_].notna()).sum())
        need = E[s_].isna(); fb_one[s_] = int(need.sum()); E.loc[need, s_] = 1.0
    unused = [g for g in E.index if (E.loc[g] == 1.0).any()]; assert (P.loc[unused].sum(1) == 0).all(), "a gene with reads got the placeholder length"
    outE = E.copy(); outE.insert(0, "SYMBOL", [ann[g][0] for g in outE.index]); outE.index.name = "GENEID"; outE.to_csv(EFFLEN)
    dropped_off = [g for g in off if g not in target]; inset = set(pc)
    nonset = [g for g in C.index if g in ann and g not in inset]                                  # primary-assembly genes outside the set: non-coding, MT, ribosomal protein
    is_rp = lambda n: bool(re.match(r"^(RPL|RPS|MRPL|MRPS)", n)) and not n.startswith("RPS6K")
    rp = [g for g in nonset if ann[g][1] in ("protein_coding", "mt_protein_coding") and is_rp(ann[g][0])]; mtg = [g for g in nonset if ann[g][1] == "mt_protein_coding"]
    nonpc = [g for g in nonset if g not in rp and g not in mtg]
    rows = []
    for s in samples:
        tot = float(C[s].sum())
        rows.append({"sample": s, "salmon_gene_counts": round(tot), "mrna_counts": round(float(P[s].sum())), "mrna_pct": round(100 * P[s].sum() / tot, 2),
                     "collapsed_copy_counts": round(float(moved[s].sum())), "collapsed_copy_pct": round(100 * moved[s].sum() / tot, 2),
                     "dropped_non_protein_coding_pct": round(100 * C.loc[nonpc, s].sum() / tot, 2), "dropped_mt_pct": round(100 * C.loc[mtg, s].sum() / tot, 2), "dropped_ribosomal_protein_pct": round(100 * C.loc[rp, s].sum() / tot, 2), "dropped_off_assembly_pct": round(100 * C.loc[dropped_off, s].sum() / tot, 3),
                     "thinned_to": depth, "genes_ge1_thinned": int((T[s] >= 1).sum()), "genes_ge5_thinned": int((T[s] >= 5).sum()), "genes_ge10_thinned": int((T[s] >= 10).sum()),
                     "efflen_from_reads": int(PW[s].notna().sum()), "efflen_from_method_mean": fb_method[s], "efflen_unweighted": fb_unw[s], "efflen_placeholder": fb_one[s]})
    write_tsv(os.path.join(RESULTS, "mrna_summary.tsv"), rows, list(rows[0].keys()))
    print(f"mRNA gene set {len(pc):,} ({sum(g in C.index for g in pc):,} in the Salmon matrix, {sum(g not in C.index for g in pc):,} only via a collapsed copy or zero); haplotype/patch copies collapsed {len(target):,} onto {moved.shape[0]:,} genes; "
          f"dropped: {len(nonpc):,} non-protein-coding + {len(mtg)} chrMT + {len(rp)} ribosomal-protein + {len(dropped_off):,} off-assembly without a unique primary match; thinned every sample to {depth:,} mRNA reads")
    print(pd.DataFrame(rows)[["sample", "salmon_gene_counts", "mrna_counts", "mrna_pct", "collapsed_copy_pct", "dropped_non_protein_coding_pct", "dropped_mt_pct", "dropped_ribosomal_protein_pct", "genes_ge5_thinned", "efflen_from_reads", "efflen_from_method_mean", "efflen_unweighted", "efflen_placeholder"]].to_string(index=False))

if __name__ == "__main__":
    main()

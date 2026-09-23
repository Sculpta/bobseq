"""Stage expression: the expression matrices the PCA and the scatters read. Per sample s and gene g, from stage `mrna` or `tables`:
  BOBseq, TruSeq (TPM_METHODS)   a = count / effective length     (reads per nucleotide of effective length)
  prime-seq, BRB-seq             a = count                        (a 3'-tag count carries no length term)
  expression = a / sum over genes x 1e6,  on the native counts and on the thinned counts (thin first, divide second).
On the native counts the BOBseq / TruSeq values equal Salmon's gene TPM renormalised to the gene set (checked here).
Writes results/expression_mrna.csv.gz, expression_mrna_thinned.csv.gz (GENEID, SYMBOL, 29 samples) and expression_units.tsv
(per sample: unit, totals, count-weighted mean effective length of the detected genes, the TPM check, what the division did to
genes < 1 kb and >= 3 kb relative to CPM)."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, MRNA, MRNA_THIN, EFFLEN, EXPR, EXPR_THIN, UNIT, TPM_METHODS, write_tsv
MIN_COUNT = 5

def convert(C, L, method):
    """counts -> per-million expression: counts / effective length for the TPM methods, counts as they are for the others."""
    a = pd.DataFrame({s: (C[s] / L[s] if method[s] in TPM_METHODS else C[s]) for s in C.columns})
    return a / a.sum() * 1e6, a.sum()

def main():
    M = pd.read_csv(MRNA, index_col=0); sym = M["SYMBOL"].fillna(""); C = M.drop(columns=["SYMBOL"]).astype(float)
    T = pd.read_csv(MRNA_THIN, index_col=0).drop(columns=["SYMBOL"]).astype(float).loc[C.index]; L = pd.read_csv(EFFLEN, index_col=0).drop(columns=["SYMBOL"]).loc[C.index]
    A = pd.read_csv(os.path.join(RESULTS, "gene_tpm_mrna.csv.gz"), index_col=0).drop(columns=["SYMBOL"]).loc[C.index]
    samples = pd.read_csv(os.path.join(RESULTS, "sample_table.tsv"), sep="\t"); method = dict(zip(samples["sample"], samples["method"]))
    E, tot = convert(C, L, method); ET, totT = convert(T, L, method)
    for X, path in ((E, EXPR), (ET, EXPR_THIN)):
        out = X.copy(); out.insert(0, "SYMBOL", sym); out.index.name = "GENEID"; out.to_csv(path)
    rows = []
    for s in C.columns:
        det = T[s] >= MIN_COUNT; row = {"sample": s, "method": method[s], "unit": UNIT[method[s]], "native_counts": int(round(C[s].sum())), "thinned_counts": int(T[s].sum()),
                                        "sum_a_native": round(float(tot[s]), 3), "sum_a_thinned": round(float(totT[s]), 3), "genes_detected_thinned": int(det.sum()),
                                        "efflen_count_weighted_mean_detected": round(float((T.loc[det, s] * L.loc[det, s]).sum() / T.loc[det, s].sum()))}
        if method[s] in TPM_METHODS:                                    # native TPM must equal Salmon's gene TPM renormalised to this gene set
            ref = A[s] / A[s].sum() * 1e6; ok = C[s] >= 10; row["max_rel_dev_vs_salmon_tpm"] = f"{float((E.loc[ok, s] / ref[ok] - 1).abs().max()):.1e}"
            hi, lo = L[s] >= 3000, L[s] < 1000            # what the division does: long genes lose, short genes gain, relative to CPM
            cpm = T[s] / T[s].sum() * 1e6; d = np.log2((ET[s] + 1) / (cpm + 1))
            row["median_log2_tpm_over_cpm_lt1kb"] = round(float(d[lo & det].median()), 2); row["median_log2_tpm_over_cpm_ge3kb"] = round(float(d[hi & det].median()), 2)
        rows.append(row)
    cols = list(rows[0].keys()) + [k for k in rows[-1].keys() if k not in rows[0]]
    cols = ["sample", "method", "unit", "native_counts", "thinned_counts", "sum_a_native", "sum_a_thinned", "genes_detected_thinned", "efflen_count_weighted_mean_detected", "max_rel_dev_vs_salmon_tpm", "median_log2_tpm_over_cpm_lt1kb", "median_log2_tpm_over_cpm_ge3kb"]
    write_tsv(os.path.join(RESULTS, "expression_units.tsv"), rows, cols)
    print(pd.DataFrame(rows)[cols].fillna("").to_string(index=False))
    print(f"\nexpression_mrna[_thinned].csv.gz: {E.shape[0]:,} genes x {E.shape[1]} samples; TPM for {', '.join(TPM_METHODS)}, CPM for the rest; every column sums to 1e6")

if __name__ == "__main__":
    main()

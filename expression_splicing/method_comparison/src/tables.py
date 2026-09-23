"""Stage tables: the inputs of every later stage from the supplementary data tables, without Salmon quant.sf files.

S7 (tximport gene counts), S8 (gene TPM) and S9 (tximport gene effective length) hold the 129 samples of the preprint (119
sample-arms + 10 public TruSeq runs) on the 19,934-gene mRNA set with the haplotype copies folded: the state stages counts + mrna
reach from quant.sf. This stage selects the 29 samples of counts.SAMPLES, renames them, thins every sample to the smallest mRNA
library (common.downsample, seed 0, samples in list order) and writes results/gene_counts_mrna.csv.gz,
gene_counts_mrna_thinned.csv.gz, gene_tpm_mrna.csv.gz, gene_efflen_mrna.csv.gz (the 67 genes without a Salmon row get the
placeholder length 1; they have no reads anywhere) and results/sample_table.tsv when absent.
The deposited tables were summarised with tximport in R, stage counts with pandas; one gene (GOLGA8M) differs by <= 1.9 counts, so
the thinned draw and the numbers downstream can differ from the quant.sf route by that rounding."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, MRNA, MRNA_THIN, EFFLEN, UNIT, table, downsample, write_tsv
from counts import SAMPLES
SEED = 0

def load(prefix):
    M = pd.read_csv(table(prefix), index_col=0); return M.drop(columns=["SYMBOL"]).astype(float), M["SYMBOL"].fillna("")

def write(df, sym, path):
    out = df.copy(); out.insert(0, "SYMBOL", sym.reindex(out.index).fillna("").to_numpy()); out.index.name = "GENEID"; out.to_csv(path)

def main():
    os.makedirs(RESULTS, exist_ok=True)
    C, sym = load("S7"); A, _ = load("S8"); L, _ = load("S9")
    names = [s[1] for s in SAMPLES]; src = {s[1]: (s[3] if s[2] == "local" else "truseq_public__" + s[3]) for s in SAMPLES}
    missing = [src[n] for n in names if src[n] not in C.columns]
    if missing: raise SystemExit(f"columns missing from S7: {missing}")
    genes = sorted(C.index)                                                  # the gene order of stage mrna (sorted gene ids)
    P = C.loc[genes, [src[n] for n in names]]; P.columns = names
    PA = A.loc[genes, [src[n] for n in names]]; PA.columns = names
    E = L.loc[genes, [src[n] for n in names]]; E.columns = names; E = E.fillna(1.0)
    write(P, sym, MRNA); write(PA, sym, os.path.join(RESULTS, "gene_tpm_mrna.csv.gz")); write(E, sym, EFFLEN)
    R = P.round().astype(int); depth = int(R.sum().min()); rng = np.random.default_rng(SEED)
    T = pd.DataFrame({s: downsample(R[s].to_numpy(), depth, rng) for s in names}, index=R.index); write(T, sym, MRNA_THIN)
    st = os.path.join(RESULTS, "sample_table.tsv")
    if not os.path.exists(st):
        rows = [{"sample": n, "method": m, "group": g, "expression_unit": UNIT[m], "source": src[n], "unit": "dedup reads" if k == "local" else "reads (no UMI)"} for m, n, k, _, g in SAMPLES]
        write_tsv(st, rows, list(rows[0].keys()))
    print(f"{len(genes):,} genes x {len(names)} samples from S7 / S8 / S9; every sample thinned to {depth:,} mRNA reads; wrote gene_counts_mrna, gene_counts_mrna_thinned, gene_tpm_mrna, gene_efflen_mrna")

if __name__ == "__main__":
    main()

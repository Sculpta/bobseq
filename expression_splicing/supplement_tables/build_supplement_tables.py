#!/usr/bin/env python3
"""How the supplementary data tables of the preprint were assembled (provenance; the tables are deposited with the preprint,
this script ran on the internal collection directories and is not runnable from the repository).

The smallest set of tables from which every expression and splicing figure can be recreated. Dense, gzip-compressed CSV.

  S1  sample_metadata.csv                   the 119 sample-arms (= expression_splicing/sample_metadata.csv)
  S2  public_truseq_runs.tsv                the 10 public TruSeq HEK293T runs of the method comparison (GEO / SRA ids, kit, layout, spots)
  S3  psi_dedup.csv.gz                      PSI, UMI-deduplicated molecules: 8,701 junctions x 119 arms (empty = junction not detected in the sample)
  S4  junction_counts_dedup.csv.gz          the molecule counts behind S3 (same rows / columns)
  S5  psi_fragments.csv.gz                  PSI, per fragment (no deduplication): 36,433 junctions x 119 arms -- the unit of screen_splicing
  S6  junction_counts_fragments.csv.gz      the fragment counts behind S5
  S7  gene_counts_salmon.csv.gz             Salmon / tximport gene counts, the 19,934-gene mRNA set x (119 arms + 10 TruSeq runs)
  S8  gene_tpm_salmon.csv.gz                Salmon gene TPM re-scaled over the 19,934 genes, same columns
  S9  gene_effective_length_salmon.csv.gz   tximport gene effective length (TPM-weighted mean transcript effective length) per sample, same columns
  S10 mrna_gene_set.tsv                     the 19,934 genes (Ensembl 113 protein-coding minus chrMT minus ribosomal-protein genes) with symbol / chromosome
  S11 gene_counts_featurecounts.csv.gz      featureCounts per Ensembl gene, all biotypes, 119 arms (the matrix of replicate_correlation)

Inputs, as they were on the analysis machine: the collection directory (metadata, counts, annotations, the PSI matrices in
MatrixMarket form) of the UMI-deduplicated build, its per-fragment counterpart, and the results of the 35-sample method
comparison for the ten public TruSeq runs.
"""
import json, os, shutil
import numpy as np, pandas as pd, scipy.io
HERE = os.path.dirname(os.path.abspath(__file__))
COLL = os.environ.get("COLLECTION_DIR", "")            # the UMI-deduplicated collection: metadata/, counts/, annotations/, psi/
READS = os.environ.get("READS_COLLECTION_DIR", "")     # the per-fragment collection: counts/, psi/
FULL = os.environ.get("METHOD_COMPARISON_FULL", "")    # results of the 35-sample method comparison (gene_counts_mrna.csv.gz, sample_table.tsv, truseq_samples.tsv)
FULL_TPM = os.environ.get("METHOD_COMPARISON_FULL_TPM", "")   # its TPM variant (gene_tpm_mrna.csv.gz, gene_efflen_mrna.csv.gz)

def dense(d, sub, name):
    """MatrixMarket sparse matrix + rownames + colnames -> dense DataFrame; explicit entries only (missing = not detected)."""
    M = scipy.io.mmread(os.path.join(d, sub, f"{name}.mtx")).tocsr()
    rn = open(os.path.join(d, sub, f"{name}.rownames.txt")).read().split("\n")[:-1]; cn = open(os.path.join(d, sub, f"{name}.colnames.txt")).read().split("\n")[:-1]
    assert M.shape == (len(rn), len(cn)), (M.shape, len(rn), len(cn))
    A = M.toarray().astype(float); A[(M.toarray() == 0)] = np.nan          # a stored zero never occurs (entries only for count >= 1)
    return pd.DataFrame(A, index=pd.Index(rn, name="junction"), columns=cn), M.nnz

def main():
    log = []
    # S1 / S2 / S10 -- copies
    shutil.copy(os.path.join(COLL, "metadata", "sample_metadata.csv"), os.path.join(HERE, "S1_sample_metadata.csv"))
    t = pd.read_csv(os.path.join(FULL, "truseq_samples.tsv"), sep="\t"); st = pd.read_csv(os.path.join(FULL, "sample_table.tsv"), sep="\t")
    t.insert(0, "column_in_gene_tables", "truseq_public__" + t["run"]); t.to_csv(os.path.join(HERE, "S2_public_truseq_runs.tsv"), sep="\t", index=False)
    gs = pd.read_csv(os.path.join(COLL, "annotations", "mrna_gene_set.tsv"), sep="\t"); gs.to_csv(os.path.join(HERE, "S10_mrna_gene_set.tsv"), sep="\t", index=False)
    # S3-S6 -- PSI and junction counts, both builds
    for build, d, tag in (("dedup", COLL, "dedup"), ("fragments", READS, "fragments")):
        psi, nnz = dense(d, "psi", "psi_matrix"); cnt, nnz2 = dense(d, "counts", "junction_counts_matrix")
        assert list(psi.index) == list(cnt.index) and list(psi.columns) == list(cnt.columns) and nnz == nnz2
        psi.to_csv(os.path.join(HERE, f"{'S3' if build == 'dedup' else 'S5'}_psi_{tag}.csv.gz"), float_format="%.3f")
        cnt.to_csv(os.path.join(HERE, f"{'S4' if build == 'dedup' else 'S6'}_junction_counts_{tag}.csv.gz"), float_format="%.0f")
        log.append(f"{build}: {psi.shape[0]:,} junctions x {psi.shape[1]} arms, {nnz:,} detected entries")
    # S7-S9 -- gene tables: the collection's 119 arms + the 10 TruSeq runs (from the 35-sample inter-sample results, same tximport rule, same gene set)
    tru = [s_ for s_ in st["sample"] if s_.startswith("TruSeq")]; ren = {s_: "truseq_public__" + src for s_, src in zip(st["sample"], st["source"]) if s_ in tru}
    for tag, coll_file, isv_dir, isv_file in (("S7_gene_counts", "salmon_gene_counts.csv.gz", FULL, "gene_counts_mrna.csv.gz"),
                                              ("S8_gene_tpm", "salmon_gene_tpm.csv.gz", FULL_TPM, "gene_tpm_mrna.csv.gz"),
                                              ("S9_gene_effective_length", "salmon_gene_length.csv.gz", FULL_TPM, "gene_efflen_mrna.csv.gz")):
        C = pd.read_csv(os.path.join(COLL, "counts", coll_file)); assert list(C.columns[:2]) == ["GENEID", "SYMBOL"] and len(C) == len(gs)
        I = pd.read_csv(os.path.join(isv_dir, isv_file), index_col=0); T = I[tru].rename(columns=ren)
        if tag == "S8_gene_tpm": T = T / T.sum() * 1e6                      # Salmon TPM re-scaled over the 19,934 genes, as the collection's table
        out = C.set_index("GENEID").join(T, how="left"); assert len(out) == len(gs) and out[list(ren.values())].isna().sum().sum() == 0
        out.to_csv(os.path.join(HERE, f"{tag}_salmon.csv.gz"), float_format="%.4f"); log.append(f"{tag}: {out.shape[0]:,} genes x {out.shape[1] - 1} samples")
    # S11 -- featureCounts, all biotypes
    F = pd.read_csv(os.path.join(COLL, "counts", "gene_counts.csv.gz")); F.to_csv(os.path.join(HERE, "S11_gene_counts_featurecounts.csv.gz"), index=False); log.append(f"S11: {F.shape[0]:,} genes x {F.shape[1] - 2} arms")
    print("\n".join(log))

if __name__ == "__main__":
    main()

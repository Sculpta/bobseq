#!/usr/bin/env python3
"""Inter-sample variability of mRNA expression per method: BOBseq's three untreated controls against every control sample of
prime-seq, BRB-seq and TruSeq, length-normalised (TPM) where the library is full-length.   run.py [--only STAGE ...] [--list]

Without --only the stages tables, expression, pca, de, capped_views and scatter run (the route from the supplementary data
tables). `counts` and `mrna` are the route from Salmon quant.sf files and replace `tables`."""
import argparse, os, subprocess, sys
P = os.path.dirname(os.path.abspath(__file__))
STAGES = [
    ("tables",       "the inputs from the supplementary data tables S7 / S8 / S9 (no quant.sf needed): the 29 samples, thinned to the smallest mRNA library -> results/gene_counts_mrna[_thinned].csv.gz, gene_tpm_mrna.csv.gz, gene_efflen_mrna.csv.gz", ["src/tables.py"]),
    ("counts",       "the quant.sf route: gene counts (tximport rule), Salmon gene TPM and gene effective lengths from <SALMON_QUANT>/<source>/quant.sf -> results/gene_{counts,tpm,efflen,efflen_unweighted}.csv.gz", ["src/counts.py"]),
    ("mrna",         "the quant.sf route: the mRNA gene set (S10), haplotype copies collapsed, every sample thinned to the smallest mRNA library -> results/gene_counts_mrna[_thinned].csv.gz, gene_efflen_mrna.csv.gz", ["src/mrna.py"]),
    ("expression",   "TPM for BOBseq and TruSeq (counts / effective length), CPM for prime-seq and BRB-seq, per million over the gene set, native and thinned -> results/expression_mrna[_thinned].csv.gz", ["src/expression.py"]),
    ("pca",          "PCA + clustered Spearman heatmap + dendrogram of the 29 samples, coloured by method (native and thinned runs); Fig. 6A = figures/pca_thinned_labs.svg (TruSeq marked by laboratory)", ["src/pca.py"]),
    ("de",           "DESeq2 (R) between methods on the thinned counts, one gene length per method as normalisation factors for BOBseq / TruSeq, 6 pairwise contrasts, twice: the full mRNA gene set -> results/de; the capped common set (7,500 genes) -> results/de_capped", ["src/de.py"]),
    ("capped_views", "the signed Venn of the capped-set DE genes vs TruSeq (Fig. 6B = figures/de_capped/venn_vs_truseq_small.svg) and the PCA on exactly the 7,500 capped genes", ["src/de_capped_views.py"]),
    ("scatter",      "gene-vs-gene expression of one method against each other method (method means, thinned matrix), outliers beyond 4x coloured and tabulated, the 15 per side with the largest deviation x expression named; Fig. S5A-C = figures/scatter_truseq_vs_methods.svg", ["src/scatter.py", "TruSeq", "BOBseq"]),
]
DEFAULT = ["tables", "expression", "pca", "de", "capped_views", "scatter"]
ap = argparse.ArgumentParser(description=__doc__); ap.add_argument("--only", nargs="+"); ap.add_argument("--list", action="store_true"); a = ap.parse_args()
if a.list:
    for n, d, _ in STAGES: print(f"  {n:<13} {d}")
    sys.exit(0)
for n, d, cmd in [s for s in STAGES if s[0] in (a.only or DEFAULT)]:
    print(f"== {n}: {d}", flush=True)
    if subprocess.run([sys.executable] + cmd, cwd=P).returncode: sys.exit(f"stage {n} failed")

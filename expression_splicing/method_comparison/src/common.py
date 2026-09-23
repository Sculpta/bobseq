"""method_comparison: BOBseq's three untreated control wells against every control sample of prime-seq (8), BRB-seq (8) and
TruSeq (10): 29 samples on the mRNA gene set, depth-matched, with the two full-length libraries LENGTH-NORMALISED: BOBseq and
TruSeq expression is TPM (counts / Salmon effective gene length, per million), prime-seq and BRB-seq stay CPM (a 3'-tag count
carries no length term). TPM is used for the PCA and the scatters; DESeq2 stays on counts and takes the same lengths as per-gene
normalisation factors. Shared facts for the stages: locations, which methods are length-normalised, the method palette, the
depth-thinning routine, small I/O helpers and the matplotlib style. The sample list is in counts.py; the gene set (the mRNA gene
set of the preprint, S10_mrna_gene_set.tsv: protein-coding minus MT-encoded and ribosomal-protein genes, haplotype copies
collapsed) is applied by mrna.py or tables.py; expression.py converts."""
import csv, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); ANALYSIS = os.path.dirname(HERE)
sys.path.insert(0, os.path.dirname(ANALYSIS))
from locations import GTF, table, open_text                                   # noqa: E402,F401
RESULTS = os.path.join(ANALYSIS, "results"); FIGURES = os.path.join(ANALYSIS, "figures")
SALMON_QUANT = os.environ.get("BOBSEQ_SALMON_QUANT", os.path.join(ANALYSIS, "data", "salmon"))   # <source>/quant.sf per sample (stage counts, the quant.sf route)
TX2GENE = os.environ.get("BOBSEQ_TX2GENE", os.path.join(ANALYSIS, "data", "homo_sapiens_core_113_38_tx2gene.tsv"))   # TXNAME GENEID SYMBOL (stage counts)
MRNA_LENGTH = os.path.join(RESULTS, "mrna_length.tsv")                        # gene_id, length_nt: see mrna_length()
MRNA = os.path.join(RESULTS, "gene_counts_mrna.csv.gz"); MRNA_THIN = os.path.join(RESULTS, "gene_counts_mrna_thinned.csv.gz")
EFFLEN = os.path.join(RESULTS, "gene_efflen_mrna.csv.gz")                    # per-sample gene effective length (nt), stage mrna
EXPR = os.path.join(RESULTS, "expression_mrna.csv.gz"); EXPR_THIN = os.path.join(RESULTS, "expression_mrna_thinned.csv.gz")   # stage expression
TPM_METHODS = ("BOBseq", "TruSeq")                # length-normalised (TPM); prime-seq and BRB-seq are CPM
UNIT = {m: ("TPM" if m in TPM_METHODS else "CPM") for m in ("BOBseq", "prime-seq", "BRB-seq", "TruSeq")}
COL = {"BOBseq": "#1f78b4", "BRB-seq": "#33a02c", "prime-seq": "#ff7f00", "TruSeq": "#e31a1c"}                 # method colours, all figures

LENGTH_EDGES = [0, 1000, 2500, 5000, 10000, 1e9]; LENGTH_LABELS = ["< 1,000", "1,000-2,500", "2,500-5,000", "5,000-10,000", "> 10,000"]   # fixed mRNA-length bins (nt), as ../length_bias
# gene class by symbol (the replication-dependent histone mRNAs, a class of the outlier tables):
HISTONE_RE = r"^(H1-[1-6]|H2AC\d+|H2BC\d+|H3C\d+|H4C\d+)$"  # replication-dependent (cluster) histone genes: stem-loop, NOT polyadenylated -- the replication-independent variants (H1-0, H1-10, H2AZ1/2, H2AX, H3-3A/B ...) are polyadenylated and excluded

def mrna_length():
    """Gene -> mRNA length (nt) for the length columns of the tables (not used by any figure): the TPM-weighted mean transcript
    length in the deepest TruSeq run (SRR33776933), i.e. the length of what was quantified, not the per-sample effective length
    the TPM divides by. Read from results/mrna_length.tsv (gene_id, length_nt); computed from that run's quant.sf and the tx2gene
    table when they are present under SALMON_QUANT / TX2GENE and the file is not; otherwise every length is NaN and the length
    columns stay empty."""
    import pandas as pd
    if os.path.exists(MRNA_LENGTH):
        t = pd.read_csv(MRNA_LENGTH, sep="\t"); return pd.Series(t.length_nt.to_numpy(), index=t.gene_id)
    quant = os.path.join(SALMON_QUANT, "SRR33776933", "quant.sf")
    if not (os.path.exists(quant) and os.path.exists(TX2GENE)):
        print("mrna_length: results/mrna_length.tsv and the SRR33776933 quant.sf are absent; mRNA lengths left empty", flush=True)
        return pd.Series(dtype=float)
    q = pd.read_csv(quant, sep="\t"); tx = pd.read_csv(TX2GENE, sep="\t")
    q["gene"] = q.Name.str.split(".").str[0].map(dict(zip(tx.TXNAME.str.split(".").str[0], tx.GENEID))); q["w"] = q.TPM + 1e-3
    L = (q.Length * q.w).groupby(q.gene).sum() / q.w.groupby(q.gene).sum()
    pd.DataFrame({"gene_id": L.index, "length_nt": L.round(1).to_numpy()}).to_csv(MRNA_LENGTH, sep="\t", index=False)
    return L

def downsample(v, n, rng):
    """Thin one integer count vector to n total reads by sampling reads without replacement (hypergeometric); unchanged if already <= n."""
    v = np.asarray(v).astype(int)
    if v.sum() <= n: return v
    return np.bincount(rng.choice(np.repeat(np.arange(len(v)), v), n, replace=False), minlength=len(v))

def write_tsv(path, rows, cols):
    """Write a list of dicts as TSV (atomic replace)."""
    with open(path + ".tmp", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    os.replace(path + ".tmp", path)

def style():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib import font_manager as _fm       # Arial metrics on a box without Arial: register its
    import dataclasses as _dc                        # metric-identical twin Liberation Sans under the name Arial,
    for _p in (q for q in _fm.findSystemFonts() if "LiberationSans-" in q):   # so the layout is Arial-exact and the
        _fm.fontManager.ttflist.append(_dc.replace(_fm.ttfFontProperty(_fm.get_font(_p)), name="Arial"))   # file says Arial
    plt.rcParams.update({"font.family": "Arial",
                         "svg.fonttype": "none", "pdf.fonttype": 42, "savefig.dpi": 400, "font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
    return plt


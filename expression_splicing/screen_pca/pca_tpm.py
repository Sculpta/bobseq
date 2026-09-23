"""24-plex expression PCA on gene TPM (the human wells of set bobseq_pe_native; Salmon on the UMI-deduplicated reads,
tximport to genes). Fig. 7A is figures/pca_no_ris500.svg.

Two runs: all 21 wells, and the 18 without the QC-failed Ris 500 mM wells.
  gene TPM (S8_gene_tpm_salmon.csv.gz; effective-length corrected) -> log2(TPM + 1) -> genes with mean TPM >= 1 ->
  the run's N_HVG most variable genes -> centred, unit-variance scaled per gene -> PCA (sklearn, full SVD).
  PC1/PC2 and PC3/PC4 scatters of the sample SCORES, coloured by condition.
  Heatmaps: sample x sample Spearman correlation of log2(TPM+1) over the run's own N_HEAT most variable genes
  (re-selected per run), hierarchically clustered on both axes (average linkage, correlation distance) with dendrograms
  (seaborn clustermap); labels coloured by condition. Dendrogram-only figure too.
Writes results/pca_scores_<run>.tsv, pca_variance_<run>.tsv, correlation_<run>.tsv, genes_<run>.tsv (the selected genes);
figures/pca_<run>.svg, heatmap_<run>.svg, dendrogram_<run>.svg.
    python3 pca_tpm.py"""
import csv, os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, seaborn as sns
from sklearn.decomposition import PCA
from scipy.stats import rankdata
# Illustrator-safe text: real <text> elements in Arial. svg.fonttype "none" stops matplotlib turning every
# character into a <use> reference to an outline, which is what Illustrator mangles; pdf.fonttype 42 does the
# same for the PDFs (Type 3 otherwise).
from matplotlib import font_manager as _fm       # Arial metrics on a box without Arial: register its
import dataclasses as _dc                        # metric-identical twin Liberation Sans under the name Arial,
for _p in (q for q in _fm.findSystemFonts() if "LiberationSans-" in q):   # so the layout is Arial-exact and the
    _fm.fontManager.ttflist.append(_dc.replace(_fm.ttfFontProperty(_fm.get_font(_p)), name="Arial"))   # file says Arial
plt.rcParams.update({"font.family": "Arial",
                     "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
SC = 0.625                 # publication size: 0.625 x the working size (half, then +25 %), fonts unchanged

ANALYSIS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(ANALYSIS))
from locations import table, SAMPLES
COL = {"Hek control": "#444444", "Ris dose 25mM": "#1f77b4", "Ris dose 500mM": "#9ecae1", "CHX dose 1ug/mL": "#fdae6b", "CHX dose 50ug/mL": "#d94801", "CADM1 1nM": "#a1d99b", "CADM1 100nM": "#31a354"}
RUNS = {"all21": lambda c: True, "no_ris500": lambda c: c != "Ris dose 500mM"}
N_HVG = 5000            # PCA: per-run most variable genes
N_HEAT = 5000           # heatmaps: most variable genes, re-selected per run

def main():
    os.makedirs(os.path.join(ANALYSIS, "results"), exist_ok=True); os.makedirs(os.path.join(ANALYSIS, "figures"), exist_ok=True)
    meta = [r for r in csv.DictReader(open(SAMPLES)) if r["set"] == "bobseq_pe_native"]
    g = pd.read_csv(table("S8"))                                       # gene TPM = sum of the transcripts' Salmon TPM
    for run, keep in RUNS.items():
        W = [r for r in meta if keep(r["condition"])]; accs = [r["acc_number"] for r in W]; names = [r["sample_name"] for r in W]; cond = [r["condition"] for r in W]
        tpm = g[accs].to_numpy(float); ok = tpm.mean(1) >= 1; L = np.log2(tpm[ok] + 1)          # TPM already sums to 1e6 per well
        hv = np.argsort(-L.var(1))[:N_HVG]; X = L[hv].T; Xs = (X - X.mean(0)) / (X.std(0) + 1e-9)
        # the selected genes (the table behind the PCA and the heatmap): id, symbol, mean TPM, variance of log2(TPM+1)
        gid, gsym = g["GENEID"].to_numpy()[ok][hv], g["SYMBOL"].fillna("").to_numpy()[ok][hv]
        pd.DataFrame({"rank": np.arange(1, len(hv) + 1), "gene_id": gid, "symbol": gsym, "mean_tpm": tpm[ok][hv].mean(1).round(3), "var_log2tpm": L[hv].var(1).round(4)}).to_csv(
            os.path.join(ANALYSIS, "results", f"genes_{run}.tsv"), sep="\t", index=False)
        p = PCA(n_components=min(10, len(W) - 1), svd_solver="full", random_state=0); S = p.fit_transform(Xs); v = p.explained_variance_ratio_
        pd.DataFrame(S, index=names, columns=[f"PC{i + 1}" for i in range(S.shape[1])]).assign(condition=cond).to_csv(os.path.join(ANALYSIS, "results", f"pca_scores_{run}.tsv"), sep="\t")
        pd.DataFrame({"pc": [f"PC{i + 1}" for i in range(len(v))], "var_ratio": v}).to_csv(os.path.join(ANALYSIS, "results", f"pca_variance_{run}.tsv"), sep="\t", index=False)
        # ---- PCA scatters
        fig, axes = plt.subplots(1, 2, figsize=(9.5 * SC, 4.6 * SC))
        for ax, (a, b) in zip(axes, ((0, 1), (2, 3))):          # PC1/PC2 and PC3/PC4
            for c in COL:
                m = np.array([x == c for x in cond])
                if m.any(): ax.scatter(S[m, a], S[m, b], s=48, color=COL[c], label=c, zorder=3, alpha=0.7, edgecolor="white", lw=0.5)
            ax.set_xlabel(f"PC{a + 1} ({100 * v[a]:.0f} %)"); ax.set_ylabel(f"PC{b + 1} ({100 * v[b]:.0f} %)"); ax.set_box_aspect(1)
        axes[1].legend(fontsize=6, frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5), handletextpad=0.2, borderpad=0.2, labelspacing=0.3)
        fig.suptitle(f"PCA, log2 TPM, top {N_HVG} variable genes — {len(W)} wells" + (" (no Ris 500)" if run == "no_ris500" else ""), fontsize=9)
        fig.tight_layout(); [fig.savefig(os.path.join(ANALYSIS, "figures", f"pca_{run}.{e}")) for e in ("svg",)]; plt.close(fig)
        # ---- clustered correlation heatmap
        # correlation on the run's own N_HEAT most variable genes (re-selected per run)
        hvh = np.argsort(-L.var(1))[:N_HEAT]; XE = L[hvh].T
        R = pd.DataFrame(np.corrcoef(rankdata(XE, axis=1)), index=names, columns=names)     # Spearman = Pearson on average ranks (ties share a rank)
        R.round(4).to_csv(os.path.join(ANALYSIS, "results", f"correlation_{run}.tsv"), sep="\t")
        colors = pd.Series([COL[c] for c in cond], index=names, name="condition")
        cg = sns.clustermap(R, method="average", metric="correlation", cmap="Reds", vmin=0, vmax=1, row_colors=colors, col_colors=colors,
                            figsize=(8.5 * SC, 8.5 * SC), dendrogram_ratio=0.12, cbar_pos=(0.02, 0.8, 0.03, 0.15), cbar_kws={"label": "Spearman ρ"}, xticklabels=True, yticklabels=True)
        for lab in cg.ax_heatmap.get_xticklabels() + cg.ax_heatmap.get_yticklabels(): lab.set_fontsize(6); lab.set_color(COL[cond[names.index(lab.get_text())]])
        cg.fig.suptitle(f"Spearman ρ, log2 TPM, top {N_HEAT} variable genes (this run) — {len(W)} wells" + (" (no Ris 500)" if run == "no_ris500" else ""), fontsize=8, y=0.995)
        [cg.savefig(os.path.join(ANALYSIS, "figures", f"heatmap_{run}.{e}")) for e in ("svg",)]; plt.close(cg.fig)
        # ---- dendrogram only (the heatmap's own tree)
        from scipy.cluster.hierarchy import dendrogram
        fig, ax = plt.subplots(figsize=(7.5 * SC, 4.6 * SC))
        dendrogram(cg.dendrogram_col.linkage, labels=names, ax=ax, leaf_rotation=90, leaf_font_size=7, color_threshold=0, above_threshold_color="#333")
        for lab in ax.get_xticklabels(): lab.set_color(COL[cond[names.index(lab.get_text())]])
        # NOTE the distance is NOT 1 - rho: clustermap(metric="correlation") on the rho matrix takes 1 - Pearson correlation
        # of two wells' ROWS of that matrix, i.e. how alike their correlation profiles are
        ax.set_ylabel("1 \u2212 correlation of the wells' \u03c1 profiles (average linkage)", fontsize=8)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.set_title(f"Clustering, log2 TPM, top {N_HEAT} variable genes — {len(W)} wells" + (" (no Ris 500)" if run == "no_ris500" else ""), fontsize=8.5, loc="left")
        fig.tight_layout(); fig.savefig(os.path.join(ANALYSIS, "figures", f"dendrogram_{run}.svg")); plt.close(fig)
        print(f"{run}: {len(W)} wells, {ok.sum()} expressed genes, PC1 {100 * v[0]:.0f} %, PC2 {100 * v[1]:.0f} %, PC3 {100 * v[2]:.0f} %")

if __name__ == "__main__":
    main()

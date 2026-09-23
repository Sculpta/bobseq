"""Stage pca: PCA and clustered Spearman heatmap of the 29 samples (BOBseq 3 / prime-seq 8 / BRB-seq 8 / TruSeq 10) on the
expression matrices of stage `expression`, two runs: `native` (counts as quantified) and `thinned` (every sample at the smallest
mRNA library). Per run: log2(x + 1) -> genes with mean x >= 1 -> the N_HVG most variable genes (variance of log2(x + 1)) ->
centred, unit-variance scaled per gene -> PCA (sklearn, full SVD, random_state 0); PC1/PC2 and PC3/PC4 scatters of the sample
scores coloured by method (figures/pca_<run>.svg) and the same with the TruSeq runs marked by laboratory, GEO series
(pca_<run>_labs.svg; Fig. 6A = pca_thinned_labs.svg). Heatmap: sample x sample Spearman rho of log2(x + 1) over the same N_HEAT
genes, average-linkage clustering on the correlation distance (seaborn clustermap), labels coloured by method; the dendrogram
alone as well. Writes results/pca_scores_<run>.tsv, pca_variance_<run>.tsv, correlation_<run>.tsv, genes_<run>.tsv."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, FIGURES, MRNA_THIN, EXPR, EXPR_THIN, COL, UNIT, TPM_METHODS, style
plt = style(); import seaborn as sns
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import dendrogram
from scipy.stats import rankdata
SC = 0.625                              # publication size: 0.625 x the working size, fonts unchanged
N_HVG = N_HEAT = 5000
ORDER = ["BOBseq", "prime-seq", "BRB-seq", "TruSeq"]

def main():
    M = pd.read_csv(EXPR, index_col=0); sym = M["SYMBOL"].fillna(""); M = M.drop(columns=["SYMBOL"]); E_thin = pd.read_csv(EXPR_THIN, index_col=0).drop(columns=["SYMBOL"]).loc[M.index]
    depth = int(pd.read_csv(MRNA_THIN, index_col=0).drop(columns=["SYMBOL"]).sum().min())
    samples = pd.read_csv(os.path.join(RESULTS, "sample_table.tsv"), sep="\t"); method = dict(zip(samples["sample"], samples["method"])); group = dict(zip(samples["sample"], samples["group"]))
    labs = sorted({group[s_] for s_ in method if method[s_] == "TruSeq"}); LAB_MARK = dict(zip(labs, ("o", "s", "^", "D")))          # TruSeq laboratory (GEO series) -> marker
    names = list(M.columns); meth = [method[s] for s in names]
    RUNS = {"native": M.to_numpy(float), "thinned": E_thin[names].to_numpy(float)}
    unit = f"TPM ({' / '.join(TPM_METHODS)}) or CPM"
    os.makedirs(FIGURES, exist_ok=True)
    for run, expr in RUNS.items():                                     # per-million expression in the method's unit (TPM or CPM), every column sums to 1e6
        ok = expr.mean(1) >= 1; L = np.log2(expr[ok] + 1)
        hv = np.argsort(-L.var(1))[:N_HVG]; X = L[hv].T; Xs = (X - X.mean(0)) / (X.std(0) + 1e-9)
        pd.DataFrame({"rank": np.arange(1, len(hv) + 1), "gene_id": M.index.to_numpy()[ok][hv], "symbol": sym.to_numpy()[ok][hv],
                      "mean_expr": expr[ok][hv].mean(1).round(3), "var_log2expr": L[hv].var(1).round(4)}).to_csv(os.path.join(RESULTS, f"genes_{run}.tsv"), sep="\t", index=False)
        p = PCA(n_components=min(10, len(names) - 1), svd_solver="full", random_state=0); S = p.fit_transform(Xs); v = p.explained_variance_ratio_
        pd.DataFrame(S, index=names, columns=[f"PC{i + 1}" for i in range(S.shape[1])]).assign(method=meth).to_csv(os.path.join(RESULTS, f"pca_scores_{run}.tsv"), sep="\t")
        pd.DataFrame({"pc": [f"PC{i + 1}" for i in range(len(v))], "var_ratio": v}).to_csv(os.path.join(RESULTS, f"pca_variance_{run}.tsv"), sep="\t", index=False)
        tag = f"{len(names)} HEK samples (BOBseq 3 / prime-seq 8 / BRB-seq 8 / TruSeq 10), the mRNA gene set\nBOBseq + TruSeq TPM, prime-seq + BRB-seq CPM" + (f"; every sample thinned to {depth:,} mRNA reads" if run == "thinned" else "; counts as quantified")
        # ---- PCA scatters
        fig, axes = plt.subplots(1, 2, figsize=(9.5 * SC, 4.6 * SC))
        for ax, (a, b) in zip(axes, ((0, 1), (2, 3))):          # PC1/PC2 and PC3/PC4
            for m in ORDER:
                sel = np.array([x == m for x in meth])
                ax.scatter(S[sel, a], S[sel, b], s=48, color=COL[m], label=m, zorder=3, alpha=0.7, edgecolor="white", lw=0.5)
            ax.set_xlabel(f"PC{a + 1} ({100 * v[a]:.0f} %)"); ax.set_ylabel(f"PC{b + 1} ({100 * v[b]:.0f} %)"); ax.set_box_aspect(1)
        axes[1].legend(fontsize=6, frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5), handletextpad=0.2, borderpad=0.2, labelspacing=0.3)
        fig.suptitle(f"PCA, log2({unit}+1), top {N_HVG} variable genes — {tag}", fontsize=7); fig.tight_layout(rect=[0, 0, 1, 0.92])
        fig.savefig(os.path.join(FIGURES, f"pca_{run}.svg")); plt.close(fig)
        # ---- the same scores with TruSeq annotated by laboratory (marker = GEO series) -> pca_<run>_labs.svg
        fig, axes = plt.subplots(1, 2, figsize=(9.5 * SC, 4.6 * SC))
        for ax, (a, b) in zip(axes, ((0, 1), (2, 3))):
            for m in ORDER:
                if m == "TruSeq":
                    for lab in labs:
                        sel = np.array([x == m and group[n] == lab for x, n in zip(meth, names)])
                        ax.scatter(S[sel, a], S[sel, b], s=48, color=COL[m], marker=LAB_MARK[lab], label=f"TruSeq {lab}", zorder=3, alpha=0.75, edgecolor="white", lw=0.5)
                else:
                    sel = np.array([x == m for x in meth]); ax.scatter(S[sel, a], S[sel, b], s=48, color=COL[m], label=m, zorder=3, alpha=0.7, edgecolor="white", lw=0.5)
            ax.set_xlabel(f"PC{a + 1} ({100 * v[a]:.0f} %)"); ax.set_ylabel(f"PC{b + 1} ({100 * v[b]:.0f} %)"); ax.set_box_aspect(1)
        axes[1].legend(fontsize=6, frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5), handletextpad=0.2, borderpad=0.2, labelspacing=0.3)
        fig.suptitle(f"PCA, log2({unit}+1), top {N_HVG} variable genes — {tag}; TruSeq marker = laboratory (GEO series)", fontsize=7); fig.tight_layout(rect=[0, 0, 1, 0.92])
        fig.savefig(os.path.join(FIGURES, f"pca_{run}_labs.svg")); plt.close(fig)
        # ---- clustered Spearman heatmap on the same genes
        XE = L[np.argsort(-L.var(1))[:N_HEAT]].T
        R = pd.DataFrame(np.corrcoef(rankdata(XE, axis=1)), index=names, columns=names)     # Spearman = Pearson on average ranks (ties share a rank)
        R.round(4).to_csv(os.path.join(RESULTS, f"correlation_{run}.tsv"), sep="\t")
        colors = pd.Series([COL[m] for m in meth], index=names, name="method")
        cg = sns.clustermap(R, method="average", metric="correlation", cmap="Reds", vmin=0, vmax=1, row_colors=colors, col_colors=colors,
                            figsize=(8.5 * SC, 8.5 * SC), dendrogram_ratio=0.12, cbar_pos=(0.02, 0.8, 0.03, 0.15), cbar_kws={"label": "Spearman ρ"}, xticklabels=True, yticklabels=True)
        for lab in cg.ax_heatmap.get_xticklabels() + cg.ax_heatmap.get_yticklabels(): lab.set_fontsize(6); lab.set_color(COL[method[lab.get_text()]])
        cg.fig.suptitle(f"Spearman ρ, log2({unit}+1), top {N_HEAT} variable genes — {tag}", fontsize=6.5, y=0.995)
        cg.savefig(os.path.join(FIGURES, f"heatmap_{run}.svg")); plt.close(cg.fig)
        # ---- dendrogram only (the heatmap's column linkage)
        fig, ax = plt.subplots(figsize=(7.5 * SC, 4.6 * SC))
        dendrogram(cg.dendrogram_col.linkage, labels=names, ax=ax, leaf_rotation=90, leaf_font_size=7, color_threshold=0, above_threshold_color="#333")
        for lab in ax.get_xticklabels(): lab.set_color(COL[method[lab.get_text()]])
        ax.set_ylabel("correlation distance (average linkage)"); ax.set_title(f"Hierarchical clustering, log2({unit}+1) of the top {N_HEAT} variable genes — {tag}", fontsize=6.5, loc="left")
        fig.tight_layout(); fig.savefig(os.path.join(FIGURES, f"dendrogram_{run}.svg")); plt.close(fig)
        within = [R.iloc[i, j] for i in range(len(names)) for j in range(i + 1, len(names)) if meth[i] == meth[j]]
        between = [R.iloc[i, j] for i in range(len(names)) for j in range(i + 1, len(names)) if meth[i] != meth[j]]
        print(f"{run}: {ok.sum():,} expressed genes; PC1 {100 * v[0]:.0f} %, PC2 {100 * v[1]:.0f} %, PC3 {100 * v[2]:.0f} %; "
              f"within-method rho {min(within):.2f}-{max(within):.2f} (median {np.median(within):.2f}), between {min(between):.2f}-{max(between):.2f} (median {np.median(between):.2f})")
        for m in ORDER:
            w = [R.iloc[i, j] for i in range(len(names)) for j in range(i + 1, len(names)) if meth[i] == meth[j] == m]
            print(f"   {m:<10} PC1 {S[[x == m for x in meth], 0].mean():+7.1f}  PC2 {S[[x == m for x in meth], 1].mean():+7.1f}  within rho {min(w):.2f} / {np.median(w):.2f} / {max(w):.2f} (min / median / max, {len(w)} pairs)")

if __name__ == "__main__":
    main()

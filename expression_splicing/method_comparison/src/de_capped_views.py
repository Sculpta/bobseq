"""Stage capped_views: two views of the capped-set DE (stage de, results/de_capped), written under results/de_capped/ and
figures/de_capped/.
1. Signed Venn of the DE genes against TruSeq (venn_vs_truseq.svg; Fig. 6B = venn_vs_truseq_small.svg, the same at 0.6 x the
   canvas): the three contrasts with TruSeq as the reference (log2FC > 0 = higher than TruSeq, DE = padj < 0.001) as one
   three-set Venn whose members are (gene, direction) pairs, so a shared region holds only genes DE in the same direction in
   every method of the region; every region printed as up over down; the direction-discordant pairs counted under the figure.
   Areas are not proportional. venn_vs_truseq.tsv: per region up / down / total and the genes, plus the discordance counts.
2. PCA on exactly the capped gene set (pca_capped_<run>.svg, pca_capped_<run>_labs.svg; pca_capped_scores_<run>.tsv,
   pca_capped_variance_<run>.tsv): the recipe of pca.py on the 7,500 genes of results/de_capped_gene_set.tsv instead of the
   most variable genes, native and thinned runs.
The capped-set TruSeq scatter is stage scatter with --capped."""
import itertools, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, FIGURES, MRNA_THIN, EXPR, EXPR_THIN, COL, TPM_METHODS, write_tsv, style
plt = style()
from sklearn.decomposition import PCA
SC = 0.625
RD, FD = os.path.join(RESULTS, "de_capped"), os.path.join(FIGURES, "de_capped")
MAX_PADJ = 0.001
VS_TRUSEQ = {"BOBseq": "bobseq_vs_truseq", "prime-seq": "primeseq_vs_truseq", "BRB-seq": "brbseq_vs_truseq"}
ORDER = ["BOBseq", "prime-seq", "BRB-seq", "TruSeq"]

def venn3(ax, sets, names, colors, up, down, small=False):
    """Three-circle SIGNED Venn (equal circles, not area-proportional): members are (gene, direction) pairs, so a shared region
    holds only genes DE in the same direction in every method of the region; each region is printed as 'up' (red) over 'down'
    (blue). Returns {region: (up_ids, down_ids)}."""
    import matplotlib.patches as mp
    centres = {names[0]: (0.0, 0.58), names[1]: (-0.5, -0.29), names[2]: (0.5, -0.29)}; r = 0.78
    if small: centres = {names[0]: (0.0, 0.6), names[1]: (-0.48, -0.27), names[2]: (0.48, -0.27)}; r = 0.92      # larger overlaps: room for the same-size labels
    for n in names: ax.add_patch(mp.Circle(centres[n], r, facecolor=colors[n], edgecolor=colors[n], alpha=0.25, lw=1.0))
    regions = {}
    for k in (1, 2, 3):
        for c in itertools.combinations(names, k):
            excl = set.intersection(*[sets[n] for n in c]) - set.union(set(), *[sets[n] for n in names if n not in c])
            regions[c] = (sorted(g for g, d in excl if d == "up"), sorted(g for g, d in excl if d == "down"))
    pos = {(names[0],): (0.0, 0.95), (names[1],): (-0.85, -0.55), (names[2],): (0.85, -0.55),
           (names[0], names[1]): (-0.47, 0.22), (names[0], names[2]): (0.47, 0.22), (names[1], names[2]): (0.0, -0.62), (names[0], names[1], names[2]): (0.0, 0.0)}
    if small: pos = {(names[0],): (0.0, 1.12), (names[1],): (-0.98, -0.62), (names[2],): (0.98, -0.62),
                     (names[0], names[1]): (-0.62, 0.42), (names[0], names[2]): (0.62, 0.42), (names[1], names[2]): (0.0, -0.84), (names[0], names[1], names[2]): (0.0, -0.1)}
    dy = 0.15 if small else 0.075
    for c, (u, d) in regions.items():
        ax.text(pos[c][0], pos[c][1] + dy, f"\u2191 {len(u):,}", ha="center", va="center", fontsize=7, color=up)
        ax.text(pos[c][0], pos[c][1] - dy, f"\u2193 {len(d):,}", ha="center", va="center", fontsize=7, color=down)
    lab = {names[0]: (0.0, 1.5), names[1]: (-1.0, -1.2), names[2]: (1.0, -1.2)}
    if small: lab = {names[0]: (0.0, 1.65), names[1]: (-1.72, -1.32), names[2]: (1.72, -1.62)}         # staggered: the two labels are wider than half the canvas
    for n in names:
        nu, nd = sum(1 for g, d in sets[n] if d == "up"), sum(1 for g, d in sets[n] if d == "down")
        ha = {names[1]: "left", names[2]: "right"}.get(n, "center") if small else "center"
        ax.text(*lab[n], f"{n} ({nu + nd:,})" if small else f"{n} ({nu + nd:,})\n\u2191 {nu:,} / \u2193 {nd:,}", ha=ha, va="center", fontsize=7, color=colors[n])
    ax.set_xlim(-1.75, 1.75); ax.set_ylim(-1.9, 1.7); ax.set_aspect("equal"); ax.axis("off")
    return regions

def venns(small=False):
    """small=True: the same Venn at 0.6 x the canvas with the same font sizes -> venn_vs_truseq_small.svg (the figure panel);
    the explanatory subtitle is shortened to fit, the numbers are identical."""
    de = {m: pd.read_csv(os.path.join(RD, c, f"de_{c}.tsv"), sep="\t").set_index("gene_id") for m, c in VS_TRUSEQ.items()}
    sym = pd.concat([d["symbol"] for d in de.values()]); sym = sym[~sym.index.duplicated()]
    names = list(VS_TRUSEQ); colors = {m: COL[m] for m in names}; UP, DN = "#d62728", "#1f77b4"
    sig = pd.DataFrame({m: d.padj < MAX_PADJ for m, d in de.items()}); fc = pd.DataFrame({m: d.log2fc for m, d in de.items()})
    sets = {m: {(g, "up" if fc.loc[g, m] > 0 else "down") for g in sig.index[sig[m]]} for m in names}     # (gene, direction) pairs
    k = 0.6 if small else 1.0; fig, ax = plt.subplots(figsize=(6.4 * SC * k, 7.4 * SC * k)); rows = []      # small = 0.6 x the canvas, fonts unchanged
    regions = venn3(ax, sets, names, colors, UP, DN, small)
    for c, (u, d) in regions.items():
        rows.append({"region": " & ".join(c), "n_methods": len(c), "up_higher_than_truseq": len(u), "down_lower_than_truseq": len(d), "total": len(u) + len(d),
                     "symbols_up": ";".join(str(x) for x in sym.reindex(u).fillna("").tolist() if x)[:4000], "symbols_down": ";".join(str(x) for x in sym.reindex(d).fillna("").tolist() if x)[:4000]})
    # direction-discordant genes: DE in both methods of a pair (any third) with opposite signs -- they sit in method-specific regions above
    disc = []
    for a, b in itertools.combinations(names, 2):
        both = sig.index[sig[a] & sig[b]]; n_disc = int((np.sign(fc.loc[both, a]) != np.sign(fc.loc[both, b])).sum())
        disc.append((a, b, len(both), n_disc)); rows.append({"region": f"discordant: {a} & {b}", "n_methods": 2, "up_higher_than_truseq": "", "down_lower_than_truseq": "", "total": n_disc,
                                                            "symbols_up": f"DE in both {a} and {b} (any third method): {len(both):,}; of these with opposite direction: {n_disc:,}", "symbols_down": ""})
    all3 = sig.index[sig.all(axis=1)]; n3 = int((np.sign(fc.loc[all3]).nunique(axis=1) > 1).sum())
    rows.append({"region": "discordant: all three", "n_methods": 3, "up_higher_than_truseq": "", "down_lower_than_truseq": "", "total": n3, "symbols_up": f"DE in all three: {len(all3):,}; not the same direction in all three: {n3:,}", "symbols_down": ""})
    n_tested = int(pd.concat([d.padj.notna() for d in de.values()], axis=1).all(axis=1).sum())
    if small:
        fig.suptitle(f"DE vs TruSeq, capped set ({n_tested:,} genes), padj < {MAX_PADJ}\nsigned Venn: \u2191 higher / \u2193 lower than TruSeq;\nshared = same direction in every method", fontsize=7, y=0.995)
        ax.text(0.0, -1.9, "opposite directions in two methods:\n" + "\n".join(f"{a} & {b} {n_disc:,} of {n:,}" for a, b, n, n_disc in disc) + f"\nall three: {n3:,} of {len(all3):,}",
                ha="center", va="top", fontsize=6, color="#444", transform=ax.transData)
        ax.set_xlim(-1.75, 1.75); ax.set_ylim(-3.3, 1.85); fig.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.84)
        fig.savefig(os.path.join(FD, "venn_vs_truseq_small.svg"), bbox_inches="tight"); plt.close(fig); return
    fig.suptitle(f"DE genes of the capped common set against TruSeq — each circle = one contrast vs TruSeq\n"
                 f"DESeq2 on the thinned counts, {n_tested:,} genes tested in all three contrasts, DE = padj < {MAX_PADJ}\n"
                 f"signed Venn: members are (gene, direction) pairs — \u2191 higher than TruSeq, \u2193 lower — so a shared region holds only genes\n"
                 f"DE in the same direction in every method of the region; a gene up in one method and down in another counts in each method's own region",
                 fontsize=7, y=0.995)
    ax.text(0.0, -1.42, "genes DE in two methods with opposite directions (counted as method-specific above):\n" +
            " · ".join(f"{a} & {b} {n_disc:,} of {n:,}" for a, b, n, n_disc in disc) + f"\nDE in all three {len(all3):,}, of which not the same direction in all three {n3:,}",
            ha="center", va="top", fontsize=6, color="#444", transform=ax.transData)
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.86); fig.savefig(os.path.join(FD, "venn_vs_truseq.svg"), bbox_inches="tight"); plt.close(fig)
    write_tsv(os.path.join(RD, "venn_vs_truseq.tsv"), rows, list(rows[0].keys()))
    print(pd.DataFrame(rows).drop(columns=["symbols_up", "symbols_down"]).to_string(index=False))

def pca_capped():
    cap = pd.read_csv(os.path.join(RESULTS, "de_capped_gene_set.tsv"), sep="\t")["gene_id"].tolist()
    M = pd.read_csv(EXPR, index_col=0).drop(columns=["SYMBOL"]).loc[cap]; E_thin = pd.read_csv(EXPR_THIN, index_col=0).drop(columns=["SYMBOL"]).loc[cap, M.columns]
    depth = int(pd.read_csv(MRNA_THIN, index_col=0).drop(columns=["SYMBOL"]).sum().min())
    samples = pd.read_csv(os.path.join(RESULTS, "sample_table.tsv"), sep="\t"); method = dict(zip(samples["sample"], samples["method"])); group = dict(zip(samples["sample"], samples["group"]))
    labs = sorted({group[s_] for s_ in method if method[s_] == "TruSeq"}); LAB_MARK = dict(zip(labs, ("o", "s", "^", "D")))
    names = list(M.columns); meth = [method[s] for s in names]; unit = f"TPM ({' / '.join(TPM_METHODS)}) or CPM"
    for run, expr in {"native": M.to_numpy(float), "thinned": E_thin.to_numpy(float)}.items():
        L = np.log2(expr + 1); X = L.T; Xs = (X - X.mean(0)) / (X.std(0) + 1e-9)
        p = PCA(n_components=min(10, len(names) - 1), svd_solver="full", random_state=0); S = p.fit_transform(Xs); v = p.explained_variance_ratio_
        pd.DataFrame(S, index=names, columns=[f"PC{i + 1}" for i in range(S.shape[1])]).assign(method=meth).to_csv(os.path.join(RD, f"pca_capped_scores_{run}.tsv"), sep="\t")
        pd.DataFrame({"pc": [f"PC{i + 1}" for i in range(len(v))], "var_ratio": v}).to_csv(os.path.join(RD, f"pca_capped_variance_{run}.tsv"), sep="\t", index=False)
        tag = f"{len(names)} HEK samples (BOBseq 3 / prime-seq 8 / BRB-seq 8 / TruSeq 10)\nBOBseq + TruSeq TPM, prime-seq + BRB-seq CPM" + (f"; every sample thinned to {depth:,} mRNA reads" if run == "thinned" else "; counts as quantified")
        for labs_fig in (False, True):
            fig, axes = plt.subplots(1, 2, figsize=(9.5 * SC, 4.6 * SC))
            for ax, (a, b) in zip(axes, ((0, 1), (2, 3))):
                for m in ORDER:
                    if labs_fig and m == "TruSeq":
                        for lab in labs:
                            sel = np.array([x == m and group[n] == lab for x, n in zip(meth, names)])
                            ax.scatter(S[sel, a], S[sel, b], s=48, color=COL[m], marker=LAB_MARK[lab], label=f"TruSeq {lab}", zorder=3, alpha=0.75, edgecolor="white", lw=0.5)
                    else:
                        sel = np.array([x == m for x in meth]); ax.scatter(S[sel, a], S[sel, b], s=48, color=COL[m], label=m, zorder=3, alpha=0.7, edgecolor="white", lw=0.5)
                ax.set_xlabel(f"PC{a + 1} ({100 * v[a]:.0f} %)"); ax.set_ylabel(f"PC{b + 1} ({100 * v[b]:.0f} %)"); ax.set_box_aspect(1)
            axes[1].legend(fontsize=6, frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5), handletextpad=0.2, borderpad=0.2, labelspacing=0.3)
            fig.suptitle(f"PCA, log2({unit}+1), the capped common set of {len(cap):,} genes (stage de) — {tag}" + ("; TruSeq marker = laboratory (GEO series)" if labs_fig else ""), fontsize=7); fig.tight_layout(rect=[0, 0, 1, 0.92])
            fig.savefig(os.path.join(FD, f"pca_capped_{run}" + ("_labs" if labs_fig else "") + ".svg")); plt.close(fig)
        print(f"pca capped {run}: PC1-4 " + " / ".join(f"{100 * x:.0f} %" for x in v[:4]))

if __name__ == "__main__":
    os.makedirs(FD, exist_ok=True); venns(); venns(small=True); pca_capped()

"""Stage scatter: gene-vs-gene expression of one FOCAL method against each of the other three, outliers coloured, the strongest named.

    python3 src/scatter.py [--capped] [FOCAL ...]    FOCAL = BOBseq (default) and / or TruSeq; --capped restricts the genes to
                                                     the capped common set of stage de and writes under results/de_capped/ and
                                                     figures/de_capped/
Fig. S5A-C = figures/scatter_truseq_vs_methods.svg (focal method TruSeq, full gene set).

Method means of log2(x + 1) on the thinned expression matrix (stage expression): x = the reference method, y = the focal method,
genes with mean expression >= MIN_EXPR in both. Outlier = |deviation from the diagonal shifted by the median of y - x| >=
OUT_LOG2 (2 = 4-fold; over-represented in the focal method red, under-represented blue). Named: the N_LABEL outliers per side
with the largest deviation x mean log2 expression, circled and listed on a rail beside the panel. The outlier sets are
intersected across the three comparisons (exclusive intersections, up and down separately) with the replication-dependent
histone mRNAs (common.HISTONE_RE) and the mRNA length (common.mrna_length, optional) as classes.
Per focal method <f> writes figures/scatter_<f>_vs_methods.svg and results/scatter_<f>_{outliers,labelled,outlier_overlap,summary}.tsv."""
import itertools, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, FIGURES, EXPR_THIN, MRNA_THIN, COL, UNIT, HISTONE_RE, mrna_length, write_tsv, style
from de import dodge
plt = style()
SC = 0.625
MIN_EXPR, OUT_LOG2, N_LABEL = 5.0, 2.0, 15
METHODS = ["BOBseq", "prime-seq", "BRB-seq", "TruSeq"]
UP, DN = "#d62728", "#1f77b4"
TAG = {"BOBseq": "bobseq", "TruSeq": "truseq", "prime-seq": "primeseq", "BRB-seq": "brbseq"}       # file-name tag of the focal method

def main(focal="BOBseq", capped=False):
    refs = [m for m in METHODS if m != focal]; tag = TAG[focal]
    RD, FD = (os.path.join(RESULTS, "de_capped"), os.path.join(FIGURES, "de_capped")) if capped else (RESULTS, FIGURES); os.makedirs(FD, exist_ok=True)
    E = pd.read_csv(EXPR_THIN, index_col=0)
    if capped: E = E.loc[pd.read_csv(os.path.join(RESULTS, "de_capped_gene_set.tsv"), sep="\t")["gene_id"]]          # the 7,500-gene common set of stage de
    sym = E["SYMBOL"].fillna(""); E = E.drop(columns=["SYMBOL"]).astype(float)
    st = pd.read_csv(os.path.join(RESULTS, "sample_table.tsv"), sep="\t"); cols = {m: list(st[st.method == m]["sample"]) for m in METHODS}
    depth = int(pd.read_csv(MRNA_THIN, index_col=0).drop(columns=["SYMBOL"]).sum().min()); ml = mrna_length(); glen = ml.reindex(E.index) if len(ml) else pd.Series(np.nan, index=E.index)
    L = {m: np.log2(E[c] + 1).mean(1) for m, c in cols.items()}; mean_expr = {m: E[c].mean(1) for m, c in cols.items()}
    hist = sym.str.match(HISTONE_RE); over, under = f"over-represented in {focal}", f"under-represented in {focal}"
    fig, axes = plt.subplots(3, 1, figsize=(8.6 * SC, 19.0 * SC), gridspec_kw={"hspace": 0.55}); rows, summary, sets, top_rows = [], [], {}, []
    lim = (0.0, float(np.ceil(max(L[m].max() for m in L))) + 0.5)
    for ax, ref in zip(axes, refs):
        ok = (mean_expr[focal] >= MIN_EXPR) & (mean_expr[ref] >= MIN_EXPR); x, y = L[ref][ok], L[focal][ok]
        off = float((y - x).median()); d = y - x - off; up, dn = d >= OUT_LOG2, d <= -OUT_LOG2
        r = float(np.corrcoef(x, y)[0, 1]); rmse = float(np.sqrt(np.mean((y - x) ** 2)))
        ax.scatter(x, y, s=3, color="#bdbdbd", lw=0, alpha=0.6, rasterized=True)
        ax.scatter(x[up], y[up], s=6, color=UP, lw=0, alpha=0.8, rasterized=True, label=f"≥ {2 ** OUT_LOG2:.0f}× {over} ({int(up.sum()):,})")
        ax.scatter(x[dn], y[dn], s=6, color=DN, lw=0, alpha=0.8, rasterized=True, label=f"≥ {2 ** OUT_LOG2:.0f}× under-represented ({int(dn.sum()):,})")
        ax.plot(lim, lim, color="#999", lw=0.7, ls="--"); ax.plot(lim, [v + off for v in lim], color="#333", lw=0.7, ls=":")
        ax.set_xlim(*lim); ax.set_ylim(*lim); ax.set_box_aspect(1)
        ax.set_xlabel(f"{ref}, mean log2({UNIT[ref]}+1)"); ax.set_ylabel(f"{focal}, mean log2({UNIT[focal]}+1)")
        ax.set_title(f"{focal} vs {ref}\n{int(ok.sum()):,} genes · r {r:.3f}\nRMSE {rmse:.2f} · median offset {off:+.2f}", fontsize=7.5, loc="left")
        # the N_LABEL outliers per side with the largest expression-weighted deviation |d| x mean log2 expression: circled, named on a
        # rail right of the panel (de.py's volcano recipe)
        score = d.abs() * (x + y) / 2; pts = []
        for side, mask, col in (("up", up, UP), ("down", dn, DN)):
            top = score[mask].sort_values(ascending=False).index[:N_LABEL]
            for k, g in enumerate(top, 1):
                pts.append((float(x[g]), float(y[g]), sym[g] or g, col))
                top_rows.append({"comparison": f"{focal} vs {ref}", "side": over if side == "up" else under, "rank": k, "gene_id": g, "symbol": sym[g],
                                 "mean_log2_reference": round(float(x[g]), 3), f"mean_log2_{tag}": round(float(y[g]), 3), "deviation_log2": round(float(d[g]), 3),
                                 "mean_log2_expression": round(float((x[g] + y[g]) / 2), 3), "score_dev_x_expr": round(float(score[g]), 2), "mrna_length_nt": int(glen[g]) if pd.notna(glen[g]) else ""})
        ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=26, facecolor="none", edgecolor="#111", lw=0.7, zorder=4)
        pts.sort(key=lambda t: -t[1]); ys = dodge([(p[1] - lim[0]) / (lim[1] - lim[0]) for p in pts], gap=1.0 / (2 * N_LABEL + 1), lo=0.005, hi=0.995)
        for (x0, y0, g, col), yl in zip(pts, ys):
            ax.annotate(g, xy=(x0, y0), xycoords="data", xytext=(1.04, yl), textcoords="axes fraction", fontsize=6.5, ha="left", va="center", color=col,
                        annotation_clip=False, arrowprops=dict(arrowstyle="-", color="#888", lw=0.5, shrinkA=0, shrinkB=2))
        ax.legend(fontsize=6, frameon=False, loc="upper left", bbox_to_anchor=(0.0, -0.13), handletextpad=0.2, borderpad=0.2, labelspacing=0.2)
        for g in d.index[up | dn]:
            rows.append({"comparison": f"{focal} vs {ref}", "gene_id": g, "symbol": sym[g], "mean_log2_reference": round(float(x[g]), 3), f"mean_log2_{tag}": round(float(y[g]), 3), "deviation_log2": round(float(d[g]), 3),
                         "direction": over if d[g] > 0 else under, "histone_replication_dependent": int(hist[g]), "mrna_length_nt": int(glen[g]) if pd.notna(glen[g]) else ""})
        sets[ref] = {"up": set(d.index[up]), "down": set(d.index[dn])}
        summary.append({"comparison": f"{focal} vs {ref}", "genes_plotted": int(ok.sum()), "pearson_r": round(r, 4), "rmse_from_identity": round(rmse, 4), "median_offset_log2": round(off, 3),
                        "outliers_over": int(up.sum()), "outliers_under": int(dn.sum()), "rule": f"mean expression >= {MIN_EXPR:g} in both; |deviation from the median-offset diagonal| >= {OUT_LOG2:g} log2"})
    n_focal = len(cols[focal]); what = "3 controls" if focal == "BOBseq" else f"{n_focal} runs, three laboratories"
    genes = f"the capped common set of {len(E):,} genes (stage de)" if capped else "the thinned matrix"
    fig.suptitle(f"Gene-vs-gene expression, {focal} ({what}) against each other method — method means on {genes} ({depth:,} reads per sample)\n"
                 f"dashed = y = x, dotted = the diagonal shifted by the median offset; red / blue = genes ≥ {2 ** OUT_LOG2:.0f}× above / below it\n"
                 f"circled and named = the {N_LABEL} per side with the largest deviation × mean log2 expression", fontsize=7, y=0.995)
    fig.subplots_adjust(left=0.12, right=0.66, bottom=0.05, top=0.9)      # lay out BEFORE the rail labels (annotation_clip=False), as de.py
    fig.savefig(os.path.join(FD, f"scatter_{tag}_vs_methods.svg"), bbox_inches="tight"); plt.close(fig)
    write_tsv(os.path.join(RD, f"scatter_{tag}_outliers.tsv"), rows, list(rows[0].keys())); write_tsv(os.path.join(RD, f"scatter_{tag}_summary.tsv"), summary, list(summary[0].keys()))
    write_tsv(os.path.join(RD, f"scatter_{tag}_labelled.tsv"), top_rows, list(top_rows[0].keys()))
    # ---- are the outliers the same genes against every reference? exclusive intersections, up and down, with the classes
    combos = [c for k in (3, 2, 1) for c in itertools.combinations(refs, k)]; orows = []
    for direction in ("up", "down"):
        S = {ref: sets[ref][direction] for ref in refs}
        for c in combos:
            inter = set.intersection(*[S[r] for r in c]) - set.union(*[S[r] for r in refs if r not in c]) if len(c) < 3 else set.intersection(*[S[r] for r in c])
            ids = sorted(inter); ln = glen.reindex(ids)
            orows.append({"direction": over if direction == "up" else under, "references": " & ".join(c), "n_references": len(c), "genes_exclusive": len(ids),
                          "histone_replication_dependent": int(hist.reindex(ids).sum()), "mrna_lt_1kb": int((ln < 1000).sum()), "mrna_1_5kb": int(((ln >= 1000) & (ln < 5000)).sum()), "mrna_gt_5kb": int((ln >= 5000).sum()),
                          "median_mrna_length_nt": int(ln.median()) if len(ids) and ln.notna().any() else "", "symbols": ";".join(s_ for s_ in sym.reindex(ids).tolist() if s_)[:2000]})
    write_tsv(os.path.join(RD, f"scatter_{tag}_outlier_overlap.tsv"), orows, list(orows[0].keys()))
    print(f"== scatter, focal {focal}" + (" (capped set)" if capped else "")); print(pd.DataFrame(summary).drop(columns=["rule"]).to_string(index=False)); print(pd.DataFrame(orows).drop(columns=["symbols"]).to_string(index=False))

if __name__ == "__main__":
    capped = "--capped" in sys.argv; args = [a for a in sys.argv[1:] if a != "--capped"]
    for focal in (args or ["BOBseq"]): main(focal, capped)

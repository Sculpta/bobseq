"""Stage barcharts: for each contrast, the top-10 junctions (results/ds_top10.tsv) as PSI bar charts over ALL conditions
except Ris 500 mM (common.BAR_ORDER): bar = pooled fragment PSI of the condition's 3 wells with a Clopper-Pearson
95 % CI, dots = wells. The contrast's own conditions are drawn full, the others faded. Every number behind a bar goes
to results/barchart_psi.tsv. figures/barcharts_<contrast>.svg."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
plt = style()

def main():
    W = [w for w in usable_wells() if w["condition"] in BAR_ORDER]; accs = [w["acc_number"] for w in W]; cond = np.array([w["condition"] for w in W])
    rn, K, T = load_counts(COUNTS_FRAG, accs); ri = {r: i for i, r in enumerate(rn)}; top = read_tsv("ds_top10.tsv"); out = []
    for tag, spec in CONTRASTS.items():
        rows = [r for r in top if r["contrast"] == tag]
        if not rows: print(f"barcharts: {tag}: nothing to plot"); continue
        own = set(spec["groups"]); ncol = 5; nrow = int(np.ceil(len(rows) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(2.4 * ncol, 2.5 * nrow + 0.6), sharey=True, squeeze=False)
        for ax in axes.flat: ax.axis("off")
        for ax, r in zip(axes.flat, rows):
            ax.axis("on"); i = ri[r["junction_id"]]
            for ci, c in enumerate(BAR_ORDER):
                m = cond == c; k, n = K[i, m].sum(), T[i, m].sum(); psi = k / n if n else np.nan; lo, hi = cp_ci(int(round(k)), int(round(n)))
                per = [K[i, j] / T[i, j] if T[i, j] else np.nan for j in np.flatnonzero(m)]; a = 0.85 if c in own else 0.3
                ax.bar(ci, psi, color=COL[c], alpha=a, width=0.72, zorder=2)
                if n: ax.errorbar(ci, psi, yerr=[[psi - lo], [hi - psi]], fmt="none", ecolor="#333", elinewidth=0.9, capsize=2.5, zorder=3)
                ax.scatter([ci + (j - (len(per) - 1) / 2) * 0.17 for j in range(len(per))], per, s=13, color="#111", zorder=4)
                out.append({"contrast": tag, "rank": r["rank"], "gene": r["gene"], "lsv": r["lsv"], "junction": r["junction"], "condition": c, "n_wells": int(m.sum()),
                            "psi": round(psi, 3) if n else "", "ci_lo": round(lo, 3) if n else "", "ci_hi": round(hi, 3) if n else "",
                            "fragments_junction": int(k), "fragments_lsv": int(n), "per_well_psi": ";".join("" if x != x else f"{x:.3f}" for x in per)})
            ax.set_xticks(range(len(BAR_ORDER))); ax.set_xticklabels([LAB[c] for c in BAR_ORDER], rotation=45, ha="right", fontsize=7); ax.set_ylim(0, 1.05); ax.set_yticks([0, 0.5, 1])
            ax.tick_params(axis="y", labelleft=True, labelsize=7)
            ax.set_title(f"{r['rank']}. {r['gene']}  ΔPSI {float(r['dpsi']):+.2f}  padj {float(r['padj']):.2g}{'' if r['hit'] == '1' else '  (ns)'}\n"
                         f"{ucsc_chrom(r['chrom'])}:{r['start']}-{r['end']} ({r['strand']})", fontsize=7, loc="left")
        for ax in axes[:, 0]: ax.set_ylabel("PSI (fragments)", fontsize=8)
        nh = sum(r["hit"] == "1" for r in rows)
        fig.suptitle(f"{spec['title']} — top {len(rows)} junctions by ΔPSI among hits (padj < {MAX_PADJ}, ΔPSI ≥ {HIT_DPSI}); {nh} hit{'s' if nh != 1 else ''}, ns = not a hit (largest-ΔPSI non-hits fill the list)\n"
                     f"bars = pooled fragment PSI per condition, 95 % Clopper-Pearson CI; dots = wells; faded = conditions outside this contrast; Ris 500 mM omitted", fontsize=8.5, x=0.01, y=0.995, ha="left", va="top")
        fig.tight_layout(rect=[0, 0, 1, 0.945], h_pad=1.4, w_pad=0.6); save(fig, f"barcharts_{tag}"); plt.close(fig)
    write_tsv(os.path.join(RESULTS, "barchart_psi.tsv"), out, list(out[0].keys()))
    print("figures: " + ", ".join(f"barcharts_{t}.svg" for t in CONTRASTS))

if __name__ == "__main__":
    main()

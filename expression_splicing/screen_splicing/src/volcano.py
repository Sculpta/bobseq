"""Stage volcano: one volcano per contrast (dPSI vs -log10 padj over every tested row), two versions:
figures/volcano_<contrast>.svg (no labels) and volcano_<contrast>_labelled.svg (the top-10 junctions labelled with the
gene symbol; labels dodged onto a rail to the right of the cloud and joined to their point by a leader line; grey = not a hit).
Hits (padj < MAX_PADJ, dPSI >= HIT_DPSI, positive) are coloured; everything else grey. Each LSV contributes a +/- pair
(members sum to 1), so every hit on the right has a mirror on the left that is not called."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
plt = style()

def dodge(desired, gap, lo, hi):
    """1-D label placement: keep each label as close to `desired` (sorted descending) as possible with >= gap between
    neighbours, inside [lo, hi]. Top-down pass, then a bottom-up push if the chain ran past `lo`."""
    y = []
    for d in desired: y.append(min(d, hi) if not y else min(d, y[-1] - gap))
    if y and y[-1] < lo:
        y[-1] = lo
        for i in range(len(y) - 2, -1, -1): y[i] = max(y[i], y[i + 1] + gap)
    return y

def volcano(tag, spec, labelled):
    rows = read_tsv(f"ds_{tag}.tsv"); top = [r for r in read_tsv("ds_top10.tsv") if r["contrast"] == tag]
    d = np.array([float(r["dpsi"]) for r in rows]); q = np.array([float(r["padj"]) for r in rows]); y = -np.log10(q)
    hit = np.array([r["hit"] == "1" for r in rows]); XL = max(0.5, np.abs(d).max() * 1.1); YL = max(1.5, y.max() * 1.12)
    fig, ax = plt.subplots(figsize=(4.6, 4.2) if not labelled else (5.7, 4.2))
    ax.scatter(d[~hit], y[~hit], s=7, color="#bdbdbd", lw=0, alpha=0.7, zorder=2, rasterized=True)
    ax.scatter(d[hit], y[hit], s=18, color=spec["color"], lw=0, alpha=0.7, zorder=3, rasterized=True)
    ax.axvline(HIT_DPSI, ls=":", color="#999", lw=0.8); ax.axvline(-HIT_DPSI, ls=":", color="#999", lw=0.8); ax.axhline(-np.log10(MAX_PADJ), ls=":", color="#999", lw=0.8)
    ax.set_xlim(-XL, XL); ax.set_ylim(0, YL)
    ax.set_xlabel("ΔPSI (top dose − control)" if spec["kind"] == "dose" else "ΔPSI (treated − control)")
    ax.set_ylabel("−log10 padj (OLS slope, BH)" if spec["kind"] == "dose" else "−log10 padj (Welch t, BH)")
    nh = int(hit.sum()); ax.set_title(f"{spec['title']}\n{len(rows)} junctions tested · {nh} hit{'s' if nh != 1 else ''} (padj < {MAX_PADJ}, ΔPSI ≥ {HIT_DPSI})", loc="left")
    if labelled and top:
        # labels on a rail in the right margin (axes fraction x = 1.04), dodged vertically, leader line to the point
        pts = sorted([(float(r["dpsi"]), -np.log10(float(r["padj"])), r["gene"], r["hit"] == "1") for r in top], key=lambda t: -t[1])
        ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=26, facecolor="none", edgecolor="#111", lw=0.7, zorder=4)   # ring the top 10
        ys = dodge([p[1] / YL for p in pts], gap=0.065, lo=0.03, hi=0.97)
        for (x0, y0, g, h), yl in zip(pts, ys):
            ax.annotate(g, xy=(x0, y0), xycoords="data", xytext=(1.04, yl), textcoords="axes fraction", fontsize=7.5, ha="left", va="center",
                        color="#111" if h else "#777", annotation_clip=False,
                        arrowprops=dict(arrowstyle="-", color="#888", lw=0.5, shrinkA=0, shrinkB=2))
        ax.text(1.04, -0.01, "labels: top 10 by ΔPSI\n(grey = not a hit)", transform=ax.transAxes, fontsize=6.5, ha="left", va="top", color="#777")
    fig.tight_layout(rect=[0, 0, 0.80 if labelled else 1, 1]); name = f"volcano_{tag}" + ("_labelled" if labelled else ""); save(fig, name); plt.close(fig); return name

def main():
    names = [volcano(tag, spec, lab) for tag, spec in CONTRASTS.items() for lab in (False, True)]
    print("figures: " + ", ".join(n + ".svg" for n in names))

if __name__ == "__main__":
    main()

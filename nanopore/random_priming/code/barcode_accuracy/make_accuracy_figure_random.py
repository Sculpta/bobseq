#!/usr/bin/env python3
"""Barcode accuracy figure: one point per dataset, four panels (bobcode chemistries).

Reads only the pipeline output (bobcode_metrics.tsv from speciesmix_bobcode_metrics.py)
and the run rules (for sequencing dates). Writes

  <out-dir>/barcode_accuracy_all_chemistries_log.pdf / .png
  <out-dir>/barcode_accuracy_index.tsv   every point: chemistry, number, run, barcode,
                                          accuracy, reads scored, 95 % Wilson CI

Layout: panels left to right poly(dT), splint, TSO 5'DBCO, TSO 3'DBCO; points ranked by
accuracy within each panel and numbered from 1 in each panel. Y axis: accuracy on a
flipped log-of-inaccuracy scale, 99.5 % (top) to 40 % (bottom); points beyond are drawn
as arrowheads at the edge. Error bars: 95 % Wilson interval on the scored reads.
Bands: > 99 % darker green (bobcode benchmark), 96-99 % light green (DRUG-seq
benchmark), 90-96 % pale yellow-orange, < 90 % very light red-orange; grey rules every
10 %, 50 % (chance) in light red. Sized 9.5 x 3.45 in to sit sideways on a portrait
page; all text >= 7 pt, Arial. Datasets the pipeline marks EXCLUDE (< 100 scored mRNA
reads) are left out.

Usage:
  python3 make_accuracy_figure.py --results results/bobcode_metrics.tsv --out-dir figures
Requires matplotlib.
"""
import argparse
import csv
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

HERE = os.path.dirname(os.path.abspath(__file__))

ORDER = ("polydT", "splint", "5'TSO", "3'TSO", "randomPriming")                 # panels, left to right
TITLE = {"splint": "splint-mediated ligation\nof bobcode to cDNA",   # narrow panels wrap
         "polydT": "poly(dT) priming bobcode",
         "5'TSO":  "TSO bobcode\n(5'DBCO)",
         "3'TSO":  "TSO bobcode (3'DBCO)",
         "randomPriming": "random-priming bobcode\n(3'TSO)"}
SHORT = {"splint": "splint", "polydT": "poly(dT)", "5'TSO": "TSO (5'DBCO)", "3'TSO": "TSO (3'DBCO)", "randomPriming": "random priming"}

DOT = "#111111"
BAND = "#D4EBDA"          # > 99 % accuracy (bobcode benchmark)
BAND_MID = "#EDF7EF"      # 96-99 % (DRUG-seq benchmark)
BAND_LOW = "#FFF8EB"      # 90-96 %
BAND_POOR = "#FFF9F7"     # < 90 %
BAND_INK, BAND_MID_INK = "#2C6B47", "#3C7A52"
LO, ACC_FLOOR = 0.5, 40.0   # y axis in inaccuracy: 0.5 % (99.5 % accuracy) .. 60 % (40 %)
FIG_W, FIG_H = 1.9, 3.45
GAP_UNITS = 1               # white gap between panels, in point slots

plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "ytick.minor.width": 0.5,
})


def wilson(k, n, z=1.96):
    """95 % Wilson score interval for k of n, as percentages."""
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h) * 100, min(1.0, c + h) * 100


def load(results, run_rules):
    """{chemistry: [(run, date, barcode, inaccuracy %, reads scored)]}, worst first."""
    dates = {r["run_folder"].split("/")[-1]: r["seq_date"]
             for r in csv.DictReader(open(run_rules), delimiter="\t")}
    data = {c: [] for c in ORDER}
    for x in csv.DictReader(open(results), delimiter="\t"):
        if x["chemistry"] not in data or x["status"].startswith("EXCLUDE") or not x["total_mRNA_accuracy_pct"]:
            continue
        data[x["chemistry"]].append((x["run"], dates[x["run"]], x["barcode"][2:],
                                     100 - float(x["total_mRNA_accuracy_pct"]), int(x["total_mRNA_scored"])))
    return [(c, sorted(data[c], key=lambda t: t[3], reverse=True)) for c in ORDER if data[c]]


def render(DATA, out):
    NS = [len(v) for _, v in DATA]
    fig, axes = plt.subplots(1, len(DATA), figsize=(FIG_W, FIG_H), sharey=True,
                             gridspec_kw=dict(width_ratios=NS, wspace=GAP_UNITS * len(NS) / sum(NS)))
    if len(DATA) == 1:
        axes = [axes]
    lo, hi = LO, 100 - ACC_FLOOR
    arrows, floors = [], []
    for pi, (ax, (name, recs)) in enumerate(zip(axes, DATA)):
        n = len(recs)
        for i, (_run, _date, _bc, v, nreads) in enumerate(recs):
            style = dict(mfc=DOT, mec=DOT, mew=0)
            if v > hi:                                   # below the 40 % floor
                ax.plot(i, hi, "v", ms=3.4, clip_on=False, zorder=3, **style)
                floors.append((name, i + 1, v))
            elif v < lo:                                 # above the 99.5 % ceiling
                ax.plot(i, lo, "^", ms=3.4, clip_on=False, zorder=3, **style)
                arrows.append((name, i + 1, v))
            else:
                ax.plot(i, v, "o", ms=3.0, zorder=3, **style)
                blo, bhi = wilson(round(v / 100 * nreads), nreads)      # interval on inaccuracy
                ax.plot([i, i], [max(blo, lo), min(bhi, hi)], lw=1.0,
                        color="black", zorder=2, solid_capstyle="butt")
        ax.axhspan(hi, 10, color=BAND_POOR, lw=0, zorder=0)
        ax.axhspan(10, 4, color=BAND_LOW, lw=0, zorder=0)
        ax.axhspan(1, 4, color=BAND_MID, lw=0, zorder=0)
        ax.axhspan(lo, 1, color=BAND, lw=0, zorder=0)
        for acc in range(90, int(ACC_FLOOR), -10):        # grey rule every 10 %; 50 % in light red
            ax.axhline(100 - acc, color=("#F2A0A0" if acc == 50 else "0.82"), lw=0.6, zorder=1)
        if pi == 0:                                       # band labels, left of the first panel
            ax.text(0.3, 3.85, "DRUG-seq benchmark (>96%)", ha="left", va="bottom",
                    fontsize=4.5, color=BAND_MID_INK, zorder=2)
            ax.text(0.3, 0.965, "bobcode benchmark (>99%)", ha="left", va="bottom",
                    fontsize=4.5, color=BAND_INK, zorder=2)
        ax.set_xlim(-0.5, n - 0.5)
        major = [i for i in range(n) if i == 0 or (i + 1) % 5 == 0 or i == n - 1]   # 1, 5, 10, ... last
        ax.set_xticks(major)
        ax.set_xticks([i for i in range(n) if i not in major], minor=True)
        ax.set_xticklabels([str(i + 1) for i in major], rotation=90, ha="center", va="top", fontsize=7.0)
        ax.set_title(TITLE[name], fontsize=7, pad=4)
        for s in (("top", "right") if pi == 0 else ("top", "right", "left")):
            ax.spines[s].set_visible(False)
        if pi:
            ax.tick_params(axis="y", which="both", length=0)
        ax.tick_params(axis="x", which="major", length=3.5, width=0.7, pad=2)
        ax.tick_params(axis="x", which="minor", length=2, width=0.6)
        ax.set_axisbelow(True)

    ax0 = axes[0]
    ax0.set_yscale("log")
    ax0.set_ylim(hi, lo)                                  # flipped: best at the top
    ax0.set_yticks([LO, 1, 10, 20, 30, 40, 50])           # 99.5, 99, 90, 80, 70, 60, 50
    ax0.set_yticks([100 - ACC_FLOOR], minor=True)         # 40 %: unlabelled (would collide with 50)
    ax0.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{100 - v:g}"))
    ax0.yaxis.set_minor_formatter(FuncFormatter(lambda v, _: ""))
    ax0.tick_params(axis="y", which="minor", length=2.5)
    ax0.tick_params(axis="y", which="major", labelsize=10)
    ax0.set_ylabel("Barcode accuracy (%)", fontsize=10)
    fig.subplots_adjust(left=0.80 / FIG_W, right=1 - 0.06 / FIG_W, top=1 - 0.50 / FIG_H, bottom=0.33 / FIG_H)
    for lst, what in ((arrows, "above the 99.5 % ceiling"), (floors, "below the 40 % floor")):
        if lst:
            print(f"{len(lst)} {what}: " + ", ".join(f"{SHORT[c]} #{r}" for c, r, _ in lst))
    fig.savefig(out + ".pdf")
    fig.savefig(out + ".png", dpi=300)
    plt.close(fig)


def write_index(DATA, path):
    with open(path, "w") as fh:
        fh.write("n\tchemistry\trun\tdate\tbarcode\taccuracy_pct\treads_scored\tci95_lo_pct\tci95_hi_pct\n")
        for name, recs in DATA:
            for i, (run, date, bc, v, nreads) in enumerate(recs):
                ilo, ihi = wilson(round(v / 100 * nreads), nreads)
                fh.write(f"{i + 1}\t{name}\t{run}\t{date}\tBC{bc}\t{100 - v:.2f}\t{nreads}\t"
                         f"{100 - ihi:.2f}\t{100 - ilo:.2f}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True, help="bobcode_metrics.tsv from speciesmix_bobcode_metrics.py")
    ap.add_argument("--run-rules", default=os.path.join(HERE, "bobcode_run_rules.tsv"))
    ap.add_argument("--out-dir", default=".")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    DATA = load(a.results, a.run_rules)
    render(DATA, os.path.join(a.out_dir, "barcode_accuracy_all_chemistries_log"))
    write_index(DATA, os.path.join(a.out_dir, "barcode_accuracy_index.tsv"))
    for name, recs in DATA:
        accs = sorted(100 - t[3] for t in recs)
        print(f"  {SHORT[name]:13s} n={len(recs):3d}  median {accs[len(accs) // 2]:.2f} %  "
              f"range {accs[0]:.2f}-{accs[-1]:.2f} %")


if __name__ == "__main__":
    main()

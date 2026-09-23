"""Stage figures. UMI-deduplicated junction reads, BOBseq PE, all wells of a
condition pooled (bars), per well (dots). SC = 0.625 publication scale, fonts unchanged.

figures/SMN_inclusion_ratio.svg   ratio of inclusion = SMN2 ex6>7 reads / SMN1 ex6>7 reads, per condition (bar, 95 % CI of a
                                  ratio of two Poisson counts, log-scale Katz interval) and per well (dots). NOTE: the ratio moves
                                  with SMN2 exon-7 inclusion AND with the SMN2/SMN1 expression ratio; it is not a PSI.
figures/SMN_junction_reads.svg    junction reads per event, stacked per condition (as CADM1's middle panel): SMN2 ex6>7 / ex7>8,
                                  SMN1 ex6>7 / ex7>8, and the shared ex6>8 skip (paralog-unassignable).
figures/SMN_junction_reads_pct.svg  ONE stacked bar per condition over the same events, as a PERCENTAGE of all exon-7
                                  junction reads of that condition: the per-paralog inclusion reads and the shared
                                  ex6>8 skip in the same bar, adding up to 100 %, so the depth difference between
                                  conditions drops out. Legends sit below their panel."""
import logging
import csv, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
# Illustrator-safe text: real <text> elements in Arial. svg.fonttype "none" stops matplotlib turning every
# character into a <use> reference to an outline, which is what Illustrator mangles; pdf.fonttype 42 does the
# same for the PDFs (Type 3 otherwise).
from matplotlib import font_manager as _fm       # Arial metrics on a box without Arial: register its
import dataclasses as _dc                        # metric-identical twin Liberation Sans under the name Arial,
for _p in (q for q in _fm.findSystemFonts() if "LiberationSans-" in q):   # so the layout is Arial-exact and the
    _fm.fontManager.ttflist.append(_dc.replace(_fm.ttfFontProperty(_fm.get_font(_p)), name="Arial"))   # file says Arial
plt.rcParams.update({"font.family": "Arial",
                     "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
SC = 0.625
# colour convention:
#   SMN2 red / SMN1 blue / unassignable purple = exon-7 inclusion; grey = junctions that skip it
RED, BLUE, GREY = "#d62728", "#1f77b4", "#999999"      # SMN2 red, SMN1 blue
AMB = ["#762a83", "#c2a5cf"]      # inclusion reads that reach no diagnostic base: purple, NOT grey (grey = skipping) and
                                  # not the blue-purple that sat next to the SMN1 blue in the stack
ORDER = ["Hek control", "Ris dose 25mM", "CHX dose 1ug/mL", "CHX dose 50ug/mL", "CADM1 1nM", "CADM1 100nM"]      # Ris 500 mM (dead wells) excluded from the figures
SHORT = {"Hek control": "control", "Ris dose 25mM": "Ris 25", "Ris dose 500mM": "Ris 500", "CHX dose 1ug/mL": "CHX 1", "CHX dose 50ug/mL": "CHX 50", "CADM1 1nM": "CADM1 1", "CADM1 100nM": "CADM1 100"}
UNIT = "dedup"

def rows(name): return list(csv.DictReader(open(os.path.join(ANALYSIS, "results", name)), delimiter="\t"))

def ratio_ci(a, b):
    """a / b for two Poisson counts; 95 % CI on the log ratio (Katz): exp(log(a/b) +- 1.96 sqrt(1/a + 1/b))."""
    if a == 0 or b == 0: return float("nan"), float("nan"), float("nan")
    r = a / b; se = np.sqrt(1 / a + 1 / b); return r, r * np.exp(-1.96 * se), r * np.exp(1.96 * se)

def legends_below(fig, axes, band=0.14):
    """One legend per panel, centred under that panel inside the band tight_layout reserved at the bottom -- so it never
    overlaps the axes or the rotated tick labels."""
    for ax in axes:
        h, l = ax.get_legend_handles_labels()
        if not h: continue
        pos = ax.get_position()
        fig.legend(h, l, loc="upper center", bbox_to_anchor=(pos.x0 + pos.width / 2, band), ncol=min(3, len(l)),
                   frameon=False, fontsize=6.5, handlelength=1.2, handletextpad=0.5, columnspacing=1.2, borderpad=0.1)


def main():
    cond = {r["condition"]: r for r in rows("SMN2_psi_conditions.tsv") if r["set"] == "bobseq_pe_native" and r["unit"] == UNIT}
    per = [r for r in rows("SMN2_psi_samples.tsv") if r["set"] == "bobseq_pe_native" and r["unit"] == UNIT]
    # ---- ratio of inclusion
    out = []
    fig, ax = plt.subplots(figsize=(6.0 * SC, 4.6 * SC))
    for i, c in enumerate(ORDER):
        a, b = int(cond[c]["ex6>7_SMN2"]), int(cond[c]["ex6>7_SMN1"]); r, lo, hi = ratio_ci(a, b)
        ax.bar(i, r if r == r else 0, color=RED, alpha=0.75, width=0.7)
        if r == r: ax.errorbar(i, r, yerr=[[r - lo], [hi - r]], fmt="none", ecolor="#333", capsize=3, lw=1)
        pts = [ratio_ci(int(w["ex6>7_SMN2"]), int(w["ex6>7_SMN1"]))[0] for w in per if w["condition"] == c]; pts = [p for p in pts if p == p]
        ax.scatter([i + (j - (len(pts) - 1) / 2) * 0.15 for j in range(len(pts))], pts, s=14, color="#111", zorder=3)
        out.append({"condition": c, "SMN2_ex6>7": a, "SMN1_ex6>7": b, "ratio": round(r, 3) if r == r else "", "ci_lo": round(lo, 3) if r == r else "", "ci_hi": round(hi, 3) if r == r else "",
                    "per_well": ";".join(f"{p:.3f}" for p in pts)})
    ax.set_xticks(range(len(ORDER))); ax.set_xticklabels([f"{SHORT[c]}\n{cond[c]['ex6>7_SMN2']} / {cond[c]['ex6>7_SMN1']}" for c in ORDER], rotation=35, ha="right", fontsize=7.5)
    ax.set_ylim(0, 2.5); ax.set_ylabel("SMN2 / SMN1 ex6>7 reads"); ax.set_xlabel("condition\nSMN2 / SMN1 ex6>7 reads", fontsize=7.5)
    ax.set_title("SMN2 : SMN1 inclusion ratio (ex6>7)\nbar = pooled wells, 95 % CI, dots = wells", fontsize=8.5, loc="left")
    fig.tight_layout(); fig.savefig(os.path.join(ANALYSIS, "figures", "SMN_inclusion_ratio.svg")); plt.close(fig)
    write_tsv(os.path.join(ANALYSIS, "results", "SMN_inclusion_ratio.tsv"), out, list(out[0].keys()))
    # ---- junction reads per event (stacked), inclusion per paralog + the shared skip; and the same as % of all
    # exon-7 junction reads of the condition (inclusion % + skip % = 100 %)
    events = [("ex6>7_SMN2", "SMN2 ex6>7", "#d62728"), ("ex7>8_SMN2", "SMN2 ex7>8", "#f4a3a3"), ("ex6>7_SMN1", "SMN1 ex6>7", "#1f77b4"), ("ex7>8_SMN1", "SMN1 ex7>8", "#9ecae1"),
              ("ex6>7_ambiguous", "ex6>7 ambiguous", AMB[0]), ("ex7>8_ambiguous", "ex7>8 ambiguous", AMB[1])]
    skip = np.array([int(cond[c]["ex6>8_SMN2"]) + int(cond[c]["ex6>8_SMN1"]) + int(cond[c]["ex6>8_ambiguous"]) for c in ORDER], float)
    incl = np.array([[int(cond[c][k]) for c in ORDER] for k, _, _ in events], float)
    total = incl.sum(0) + skip
    # ---- counts: inclusion (stacked, per paralog) and the shared skip
    fig, axes = plt.subplots(1, 2, figsize=(11 * SC, 5.2 * SC), gridspec_kw={"width_ratios": [1.6, 1]})
    ax = axes[0]; bottom = np.zeros(len(ORDER))
    for row, (k, lab, col) in zip(incl, events):
        ax.bar(range(len(ORDER)), row, bottom=bottom, color=col, width=0.7, label=lab); bottom += row
    ax.set_xticks(range(len(ORDER))); ax.set_xticklabels([SHORT[c] for c in ORDER], rotation=35, ha="right", fontsize=7.5)
    ax.set_ylabel("junction reads"); ax.set_title("inclusion junction reads\n(into / out of exon 7)", fontsize=8.5, loc="left")
    ax = axes[1]
    ax.bar(range(len(ORDER)), skip, color=GREY, width=0.7, label="ex6>8 (shared)")
    ax.set_xticks(range(len(ORDER))); ax.set_xticklabels([SHORT[c] for c in ORDER], rotation=35, ha="right", fontsize=7.5)
    ax.set_ylabel("junction reads"); ax.set_title("skip reads\n(ex6>8, shared)", fontsize=8.5, loc="left")
    fig.tight_layout(rect=[0, 0.16, 1, 1]); legends_below(fig, axes)
    fig.savefig(os.path.join(ANALYSIS, "figures", "SMN_junction_reads.svg")); plt.close(fig)
    # ---- every event in ONE bar per condition, as a % of that condition's exon-7 junction reads (adds to 100 %)
    fig, ax = plt.subplots(figsize=(7.6 * SC, 5.4 * SC)); bottom = np.zeros(len(ORDER)); den = np.where(total > 0, total, 1)
    for row, (k, lab, col) in zip(incl, events):
        v = 100 * row / den; ax.bar(range(len(ORDER)), v, bottom=bottom, color=col, width=0.7, label=lab); bottom += v
    ax.bar(range(len(ORDER)), 100 * skip / den, bottom=bottom, color=GREY, width=0.7, label="ex6>8 skip (shared)")
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([f"{SHORT[c]}\n{int(t)}" for c, t in zip(ORDER, total)], rotation=35, ha="right", fontsize=7.5)
    ax.set_xlabel("condition\ntotal exon-7 junction reads", fontsize=7.5)
    ax.set_ylim(0, 100); ax.set_ylabel("% of exon-7 junction reads")
    ax.set_title("SMN exon 7 — every event as % of the condition's inclusion + skip\nreads (dedup); each bar adds to 100 %", fontsize=8.5, loc="left")
    fig.tight_layout(rect=[0, 0.22, 1, 1]); legends_below(fig, [ax], band=0.21)
    fig.savefig(os.path.join(ANALYSIS, "figures", "SMN_junction_reads_pct.svg")); plt.close(fig)
    print("figures: SMN_inclusion_ratio.svg, SMN_junction_reads.svg, SMN_junction_reads_pct.svg")

if __name__ == "__main__":
    main()

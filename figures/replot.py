#!/usr/bin/env python3
"""Redraw every panel of the preprint from its values table alone.

Usage:  python3 replot.py [--out DIR] [--only NAME ...]

Reads the tables in values/ (next to this script) and writes one SVG + PNG per panel into
--out (default: ./replotted). No pipeline, no BAMs, no cloud access: the values tables
are the complete input. Panel types and the table each one uses:

  depth curve (mean + 95% t-interval)   p1, p2, p5 per-sample junction panels
  pooled curve                          p5 pooled junction panels
  dots per method (bar = median)        p0c, p0d2, p3b, p6b, p7e
  stacked composition bar               p3
  coverage profile vs distance          p4m2 (native and 50 nt), S14 (all 24 samples)
  read fate (mean + 95% t-interval)     p0a
  heatmap                               p4i2
  stacked QC figure                     stacked_qc_all_samples
  gene coverage tracks                  gene examples (<gene>_depth.tsv + gene_models.tsv), when those tables are present

Requires numpy, scipy, matplotlib. Text in the SVGs stays editable (svg.fonttype none).
The palette is at the top; change it here to restyle everything at once.
"""
import argparse
import csv
import glob
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = HERE   # this script lives in figures/; the tables are in figures/values/
METHODS = ["DRUG-seq", "prime-seq", "BOBseq"]
COLOR = {"DRUG-seq": "#A8325A", "prime-seq": "#8E4E00", "BOBseq": "#2E6B32"}
COMPOSITION = [  # class, label, colour (order of the stacked bar, bottom to top)
    ("mRNA", "mRNA (protein-coding)", "#4C9A4C"),
    ("ribosomal-protein", "ribosomal-protein mRNA", "#6E6E6E"),
    ("exonic other biotype", "exonic, other biotype", "#9E9E9E"),
    ("intronic", "intronic", "#C6C6C6"),
    ("intergenic", "intergenic", "#B5379B"),
    ("rRNA", "rRNA", "#E6B422"),
    ("mitochondrial", "mitochondrial", "#2260A9"),
]
DEPTH_TICKS = ([50000, 100000, 250000, 500000, 1000000, 2000000], ["50k", "100k", "250k", "500k", "1M", "2M"])
MOLECULE_TICKS = ([25000, 50000, 100000, 250000, 500000, 1000000, 2000000], ["25k", "50k", "100k", "250k", "500k", "1M", "2M"])

plt.rcParams.update({
    "svg.fonttype": "none", "font.family": ["Arial", "Nimbus Sans", "DejaVu Sans"], "font.size": 7,
    "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.linewidth": 0.6, "xtick.direction": "out", "ytick.direction": "out",
    "axes.spines.top": False, "axes.spines.right": False, "text.color": "#231F20",
    "axes.edgecolor": "#231F20", "axes.labelcolor": "#231F20", "xtick.color": "#231F20", "ytick.color": "#231F20",
})


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def find_table(name):
    """A values table by file name in figures/values/."""
    path = os.path.join(FIG, "values", name)
    if not os.path.exists(path):
        raise FileNotFoundError(name)
    return path


def save(fig, out_dir, name):
    fig.tight_layout(pad=0.4)
    for ext, kwargs in (("svg", {}), ("png", {"dpi": 200})):
        fig.savefig(os.path.join(out_dir, f"{name}.{ext}"), bbox_inches="tight", pad_inches=0.03, **kwargs)
    plt.close(fig)
    svg = os.path.join(out_dir, f"{name}.svg")
    text = open(svg).read()
    open(svg, "w").write(re.sub(r"font-family:\s*[^;\"]+", "font-family: Arial", text))
    print("wrote", name)


def method_legend(ax, counts, loc="upper left"):
    handles = [Line2D([], [], color=COLOR[m], lw=1.3) for m in METHODS]
    labels = [f"{m} (n={counts[m]})" for m in METHODS]
    ax.legend(handles=handles, labels=labels, loc=loc, frameon=False, handlelength=1.6, borderaxespad=0.3, labelspacing=0.3)


# ---------------------------------------------------------------- depth curves
def band(ax, x, y, colour):
    """Thick mean line with a 95% t-interval ribbon; depths with fewer than 3 samples are skipped."""
    depths = [d for d in sorted(set(x)) if (x == d).sum() >= 3]
    mean = np.array([y[x == d].mean() for d in depths])
    n = np.array([(x == d).sum() for d in depths])
    half = np.array([stats.t.ppf(0.975, k - 1) * y[x == d].std(ddof=1) / np.sqrt(k) for d, k in zip(depths, n)])
    ax.fill_between(depths, mean - half, mean + half, color=colour, alpha=0.3, lw=0, zorder=2)
    ax.plot(depths, mean, color=colour, lw=1.3, zorder=3)


def depth_curve(table, depth_col, y_col, xlabel, ylabel, ticks, out_dir, name, log_y=False, legend_loc="upper left"):
    rows = read_tsv(find_table(table))
    fig, ax = plt.subplots(figsize=(2.6, 2.2))
    counts = {}
    for method in METHODS:
        sub = [r for r in rows if r["method"] == method and r[y_col] != ""]
        x = np.array([float(r[depth_col]) for r in sub])
        y = np.array([float(r[y_col]) for r in sub])
        band(ax, x, y, COLOR[method])
        counts[method] = len({r["sample"] for r in sub})
    ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
    else:
        ax.set_ylim(0, None)
    ax.set_xticks(ticks[0])
    ax.set_xticklabels(ticks[1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    method_legend(ax, counts, legend_loc)
    save(fig, out_dir, name)


def pooled_curve(table, y_col, ylabel, out_dir, name):
    rows = read_tsv(find_table(table))
    fig, ax = plt.subplots(figsize=(2.8, 2.3))
    counts = {}
    for method in METHODS:
        sub = sorted([r for r in rows if r["method"] == method], key=lambda r: int(r["depth"]))
        ax.plot([int(r["depth"]) for r in sub], [float(r[y_col]) for r in sub], color=COLOR[method], lw=1.3, marker="o", ms=2.2)
        counts[method] = int(sub[0]["samples_pooled"])
    ax.set_xscale("log")
    ax.set_xticks([1e5, 1e6, 1e7, 1e8])
    ax.set_xticklabels(["100k", "1M", "10M", "100M"])
    ax.set_ylim(0, None)
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v / 1e3:g}k" if v >= 1e3 else f"{v:g}"))
    ax.set_xlabel("deduplicated molecules\n(all samples of a method pooled, subsampled)")
    ax.set_ylabel(ylabel)
    method_legend(ax, counts)
    save(fig, out_dir, name)


# ---------------------------------------------------------------- dots per method
def strip(table, y_col, ylabel, out_dir, name, log_y=False, scale=1.0, ylim=None, yticks=None):
    rows = read_tsv(find_table(table))
    fig, ax = plt.subplots(figsize=(2.2, 2.2))
    for i, method in enumerate(METHODS):
        values = np.array([float(r[y_col]) * scale for r in rows if r["method"] == method and r[y_col] != ""])
        jitter = np.random.default_rng(1).normal(0, 0.07, len(values))
        ax.scatter(i + jitter, values, s=7, color=COLOR[method], alpha=0.8, lw=0, zorder=2)
        ax.hlines(np.median(values), i - 0.28, i + 0.28, color="#231F20", lw=1.0, zorder=3)
    if log_y:
        ax.set_yscale("log")
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    if yticks:
        ax.set_yticks(yticks[0])
        ax.set_yticklabels(yticks[1])
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_xticks(range(len(METHODS)))
    ax.set_xticklabels(METHODS, rotation=30, ha="right")
    ax.set_xlim(-0.6, len(METHODS) - 0.4)
    ax.set_ylabel(ylabel)
    save(fig, out_dir, name)


# ---------------------------------------------------------------- composition
def composition(out_dir, name="p3_read_composition_native"):
    rows = {r["method"]: r for r in read_tsv(find_table("p3_read_composition_native.tsv"))}
    fig, ax = plt.subplots(figsize=(2.6, 2.2))
    left = np.zeros(len(METHODS))
    for key, label, colour in COMPOSITION:
        values = np.array([float(rows[m][key]) for m in METHODS])
        ax.barh(range(len(METHODS)), values, left=left, color=colour, height=0.7, lw=0, label=label)
        left += values
    ax.set_yticks(range(len(METHODS)))
    ax.set_yticklabels(METHODS)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of mapped reads (multimappers included)")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2, frameon=False, fontsize=5.5)
    save(fig, out_dir, name)


# ---------------------------------------------------------------- coverage profiles
def read_fate(out_dir, name="p0a_read_fate_native"):
    """p0a: % of a sample's reads into STAR surviving to mapped, uniquely mapped and filtered; mean over the samples of a
    method with a 95% t-interval ribbon (from the per-sample table, one list of four values per sample)."""
    rows = read_tsv(find_table("p0a_read_fate_native_per_sample.tsv"))
    steps = ["reads into STAR", "mapped", "uniquely mapped", "filtered"]
    fig, ax = plt.subplots(figsize=(2.6, 2.2))
    counts = {}
    for method in METHODS:
        vals = np.array([[float(v) for v in re.findall(r"[-0-9.eE+]+", r["value"])] for r in rows if r["method"] == method])
        counts[method] = len(vals)
        mean = vals.mean(axis=0)
        half = stats.t.ppf(0.975, len(vals) - 1) * vals.std(axis=0, ddof=1) / np.sqrt(len(vals))
        x = np.arange(len(steps))
        ax.fill_between(x, mean - half, mean + half, color=COLOR[method], alpha=0.3, lw=0, zorder=2)
        ax.plot(x, mean, color=COLOR[method], lw=1.3, zorder=3)
    ax.set_xticks(range(len(steps)))
    ax.set_xticklabels(steps, rotation=30, ha="right")
    ax.set_ylim(0, 100)
    ax.set_ylabel("% of reads into STAR\n(mean over samples, 95% interval)")
    method_legend(ax, counts, loc="lower left")
    save(fig, out_dir, name)


def profile_p4m2(table, out_dir, name, note=None):
    rows = read_tsv(find_table(table))
    x = np.array([float(r["nt_from_polyA_site"]) for r in rows])
    fig, ax = plt.subplots(figsize=(2.6, 2.2))
    for method in METHODS:
        y = np.array([float(r[f"{method}_mean_over_genes_smoothed25"]) for r in rows])
        ax.plot(x, y, color=COLOR[method], lw=1.2)
    ax.axhline(1, color="#9a9a9a", lw=0.6, ls=(0, (2, 2)))
    ax.set_xlim(3000, 0)
    ax.set_ylim(0, None)
    ax.set_xlabel("distance from the poly(A) site (nt)")
    ax.set_ylabel("coverage relative to gene mean\n(mean over genes, transcripts >= 1 kb,\nribosomal-protein genes excluded)" + (f"\n{note}" if note else ""))
    save(fig, out_dir, name)


def profile_all24(out_dir, name="S14_coverage_vs_nt_from_polyA_all24_native"):
    rows = read_tsv(find_table("polyA_profiles_all24_native.tsv"))
    groups = {"Hek control": "#2E6B32", "CADM1 1nM": "#4C78A8", "CADM1 100nM": "#1F4E79", "CHX dose 1ug mL": "#E6B422",
              "CHX dose 50ug mL": "#B8860B", "Ris dose 25mM": "#B5379B", "Ris dose 500mM": "#8C8C8C"}
    fig, ax = plt.subplots(figsize=(2.8, 2.3))
    x = np.arange(3000)
    for r in rows:
        y = np.array([float(r[f"nt{d}"]) for d in range(3000)])
        if r["species"] == "mouse":
            ax.plot(x, y, color="#7A4FA3", lw=0.9)
            continue
        group = re.sub(r"[-_]\d+$", "", r["sample"]).replace("_", " ")
        ax.plot(x, y, color=groups.get(group, "#444"), lw=0.9, ls=(0, (3, 2)) if group == "Ris dose 500mM" else "-")
    ax.axhline(1, color="#9a9a9a", lw=0.6, ls=(0, (2, 2)))
    ax.set_xlim(3000, 0)
    ax.set_ylim(0, 8.5)
    ax.set_xlabel("distance from the poly(A) site (nt)")
    ax.set_ylabel("coverage / gene mean\n(genes >= 1 kb, RP excluded)")
    handles = [Line2D([], [], color=c, lw=1.2) for c in list(groups.values()) + ["#7A4FA3"]]
    ax.legend(handles=handles, labels=list(groups) + ["RAW 264.7 (mouse)"], loc="upper left", frameon=False, fontsize=5.5, handlelength=1.2)
    save(fig, out_dir, name)


def heatmap(out_dir, name="p4i2_coverage_heatmaps_all_reads_1kb_noRP"):
    rows = read_tsv(find_table("p4i2_coverage_heatmaps_all_reads_1kb_noRP.tsv"))
    bins = [c for c in rows[0] if c.startswith("bin")]
    genes = sorted({r["gene"] for r in rows}, key=lambda g: float(next(r["transcript_length_nt"] for r in rows if r["gene"] == g)))
    fig, axes = plt.subplots(1, len(METHODS), figsize=(1.6 * len(METHODS), 3.0), sharey=True)
    for ax, method in zip(axes, METHODS):
        by_gene = {r["gene"]: [float(r[b]) for b in bins] for r in rows if r["method"] == method}
        ax.imshow(np.array([by_gene[g] for g in genes]), aspect="auto", cmap="viridis", vmin=0, vmax=1, interpolation="nearest")
        ax.set_title(method, fontsize=7, color=COLOR[method])
        ax.set_xticks([0, len(bins) - 1])
        ax.set_xticklabels(["5'", "3'"])
    axes[0].set_ylabel(f"{len(genes):,} genes (short to long)")
    axes[0].set_yticks([])
    save(fig, out_dir, name)


# ---------------------------------------------------------------- stacked QC figure
STACK_COLUMNS = [  # column key in the table, label, colour
    ("DRUG-seq", "DRUG-seq (n=24)", COLOR["DRUG-seq"]), ("prime-seq", "prime-seq (n=8)", COLOR["prime-seq"]),
    ("Hek control", "HEK control", COLOR["BOBseq"]), ("CADM1 1nM", "CADM1 1 nM", COLOR["BOBseq"]), ("CADM1 100nM", "CADM1 100 nM", COLOR["BOBseq"]),
    ("CHX dose 1ug mL", "CHX 1 ug/mL", COLOR["BOBseq"]), ("CHX dose 50ug mL", "CHX 50 ug/mL", COLOR["BOBseq"]), ("Ris dose 25mM", "Ris 25 nM", COLOR["BOBseq"]),
    ("Ris dose 500mM", "Ris 500 nM", "#8C8C8C"), ("RAW 264.7", "RAW control", "#7A4FA3"),
]
STACK_XPOS = [0, 1] + [2.35 + i for i in range(7)] + [9.7]
STACK_LIMITS = {  # y range per metric; metrics with "(log)" in the name use a log axis
    "reads into STAR (log)": (1e5, 3e7), "filtered reads after de-duplication (UMI, log)": (3e4, 3e6), "aligned bases (Gb, log)": (0.003, 3),
    "genes at 250k filtered reads (no ribosomal-protein genes)": (0, 13000), "intronic reads (% of mapped)": (0, 30),
    "mitochondrial reads (% of mapped)": (0, 20), "insert length median (nt)": (0, 300), "5'-3' balance (2 x centroid)": (0, 2),
}


def swarm_offsets(y, log, ylim):
    """Replicate 1, 2, 3 left to right for triplicates; a compact symmetric swarm for larger columns."""
    k = len(y)
    if k <= 3:
        return np.array([-0.2, 0.0, 0.2])[:k]
    t = np.log10(y) if log else np.asarray(y, float)
    lo, hi = (np.log10(ylim[0]), np.log10(ylim[1])) if log else ylim
    yn = (t - lo) / (hi - lo)
    offsets = np.zeros(k)
    placed = []
    candidates = [0.0] + [sign * step * 0.06 for step in range(1, 7) for sign in (-1, 1)]
    median = np.nanmedian(yn)
    for idx in np.argsort(np.where(np.isnan(yn), np.inf, np.abs(yn - median)), kind="stable"):
        if np.isnan(yn[idx]):
            continue
        free = [c for c in candidates if all(((c - px) / 0.06) ** 2 + ((yn[idx] - py) / 0.06) ** 2 >= 1 for px, py in placed)]
        offsets[idx] = free[0] if free else max(candidates, key=lambda c: min(((c - px) / 0.06) ** 2 + ((yn[idx] - py) / 0.06) ** 2 for px, py in placed))
        placed.append((offsets[idx], yn[idx]))
    missing = np.flatnonzero(np.isnan(yn))
    offsets[missing] = (np.arange(len(missing)) - (len(missing) - 1) / 2) * 0.12
    return offsets


def stacked_qc(out_dir, name="stacked_qc_all_samples"):
    rows = read_tsv(find_table("stacked_qc_all_samples.tsv"))
    metrics = list(dict.fromkeys(r["metric"] for r in rows))
    fig, axes = plt.subplots(len(metrics), 1, figsize=(4.7, 0.7 * len(metrics) + 1.7), sharex=True, gridspec_kw={"hspace": 0.3})
    for ax, metric in zip(axes, metrics):
        log = "(log)" in metric
        ylim = STACK_LIMITS.get(metric, (0, 100))
        for x0, (key, _, colour) in zip(STACK_XPOS, STACK_COLUMNS):
            cells = [r for r in rows if r["metric"] == metric and r["column"] == key]
            y = np.array([float(r["value"]) if r["value"] != "" else np.nan for r in cells])
            if metric.startswith("insert length") and np.isnan(y).all():
                ax.text(x0, np.mean(ylim), "single-end", ha="center", va="center", fontsize=5.5, color="#8a8a8a", rotation=90)
                continue
            x = x0 + swarm_offsets(y, log, ylim)
            ok = ~np.isnan(y)
            ax.scatter(x[ok], y[ok], s=4.5 if len(y) > 3 else 8, color=colour, lw=0, alpha=0.85, zorder=3)
            if ok.any():
                mean = float(np.exp(np.log(y[ok]).mean())) if log else float(y[ok].mean())
                ax.hlines(mean, x0 - 0.38, x0 + 0.38, color="#231F20", lw=0.8, zorder=4)
            if (~ok).any() and "250k" in metric:
                ax.scatter(x[~ok], np.full((~ok).sum(), ylim[0]), marker="x", s=9, color="#9a9a9a", lw=0.6, zorder=5, clip_on=False)
        if log:
            ax.set_yscale("log")
        ax.set_ylim(*ylim)
        ax.set_xlim(-0.7, STACK_XPOS[-1] + 0.7)
        ax.set_ylabel(metric, fontsize=6.3, rotation=0, ha="right", va="center", labelpad=6)
        ax.tick_params(axis="x", bottom=False)
        if metric.startswith("5'-3'"):
            ax.axhline(1, color="#9a9a9a", lw=0.5, ls=(0, (3, 2)), zorder=0)
        for xs in (1.68, 9.03):
            ax.axvline(xs, color="#D0D0D0", lw=0.5, zorder=0)
    axes[-1].set_xticks(STACK_XPOS)
    axes[-1].set_xticklabels([c[1] for c in STACK_COLUMNS], fontsize=6.5, rotation=90)
    for label, column in zip(axes[-1].get_xticklabels(), STACK_COLUMNS):
        label.set_color(column[2])
    fig.text(0.5, 0.0, "one dot per sample, line = mean (geometric mean on log rows); grey cross = sample with fewer than 250k filtered reads, not subsampled;\n"
             "Ris 500 nM = input failures; mouse rRNA includes the reads the aligner places on the human rDNA copies", ha="center", va="top", fontsize=5.5, color="#555")
    save(fig, out_dir, name)


# ---------------------------------------------------------------- gene coverage tracks
def gene_tracks(gene, out_dir):
    models = {r["gene"]: r for r in read_tsv(find_table("gene_models.tsv"))}
    model = models[gene]
    depth = read_tsv(find_table(f"{gene}_depth.tsv"))
    pos = np.array([int(r["pos_hg38"]) for r in depth])
    tracks = {m: np.array([float(r[m]) for r in depth]) for m in METHODS}
    starts = [int(x) for x in model["exon_starts_1based"].split(",")]
    ends = [int(x) for x in model["exon_ends"].split(",")]
    reverse = model["strand"] == "-"
    gene_colour = "#2F6DB5"

    def draw(x, data, exon_x, ticks, tick_labels, xlim, name, joined):
        fig, axes = plt.subplots(len(METHODS) + 1, 1, figsize=(6.8, 0.95 * len(METHODS) + 0.75),
                                 gridspec_kw={"height_ratios": [1] * len(METHODS) + [0.62], "hspace": 0.12}, sharex=True)
        for ax, method in zip(axes, METHODS):
            y = data[method]
            peak = int(y.max())
            ax.fill_between(x, 0, y, step="mid", color=COLOR[method], lw=0)
            ax.set_ylim(0, max(peak, 1) * 1.3)
            ax.set_yticks([peak])
            ax.set_yticklabels([f"{peak:,}"], fontsize=9)
            for side in ("top", "right", "bottom"):
                ax.spines[side].set_visible(False)
            ax.tick_params(axis="x", bottom=False, labelbottom=False)
            ax.text(0.015, 0.97, method, color=COLOR[method], fontsize=9.5, transform=ax.transAxes, va="top", ha="left",
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.2))
        ax = axes[-1]
        ax.set_ylim(-1.35, 1)
        ax.plot([exon_x[0][0], exon_x[-1][1]], [0, 0], color=gene_colour, lw=0.9)
        for a, b in exon_x:
            ax.add_patch(plt.Rectangle((min(a, b), -0.42), abs(b - a), 0.84, color=gene_colour, lw=0))
        if joined:
            for a, _ in exon_x[1:]:
                ax.plot([a, a], [-0.42, 0.42], color="white", lw=1.1)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_position(("data", -1.0))
        ax.set_yticks([])
        ax.set_xticks(ticks)
        ax.set_xticklabels(tick_labels, fontsize=7)
        ax.text(-0.012, 0.69, gene, color=gene_colour, fontsize=9.5, transform=ax.transAxes, ha="right", va="center")
        ax.text(1.012, 0.69, "3p", color=gene_colour, fontsize=9.5, transform=ax.transAxes, ha="left", va="center")
        ax.text(-0.012, 0.13, "hg38", fontsize=8.5, transform=ax.transAxes, ha="right", va="center")
        axes[0].set_xlim(*xlim)
        save(fig, out_dir, name)

    span = pos[-1] - pos[0]
    pad = 0.03 * span
    step = next(v for v in (500, 1000, 2000, 5000, 10000, 20000) if span / v <= 6)
    ticks = list(range((pos[0] // step + 1) * step, pos[-1], step))
    xlim = (pos[-1] + pad, pos[0] - pad) if reverse else (pos[0] - pad, pos[-1] + pad)
    draw(pos, tracks, list(zip(starts, ends)), ticks, [f"{t:,}" for t in ticks], xlim, f"{gene}_with_introns", False)

    exons = list(zip(starts, ends))
    if reverse:
        exons = exons[::-1]
    index = {p: i for i, p in enumerate(pos)}
    pieces = {m: [] for m in METHODS}
    for a, b in exons:
        sel = slice(index[a], index[b] + 1)
        for m in METHODS:
            piece = tracks[m][sel]
            pieces[m].append(piece[::-1] if reverse else piece)
    joined = {m: np.concatenate(pieces[m]) for m in METHODS}
    lengths = [b - a + 1 for a, b in exons]
    offsets = np.concatenate([[0], np.cumsum(lengths)])
    exon_x = [(offsets[i], offsets[i + 1]) for i in range(len(exons))]
    total = offsets[-1]
    chosen = [0]
    for i in range(1, len(exons)):
        if offsets[i] - offsets[chosen[-1]] >= total / 7:
            chosen.append(i)
    labels = [f"{(exons[i][1] if reverse else exons[i][0]):,}" for i in chosen]
    draw(np.arange(total) + 0.5, joined, exon_x, [offsets[i] for i in chosen], labels, (-0.03 * total, total * 1.03), f"{gene}_exons_only", True)


# ---------------------------------------------------------------- panel registry
def all_panels(out_dir):
    plan = {
        "p1_molecules_vs_depth_native": lambda: depth_curve("p1_molecules_vs_depth_native.tsv", "depth_filtered_reads", "molecules_UCI", "filtered reads per sample (subsampled)", "molecules per sample (UMI)", DEPTH_TICKS, out_dir, "p1_molecules_vs_depth_native", legend_loc="upper left"),
        "p1_molecules_vs_depth_log_native": lambda: depth_curve("p1_molecules_vs_depth_native.tsv", "depth_filtered_reads", "molecules_UCI", "filtered reads per sample (subsampled)", "molecules per sample (UMI)", DEPTH_TICKS, out_dir, "p1_molecules_vs_depth_log_native", log_y=True),
        "p2_genes_vs_depth_native": lambda: depth_curve("p2_genes_vs_depth_native.tsv", "depth_filtered_reads", "genes", "filtered reads per sample (subsampled)", "genes per sample\n(ribosomal-protein genes excluded)", DEPTH_TICKS, out_dir, "p2_genes_vs_depth_native", legend_loc="lower right"),
        "p5_junctions_vs_depth_min3_native": lambda: depth_curve("p5_junctions_vs_depth_min3_native.tsv", "depth", "junctions_expected", "deduplicated molecules per sample\n(subsampled)", "unique junctions (>= 3 molecules)", MOLECULE_TICKS, out_dir, "p5_junctions_vs_depth_min3_native"),
        "p5_annotated_junctions_vs_depth_min3_native": lambda: depth_curve("p5_annotated_junctions_vs_depth_min3_native.tsv", "depth", "annotated_junctions_expected", "deduplicated molecules per sample\n(subsampled)", "unique annotated junctions\n(>= 3 molecules per junction)", MOLECULE_TICKS, out_dir, "p5_annotated_junctions_vs_depth_min3_native"),
        "p5_pooled_annotated_junctions_vs_depth_min3_native": lambda: pooled_curve("p5_pooled_annotated_junctions_vs_depth_min3_native.tsv", "annotated_junctions_expected", "unique annotated junctions\n(>= 3 molecules per junction, pooled)", out_dir, "p5_pooled_annotated_junctions_vs_depth_min3_native"),
        "p5_pooled_all_junctions_vs_depth_min3_native": lambda: pooled_curve("p5_pooled_all_junctions_vs_depth_min3_native.tsv", "all_junctions_expected", "unique junctions, annotated or not\n(>= 3 molecules per junction, pooled)", out_dir, "p5_pooled_all_junctions_vs_depth_min3_native"),
        "p0c_reads_after_dedup_native": lambda: strip("p0c_p0d_dedup_per_sample_native.tsv", "molecules_UCI", "filtered reads per sample\nafter de-duplication (molecules, UMI)", out_dir, "p0c_reads_after_dedup_native", log_y=True, yticks=([1e5, 2e5, 5e5, 1e6, 2e6], ["100k", "200k", "500k", "1M", "2M"])),
        "p0d2_duplicate_rate_at_250k_native": lambda: strip("p0c_p0d_dedup_per_sample_native.tsv", "duplicate_rate_pct_at_250k", "duplicate rate (%) at 250k filtered reads\n1 - molecules / reads, matched depth", out_dir, "p0d2_duplicate_rate_at_250k_native", ylim=(0, 100)),
        "p3b_intronic_per_sample_native": lambda: strip("p3b_intronic_per_sample_native.tsv", "intronic_pct_of_mapped", "intronic reads (% of mapped reads)\nper sample", out_dir, "p3b_intronic_per_sample_native", ylim=(0, None)),
        "p6b_picard_balance_native": lambda: strip("p6b_picard_balance_native.tsv", "balance_2x_centroid", "5' to 3' balance (2 x centroid)", out_dir, "p6b_picard_balance_native", ylim=(0, 2)),
        "p7e_aligned_bases_total_native": lambda: strip("p7e_p7f_aligned_bases_per_sample_native.tsv", "total_aligned_bases", "aligned bases per sample (Gb)\nas sequenced", out_dir, "p7e_aligned_bases_total_native", log_y=True, scale=1e-9),
        "p3_read_composition_native": lambda: composition(out_dir),
        "p0a_read_fate_native": lambda: read_fate(out_dir),
        "p4m2_coverage_vs_nt_from_polyA_mean_1kb_noRP": lambda: profile_p4m2("p4m2_coverage_vs_nt_from_polyA_mean_1kb_noRP.tsv", out_dir, "p4m2_coverage_vs_nt_from_polyA_mean_1kb_noRP"),
        "p4m2_coverage_vs_nt_from_polyA_mean_1kb_noRP_50nt": lambda: profile_p4m2("p4m2_coverage_vs_nt_from_polyA_mean_1kb_noRP_50nt.tsv", out_dir, "p4m2_coverage_vs_nt_from_polyA_mean_1kb_noRP_50nt", note="(all methods: one 50-nt read)"),
        "S14_coverage_vs_nt_from_polyA_all24_native": lambda: profile_all24(out_dir),
        "p4i2_coverage_heatmaps_all_reads_1kb_noRP": lambda: heatmap(out_dir),
        "stacked_qc_all_samples": lambda: stacked_qc(out_dir),
    }
    try:   # the gene coverage examples are not shipped with the repository; regenerate them with benchmark/scripts/gene_coverage_tracks.py
        genes = [r["gene"] for r in read_tsv(find_table("gene_models.tsv"))]
    except FileNotFoundError:
        genes = []
    for gene in genes:
        plan[f"gene_{gene}"] = (lambda g=gene: gene_tracks(g, out_dir))
    return plan


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default=os.path.join(HERE, "replotted"))
    parser.add_argument("--only", nargs="*", help="panel names (see --list)")
    parser.add_argument("--list", action="store_true", help="list the panel names and exit")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    plan = all_panels(args.out)
    if args.list:
        print("\n".join(plan))
        raise SystemExit
    for name, draw in plan.items():
        if args.only and name not in args.only:
            continue
        draw()

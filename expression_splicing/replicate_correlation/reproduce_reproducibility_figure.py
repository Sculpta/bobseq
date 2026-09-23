#!/usr/bin/env python3
"""Reproduce the within-method replicate-reproducibility figures (BOBseq vs prime-seq vs DRUG-seq).

Two panels of a chosen inter-sample metric on an mRNA gene set:
  P1  native sequencing depth, genes with >= 5 counts in both samples of a pair.
  P2  all libraries down-sampled to the smallest library (mRNA-assigned counts),
      metric over each pair's top-N most-expressed co-detected genes.

Metric (--metric):
  pearson   per-pair Pearson correlation of log2(CPM+1)   (higher = more reproducible)
  rmse      per-pair RMSE from the identity line y=x on log2(CPM+1)  (lower = more reproducible)

Reproducibility is computed within method: BOBseq only between replicates of the same
condition (pooled across conditions); prime-seq and DRUG-seq are each a single condition.
Significance is a two-sided Mann-Whitney U on the PER-SAMPLE mean metric (the ~independent
unit), because replicate pairs are not independent.

Inputs (see README.md and datasets.csv):
  --counts   gene x sample integer count matrix (featureCounts), CSV or CSV.GZ.
             First two columns: GENEID, SYMBOL; remaining columns one per sample.
             Sample columns are prefixed by method set, e.g.
               bobseq_pe_native__<condition>-<rep>   (BOBseq)
               primeseq_native__<sample>             (prime-seq)
               drugseq_native__<well>                (DRUG-seq)
  --annot    gene annotation: Ensembl GTF (.gtf/.gtf.gz) OR a TSV with columns
             gene_id, biotype, chromosome, gene_name.

Usage:
  python reproduce_reproducibility_figure.py --counts data/gene_counts.csv.gz \
      --annot data/Homo_sapiens.GRCh38.113.gtf.gz --metric pearson --out figure_pearson
  python reproduce_reproducibility_figure.py --counts data/gene_counts.csv.gz \
      --annot data/Homo_sapiens.GRCh38.113.gtf.gz --metric rmse    --out figure_rmse

Outputs: <out>.svg (editable text), <out>.pdf, <out>.png, <out>.tsv (summary + p-values).
"""
import argparse
import gzip
import itertools
import re
import copy
from collections import defaultdict
import numpy as np
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

# ----------------------------------------------------------------------------------------
# Configuration: sample selection (edit column prefixes here to match your count matrix)
# ----------------------------------------------------------------------------------------
BOB_PREFIX = "bobseq_pe_native__"
PRIME_PREFIX = "primeseq_native__"
DRUG_PREFIX = "drugseq_native__"

# BOBseq replicate pairs are formed WITHIN each condition (triplicates), then pooled.
# Ris-500 (QC-failed input) is intentionally excluded.
BOB_CONDITIONS = [
    ["Hek_control_1", "Hek_control_2", "Hek_control_3"],
    ["CADM1_1nM-1", "CADM1_1nM-2", "CADM1_1nM-3"],
    ["CADM1_100nM-1", "CADM1_100nM-2", "CADM1_100nM-3"],
    ["CHX_dose_1ug_mL-1", "CHX_dose_1ug_mL-2", "CHX_dose_1ug_mL-3"],
    ["CHX_dose_50ug_mL-1", "CHX_dose_50ug_mL-2", "CHX_dose_50ug_mL-3"],
    ["Ris_dose_25mM-1", "Ris_dose_25mM-2", "Ris_dose_25mM-3"],
]
CELL = {"BOBseq": "HEK", "prime-seq": "HEK293T", "drug-seq": "U-2 OS"}
BOB_C, PRIME_C, DRUG_C = "#2C9A4F", "#0292D7", "#D42E97"
P2_TOP_N = 4168          # top genes/pair in P2 (max with no dropout at the full-set floor; see README)
RIBO_PROT = re.compile(r"^(RPL|RPS|MRPL|MRPS)")

# per-metric display settings: (axis label, fixed panel y-limits or None -> auto, short label, direction)
METRICS = {
    "pearson": dict(ylabel="Pearson r", ylims=[(0.88, 1.00), (0.50, 0.70)], short="r",
                    note="Pearson correlation of log2(CPM+1); higher = more reproducible"),
    "rmse": dict(ylabel="RMSE from y = x", ylims=[None, None], short="RMSE",
                 note="RMSE from identity (y=x) on log2(CPM+1); lower = more reproducible"),
}


# ----------------------------------------------------------------------------------------
# I/O
# ----------------------------------------------------------------------------------------
def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def load_counts(path):
    with _open(path) as f:
        header = f.readline().rstrip("\n").split(",")
        cols = header[2:]
        gene_ids, mat = [], []
        for line in f:
            p = line.rstrip("\n").split(",")
            gene_ids.append(p[0])
            mat.append([int(float(x)) for x in p[2:]])
    M = np.array(mat, dtype=np.int64)
    return np.array(gene_ids), cols, {c: M[:, i] for i, c in enumerate(cols)}


def load_annotation(path):
    """Return dicts biotype/chrom/symbol keyed by gene_id, from a TSV or an Ensembl GTF."""
    bt, ch, nm = {}, {}, {}
    if str(path).endswith((".gtf", ".gtf.gz")):
        gid_re = re.compile(r'gene_id "([^"]+)"')
        bio_re = re.compile(r'gene_biotype "([^"]+)"')
        nam_re = re.compile(r'gene_name "([^"]+)"')
        with _open(path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) < 9 or fields[2] != "gene":
                    continue
                g = gid_re.search(fields[8])
                if not g:
                    continue
                gid = g.group(1)
                b = bio_re.search(fields[8]); n = nam_re.search(fields[8])
                bt[gid] = b.group(1) if b else ""
                ch[gid] = fields[0]
                nm[gid] = n.group(1) if n else ""
    else:  # TSV: gene_id, biotype, chromosome, gene_name
        with _open(path) as f:
            first = f.readline()
            if not first.lower().startswith("gene_id"):
                f.seek(0)
            for line in f:
                gid, biotype, chrom, sym = (line.rstrip("\n").split("\t") + ["", "", "", ""])[:4]
                bt[gid] = biotype; ch[gid] = chrom; nm[gid] = sym
    return bt, ch, nm


def mrna_mask(gene_ids, bt, ch, nm):
    """protein-coding, minus ribosomal-protein genes, minus mitochondrial (chr MT)."""
    def keep(g):
        if bt.get(g) != "protein_coding":
            return False
        if ch.get(g) in ("MT", "chrMT", "HUMAN_MT"):
            return False
        s = nm.get(g, "")
        if RIBO_PROT.match(s) and not s.startswith("RPS6K"):
            return False
        return True
    return np.array([keep(g) for g in gene_ids])


# ----------------------------------------------------------------------------------------
# Analysis
# ----------------------------------------------------------------------------------------
def downsample(vec, n, rng):
    v = np.asarray(vec).astype(int)
    if v.sum() <= n:
        return v.astype(float)
    idx = rng.choice(np.repeat(np.arange(len(v)), v), n, replace=False)
    return np.bincount(idx, minlength=len(v)).astype(float)


def cpm_log(cnt_dict, cols):
    cpm = {c: cnt_dict[c] / cnt_dict[c].sum() * 1e6 for c in cols}
    return cpm, {c: np.log2(cpm[c] + 1.0) for c in cols}


def reduce_metric(la, lb, metric):
    if metric == "pearson":
        return float(np.corrcoef(la, lb)[0, 1])
    d = la - lb
    return float(np.sqrt(np.mean(d * d)))


def per_sample_means(triples):
    d = defaultdict(list)
    for a, b, v in triples:
        d[a].append(v); d[b].append(v)
    return np.array([np.mean(v) for v in d.values()])


def pairwise_block(cols, value_fn, within=None):
    pairs = ([p for grp in within for p in itertools.combinations(grp, 2)]
             if within else list(itertools.combinations(cols, 2)))
    trip, ngenes = [], []
    for a, b in pairs:
        v, ng = value_fn(a, b)
        trip.append((a, b, v)); ngenes.append(ng)
    return (np.array([t[2] for t in trip]), np.array(ngenes), per_sample_means(trip))


def value_threshold(cnt, logc, thr, metric):
    def f(a, b):
        m = (cnt[a] >= thr) & (cnt[b] >= thr)
        return reduce_metric(logc[a][m], logc[b][m], metric), int(m.sum())
    return f


def value_topN(cpm, logc, N, metric):
    def f(a, b):
        top = np.argsort((cpm[a] + cpm[b]) * 0.5)[::-1][:N]
        return reduce_metric(logc[a][top], logc[b][top], metric), int(len(top))
    return f


BOB_CONDITIONS_RESOLVED = None


def build_groups(value_fn, bob_cols, prime_cols, drug_cols):
    out = []
    for meth, cols, within, col, lab in [
        ("BOBseq", bob_cols, BOB_CONDITIONS_RESOLVED, BOB_C, "BOBseq\n(within-cond.)"),
        ("prime-seq", prime_cols, None, PRIME_C, "prime-seq"),
        ("drug-seq", drug_cols, None, DRUG_C, "drug-seq\n(all 24)"),
    ]:
        v, ng, ps = pairwise_block(cols, value_fn, within=within)
        out.append(dict(method=meth, label=lab, color=col, cols=cols, v=v, ngenes=ng, per_sample=ps))
    return out


# ----------------------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------------------
def _stars(p):
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"


def _fmt_reads(x):
    return f"{x/1e6:.2f}M" if x >= 1e6 else (f"{x/1e3:.0f}K" if x >= 1e3 else f"{x:.0f}")


POS = [0.0, 0.55, 1.10]


def draw_panel(ax, groups, title, ylim, ylabel, short, lib, show_labels, stats_rows):
    vals = [g["v"] for g in groups]
    if ylim is None:  # auto (rmse): headroom at top for the median row
        lo = min(v.min() for v in vals); hi = max(v.max() for v in vals); rng = hi - lo
        ylim = (lo - rng * 0.08, hi + rng * 0.30)
    T = ax.get_xaxis_transform()
    bp = ax.boxplot(vals, positions=POS, widths=0.40, showfliers=False, patch_artist=True,
                    medianprops=dict(color="#111", lw=1.6), whiskerprops=dict(color="#888", lw=1.0),
                    capprops=dict(color="#888", lw=1.0), boxprops=dict(lw=0))
    for patch, g in zip(bp["boxes"], groups):
        patch.set_facecolor(g["color"]); patch.set_alpha(0.18)
    for k, g in enumerate(groups):
        jit = (np.random.default_rng(4).random(len(g["v"])) - 0.5) * 0.18
        ax.scatter(np.full(len(g["v"]), POS[k]) + jit, g["v"], s=14, color=g["color"], alpha=0.35, lw=0, zorder=3)
        ax.text(POS[k], 0.95, f"{np.median(g['v']):.3f}", transform=T, ha="center", va="center",
                fontsize=8.0, color=g["color"], fontweight="bold")
    ax.set_ylim(*ylim); ax.set_xlim(-0.30, 1.40)
    ax.set_xticks(POS); ax.set_xticklabels([g["label"] for g in groups], fontsize=7.6)
    ax.set_ylabel(ylabel)
    for i, (a, b) in enumerate(itertools.combinations(range(len(groups)), 2)):
        p = float(mannwhitneyu(groups[a]["per_sample"], groups[b]["per_sample"], alternative="two-sided").pvalue)
        stats_rows.append((title.split("·")[0].strip(), groups[a]["method"], groups[b]["method"], p))
        y = 1.03 + i * 0.07
        ax.plot([POS[a], POS[a], POS[b], POS[b]], [y - 0.017, y, y, y - 0.017],
                transform=T, color="#444", lw=0.9, clip_on=False, zorder=6)
        ax.text((POS[a] + POS[b]) / 2, y + 0.006, _stars(p), transform=T, ha="center", va="bottom",
                fontsize=7.6, color="#333", clip_on=False)
    ax.text(0.0, 1.27, title, transform=ax.transAxes, ha="left", va="bottom", fontsize=7.8)
    rows = ["cell line", "samples (n)", "pairs (n)", "reads range", "genes/comp.", f"{short} IQR"]
    ry = [-0.15, -0.195, -0.24, -0.285, -0.33, -0.375]
    if show_labels:
        for lbl, y in zip(rows, ry):
            ax.text(-0.02, y, lbl, transform=ax.transAxes, ha="right", va="center", fontsize=5.6, color="#777", clip_on=False)
    for k, g in enumerate(groups):
        rmin = min(lib[c] for c in g["cols"]); rmax = max(lib[c] for c in g["cols"])
        iqr = float(np.subtract(*np.percentile(g["v"], [75, 25])))
        vals_col = [CELL[g["method"]], f"{len(g['cols'])}", f"{len(g['v'])}",
                    f"{_fmt_reads(rmin)}–{_fmt_reads(rmax)}", f"{g['ngenes'].mean():,.0f}", f"{iqr:.3f}"]
        for v, y in zip(vals_col, ry):
            ax.text(POS[k], y, v, transform=T, ha="center", va="center", fontsize=5.5, color="#333", clip_on=False)


# --- make matplotlib SVG text editable in Illustrator (real <text>, inlined <use>) ------
SVG_NS, XLINK = "http://www.w3.org/2000/svg", "http://www.w3.org/1999/xlink"
def _style_dict(s):
    d = {}
    for part in (s or "").split(";"):
        if ":" in part:
            k, v = part.split(":", 1); d[k.strip()] = v.strip()
    return d
def make_svg_editable(path):
    ET.register_namespace("", SVG_NS); ET.register_namespace("xlink", XLINK)
    tree = ET.parse(path); root = tree.getroot()
    defs = {e.get("id"): e for e in root.iter() if e.get("id")}
    parents = {c: p for p in root.iter() for c in p}
    for use in list(root.iter(f"{{{SVG_NS}}}use")):
        href = use.get(f"{{{XLINK}}}href") or use.get("href")
        if not href or not href.startswith("#") or href[1:] not in defs:
            continue
        new = copy.deepcopy(defs[href[1:]]); new.attrib.pop("id", None)
        tr = f"translate({use.get('x', '0')}, {use.get('y', '0')})"
        if new.get("transform"):
            tr += " " + new.get("transform")
        new.set("transform", tr)
        m = _style_dict(new.get("style")); m.update(_style_dict(use.get("style")))
        if m:
            new.set("style", ";".join(f"{k}:{v}" for k, v in m.items()))
        p = parents[use]; p.insert(list(p).index(use), new); p.remove(use)
    tree.write(path, xml_declaration=True, encoding="utf-8")


# ----------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True, help="gene x sample count matrix (CSV/CSV.GZ)")
    ap.add_argument("--annot", required=True, help="Ensembl GTF (.gtf/.gtf.gz) or TSV (gene_id,biotype,chrom,gene_name)")
    ap.add_argument("--metric", choices=["pearson", "rmse"], default="pearson", help="reproducibility metric")
    ap.add_argument("--out", default=None, help="output file stem (default: figure_<metric>)")
    args = ap.parse_args()
    out = args.out or f"figure_{args.metric}"
    M = METRICS[args.metric]

    gene_ids, cols, counts = load_counts(args.counts)
    bt, ch, nm = load_annotation(args.annot)
    mask = mrna_mask(gene_ids, bt, ch, nm)
    print(f"genes: {len(gene_ids):,} | mRNA set: {mask.sum():,} | metric: {args.metric}")

    counts = {c: counts[c][mask].astype(np.int64) for c in cols}
    lib = {c: int(counts[c].sum()) for c in cols}

    global BOB_CONDITIONS_RESOLVED
    BOB_CONDITIONS_RESOLVED = [[BOB_PREFIX + s for s in grp] for grp in BOB_CONDITIONS]
    bob_cols = [c for grp in BOB_CONDITIONS_RESOLVED for c in grp]
    prime_cols = sorted(c for c in cols if c.startswith(PRIME_PREFIX))
    drug_cols = sorted(c for c in cols if c.startswith(DRUG_PREFIX))
    missing = [c for c in bob_cols if c not in counts]
    if missing:
        raise SystemExit(f"missing BOBseq columns in matrix: {missing}")
    print(f"BOBseq {len(bob_cols)} | prime-seq {len(prime_cols)} | drug-seq {len(drug_cols)}")

    alls = bob_cols + prime_cols + drug_cols
    # P1: native depth, gene >= 5 in both
    _, log1 = cpm_log({c: counts[c].astype(float) for c in alls}, alls)
    P1 = build_groups(value_threshold({c: counts[c].astype(float) for c in cols}, log1, 5, args.metric),
                      bob_cols, prime_cols, drug_cols)
    # P2: depth-match all to the full-set minimum; top-N genes/pair
    floor = min(lib[c] for c in alls)
    rng = np.random.default_rng(0)
    cnt2 = {c: downsample(counts[c], floor, rng) for c in alls}
    cpm2, log2 = cpm_log(cnt2, alls)
    P2 = build_groups(value_topN(cpm2, log2, P2_TOP_N, args.metric), bob_cols, prime_cols, drug_cols)
    print(f"P2 depth floor (full-set min, mRNA-assigned counts) = {floor:,}")

    panels = [("P1 · native depth · gene ≥5 in both", P1, M["ylims"][0]),
              (f"P2 · all 24 drug · depth-matched to {floor:,} · top {P2_TOP_N:,} genes/pair", P2, M["ylims"][1])]

    plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "font.size": 9, "svg.fonttype": "none",
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 7.0))
    fig.subplots_adjust(left=0.17, right=0.985, bottom=0.30, top=0.80, wspace=0.42)
    stats_rows = []
    for j, (ax, (title, groups, ylim)) in enumerate(zip(axes, panels)):
        draw_panel(ax, groups, title, ylim, M["ylabel"], M["short"], lib, show_labels=(j == 0), stats_rows=stats_rows)
    fig.text(0.01, 0.02,
             f"Within-method inter-sample variability = {M['note']}. "
             "mRNA set (protein-coding − ribosomal-protein − mitochondrial); BOBseq = within-condition pairs pooled.\n"
             f"Brackets = two-sided Mann-Whitney U on per-sample mean {M['short']} (~independent unit): "
             "** p<0.01  **** p<1e-4. · reads range = native mRNA-assigned per-sample depth.",
             ha="left", va="bottom", fontsize=6.3, color="#333")
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{out}.{ext}", bbox_inches="tight", pad_inches=0.06, dpi=600)
    plt.close(fig)
    make_svg_editable(f"{out}.svg")

    with open(f"{out}.tsv", "w") as fh:
        fh.write(f"# metric={args.metric}\n")
        fh.write(f"panel\tmethod\tcell_line\tn_samples\tn_pairs\treads_native_min\treads_native_max\t"
                 f"genes_per_comparison\tmedian_{M['short']}\t{M['short']}_IQR\n")
        for title, groups, _ in panels:
            tag = title.split("·")[0].strip()
            for g in groups:
                rmin = min(lib[c] for c in g["cols"]); rmax = max(lib[c] for c in g["cols"])
                iqr = float(np.subtract(*np.percentile(g["v"], [75, 25])))
                fh.write(f"{tag}\t{g['method']}\t{CELL[g['method']]}\t{len(g['cols'])}\t{len(g['v'])}\t"
                         f"{rmin}\t{rmax}\t{g['ngenes'].mean():.0f}\t{np.median(g['v']):.4f}\t{iqr:.4f}\n")
        fh.write(f"\n# pairwise Mann-Whitney U (two-sided) on per-sample mean {M['short']}\npanel\tgroup_A\tgroup_B\tp_value\n")
        for tag, a, b, p in stats_rows:
            fh.write(f"{tag}\t{a}\t{b}\t{p:.3e}\n")
    print(f"wrote {out}.svg/.pdf/.png/.tsv")


if __name__ == "__main__":
    main()

"""screen_splicing: differential splicing of the BOBseq 24-plex (bobseq_pe_native), five contrasts.

Shared facts for every stage: locations, the wells, the contrasts, the condition palette and order, the count loader and the
small statistics helpers. Nothing here hard-codes a number from a previous run."""
import logging
import csv, os, sys
import numpy as np
from scipy.stats import beta

HERE = os.path.dirname(os.path.abspath(__file__)); ANALYSIS = os.path.dirname(HERE)
sys.path.insert(0, os.path.dirname(ANALYSIS))
from locations import table, SAMPLES
COUNTS_MOL = "S4"       # S4_junction_counts_dedup.csv.gz: junction counts in UMI-deduplicated molecules
COUNTS_FRAG = "S6"      # S6_junction_counts_fragments.csv.gz: junction counts per fragment (no deduplication), the differential-splicing unit
RESULTS = os.path.join(ANALYSIS, "results"); FIGURES = os.path.join(ANALYSIS, "figures")
SET = "bobseq_pe_native"

# --- conditions --------------------------------------------------------------------------------------
ORDER = ["Hek control", "Ris dose 25mM", "Ris dose 500mM", "CHX dose 1ug/mL", "CHX dose 50ug/mL", "CADM1 1nM", "CADM1 100nM"]
BAR_ORDER = [c for c in ORDER if c != "Ris dose 500mM"]            # what the bar charts show: every condition but the dead Ris 500 mM wells
COL = {"Hek control": "#444444", "Ris dose 25mM": "#1f77b4", "Ris dose 500mM": "#9ecae1", "CHX dose 1ug/mL": "#fdae6b",
       "CHX dose 50ug/mL": "#d94801", "CADM1 1nM": "#a1d99b", "CADM1 100nM": "#31a354"}
LAB = {"Hek control": "control", "Ris dose 25mM": "Ris 25 mM", "Ris dose 500mM": "Ris 500 mM", "CHX dose 1ug/mL": "CHX 1 µg/mL",
       "CHX dose 50ug/mL": "CHX 50 µg/mL", "CADM1 1nM": "CADM1 1 nM", "CADM1 100nM": "CADM1 100 nM"}

# --- contrasts --------------------------------------------------------------------------------------------
# kind "pairwise": Welch t-test on per-well PSI, treated vs control.  kind "dose": OLS slope of per-well PSI on dose
# rank 0/1/2 (t-test on the slope, 7 residual df); effect size = dPSI(top dose - control) from pooled counts in both.
CONTRASTS = {
    "chx_dose":   dict(kind="dose",     groups={"Hek control": 0, "CHX dose 1ug/mL": 1, "CHX dose 50ug/mL": 2},
                       title="CHX dose trend (control / 1 / 50 µg/mL)", color=COL["CHX dose 50ug/mL"]),
    "chx_high":   dict(kind="pairwise", groups={"Hek control": 0, "CHX dose 50ug/mL": 1},
                       title="CHX 50 µg/mL vs control", color=COL["CHX dose 50ug/mL"]),
    "ris_low":    dict(kind="pairwise", groups={"Hek control": 0, "Ris dose 25mM": 1},
                       title="Ris 25 mM vs control", color=COL["Ris dose 25mM"]),
    "cadm1_dose": dict(kind="dose",     groups={"Hek control": 0, "CADM1 1nM": 1, "CADM1 100nM": 2},
                       title="CADM1-ASO dose trend (control / 1 / 100 nM)", color=COL["CADM1 100nM"]),
    "cadm1_high": dict(kind="pairwise", groups={"Hek control": 0, "CADM1 100nM": 1},
                       title="CADM1-ASO 100 nM vs control", color=COL["CADM1 100nM"]),
}
# prefilter (condition-blind) and hit rule
MIN_MEAN_TOTAL, MIN_VAR = 100, 0.005        # LSV side >= 100 fragments on average and > 0 in every well; PSI variance >= 0.005
MAX_PADJ, HIT_DPSI = 0.10, 0.10             # hit = padj < 0.1 and dPSI >= 0.1, positive direction only
N_TOP = 10

# --- wells ---------------------------------------------------------------------------------------------
# Controls = Hek_control_1/2/3: Hek_control_2 is in_benchmark = no (T-depleted UMI oligo) but is a normal HEK293T well by
# every other measure and is used as the third control. The Ris 500 mM wells (input failure) are never used.
DS_EXTRA = {"Hek_control_2"}

def wells():
    """The 21 human wells of the 24-plex, in condition order, with the sample-sheet columns."""
    with open(SAMPLES) as fh:
        w = [r for r in csv.DictReader(fh) if r["set"] == SET]
    return sorted(w, key=lambda r: (ORDER.index(r["condition"]), r["sample_name"]))

def usable_wells():
    return [w for w in wells() if w["in_benchmark"] == "yes" or w["sample_name"] in DS_EXTRA]

def contrast_wells(tag):
    return [w for w in usable_wells() if w["condition"] in CONTRASTS[tag]["groups"]]

# --- counts ---------------------------------------------------------------------------------------------
def load_counts(which, accs):
    """Junction counts of the supplementary table `which` (COUNTS_FRAG or COUNTS_MOL) for the wells `accs`:
    (row ids, K rows x wells, T = per-row LSV total). An empty cell of the table is a count of 0."""
    import gzip
    with gzip.open(table(which), "rt") as fh:
        rd = csv.reader(fh); head = next(rd); idx = [head.index(a) for a in accs]; rn, K = [], []
        for row in rd:
            rn.append(row[0]); K.append([float(row[i]) if row[i] else 0.0 for i in idx])
    K = np.array(K); lsv = np.array([r.split("|")[0] for r in rn]); u, inv = np.unique(lsv, return_inverse=True)
    T = np.zeros((len(u), len(accs))); np.add.at(T, inv, K)
    return rn, K, T[inv]


def junction_ids(which):
    """The row ids of a junction table, in table order."""
    import gzip
    with gzip.open(table(which), "rt") as fh:
        next(fh); return [line.split(",", 1)[0] for line in fh]

def parse_row(rowname):
    """'chr19:CLASRP:t:45068419|19:45068054-45068419:+' -> dict(lsv, gene, side, anchor, junction, chrom, start, end, strand)."""
    lsv, j = rowname.split("|"); ch, gene, side, anchor = lsv.split(":"); c, coords, strand = j.split(":"); s, e = map(int, coords.split("-"))
    return dict(lsv=lsv, gene=gene, side=side, anchor=int(anchor), junction=j, chrom=c, start=s, end=e, strand=strand)

def ucsc_chrom(c):
    return "chrM" if c == "MT" else (c if c.startswith("chr") else f"chr{c}")

# --- statistics ------------------------------------------------------------------------------------------
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); r = np.empty(n)
    r[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(r, 1)

def cp_ci(k, n, alpha=0.05):
    """Clopper-Pearson interval for k of n; (nan, nan) at n = 0."""
    if n == 0: return (float("nan"), float("nan"))
    return (0.0 if k == 0 else beta.ppf(alpha / 2, k, n - k + 1), 1.0 if k == n else beta.ppf(1 - alpha / 2, k + 1, n - k))

# --- I/O -------------------------------------------------------------------------------------------------
def write_tsv(path, rows, cols):
    with open(path + ".tmp", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    os.replace(path + ".tmp", path)

def read_tsv(name):
    with open(os.path.join(RESULTS, name)) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

class Stats:
    """results/screen_splicing_stats.txt: the numbers quoted in README.md. One section per stage; a re-run replaces its
    own section and leaves the others untouched."""
    PATH = os.path.join(RESULTS, "screen_splicing_stats.txt")
    HEAD = "# screen_splicing -- the numbers quoted in README.md\n"

    def __init__(self, section): self.section, self.lines = section, []
    def add(self, line): self.lines.append(line)
    def write(self):
        sections, order, cur = {}, [], None
        if os.path.exists(self.PATH):
            for line in open(self.PATH).read().splitlines():
                if line.startswith("## "): cur = line[3:]; order.append(cur); sections[cur] = []
                elif cur is not None and line.strip(): sections[cur].append(line)
        if self.section not in sections: order.append(self.section)
        sections[self.section] = self.lines
        with open(self.PATH, "w") as fh:
            fh.write(self.HEAD)
            for name in order: fh.write(f"\n## {name}\n" + "\n".join(sections[name]) + "\n")

# --- figures -----------------------------------------------------------------------------------------------
def style():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib import font_manager as _fm       # Arial metrics on a box without Arial: register its
    import dataclasses as _dc                        # metric-identical twin Liberation Sans under the name Arial,
    for _p in (q for q in _fm.findSystemFonts() if "LiberationSans-" in q):   # so the layout is Arial-exact and the
        _fm.fontManager.ttflist.append(_dc.replace(_fm.ttfFontProperty(_fm.get_font(_p)), name="Arial"))   # file says Arial
    plt.rcParams.update({"font.family": "Arial",
                         "svg.fonttype": "none", "pdf.fonttype": 42, "savefig.dpi": 400, "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titlesize": 9, "axes.labelsize": 9, "legend.frameon": False, "savefig.facecolor": "white"})
    return plt

def save(fig, name, dpi=150):
    """figures/<name>.svg"""
    os.makedirs(FIGURES, exist_ok=True)
    fig.savefig(os.path.join(FIGURES, f"{name}.svg"))
    return os.path.join(FIGURES, f"{name}.svg")

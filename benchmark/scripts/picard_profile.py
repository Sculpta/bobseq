#!/usr/bin/env python3
"""Picard gene-body coverage panels for the preprint (no fragment-length rule).
Per benchmark sample (per-sample human BAMs in benchmark_uniform/per_sample_bams: primary mapped records of the benchmark alignments,
raw reads, native BOBseq = both mates): Picard CollectRnaSeqMetrics (Picard 3.5.0), STRAND_SPECIFICITY NONE,
MINIMUM_LENGTH 1000, refFlat = protein-coding Ensembl 113 transcripts with ribosomal-protein and mitochondrial genes removed
(bare chromosome names, matching these BAMs). Per set:
  p6_picard_profile_<tag>        mean over samples of the per-sample 101-bin normalised coverage, 95% t-interval ribbon
  p6b_picard_balance_<tag>       5'->3' balance = 2 x coverage centroid per sample (1 = even), points + median; _box = with the IQR box
  p6_picard_numbers_<tag>.json   per-sample balances, mean profiles, sample lists, settings
Metrics are cached in <fig set>/coverage_architecture/picard_metrics/<tag>/ (re-run only for missing files).
BM_COVSET=native|50nt selects the set; tag = native | matched."""

from settings import CODE, WORK, COVERAGE, GTF, REFFLAT, RUN_JSON
import os, re, sys, json, subprocess, collections, numpy as np, matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

P = WORK
PB = f'{P}/benchmark_uniform/per_sample_bams'
C = CODE
SET = os.environ.get('BM_COVSET', 'native')
TAG = 'native' if SET == 'native' else 'matched'
FIG = 'preprint_figures_native' if SET == 'native' else f'preprint_figures_{SET}'
OUT = f'{P}/{FIG}/main_figure_panels'
MET = f'{P}/{FIG}/coverage_architecture/picard_metrics/{TAG}'
os.makedirs(OUT, exist_ok=True)
os.makedirs(MET, exist_ok=True)
from palette import COL, METH

REFFLAT = f'{COVERAGE}/picard_refflat_pc_noRP_noMT_bare.txt'
PICARD = 'picard'
NPAR = int(os.environ.get('BM_PICARD_PAR', 5))
RIBO = re.compile(r'^(RP[LS]|RPLP|MRP[LS])[^K]*$', re.I)
# ---- refFlat (built once) ----
if not os.path.exists(REFFLAT):
    attr = re.compile(r'(\S+) "([^"]*)"')
    tx = {}
    ex = collections.defaultdict(list)
    cds = collections.defaultdict(list)
    for l in open(GTF):
        if l[0] == '#' or not l.startswith('HUMAN_'):
            continue
        f = l.rstrip('\n').split('\t')
        if f[2] not in ('exon', 'CDS'):
            continue
        a = dict(attr.findall(f[8]))
        if a.get('transcript_biotype') != 'protein_coding':
            continue
        chrom = f[0][6:]
        if chrom == 'MT':
            continue
        g = a.get('gene_name', a.get('gene_id', ''))
        t = a['transcript_id']
        if RIBO.match(g):
            continue
        tx[t] = (g, chrom, f[6])
        (ex if f[2] == 'exon' else cds)[t].append((int(f[3]) - 1, int(f[4])))
    with open(REFFLAT, 'w') as o:
        for t, (g, chrom, strand) in tx.items():
            e = sorted(ex[t])
            c = sorted(cds[t])
            cs, ce = (c[0][0], c[-1][1]) if c else (e[-1][1], e[-1][1])
            o.write(
                f'{g}\t{t}\t{chrom}\t{strand}\t{e[0][0]}\t{e[-1][1]}\t{cs}\t{ce}\t{len(e)}\t{",".join(str(s) for s, _ in e)},\t{",".join(str(x) for _, x in e)},\n'
            )
    print('refFlat written:', sum(1 for _ in open(REFFLAT)), 'transcripts')
# ---- sample list: competitors from per_sample_bams/metadata.tsv (in_benchmark), BOBseq = the benchmark wells (wells_keep minus units_excluded) ----
SETNAME = 'native' if SET == 'native' else 'matched_50nt'
samples = []
for l in open(f'{PB}/metadata.tsv'):
    f = l.rstrip('\n').split('\t')
    if (
        f[0] in ('drugseq', 'brbseq', 'primeseq')
        and {'drugseq': 'DRUG-seq', 'brbseq': 'BRB-seq', 'primeseq': 'prime-seq'}[f[0]] in METH
        and f[1] == SETNAME
        and f[8] == 'yes'
    ):
        samples.append(
            (
                {'drugseq': 'DRUG-seq', 'brbseq': 'BRB-seq', 'primeseq': 'prime-seq'}[f[0]],
                f[2],
                f[3],
            )
        )
RUN = json.load(open(RUN_JSON))
ARM = f'{P}/benchmark_uniform/bob57_24plex_pe_native'
excl = {l.split('\t')[0] for i, l in enumerate(open(f'{ARM}/units_excluded.tsv')) if i > 0}
bdir = 'bobseq_pe_native' if SET == 'native' else 'bobseq_50nt'
for w in (l.strip() for l in open(f'{ARM}/wells_keep.txt') if l.strip()):
    if w in excl or RUN['tso_species_map'][w] != 'human':
        continue
    n = RUN['bobcode_labels'][w].replace(' ', '_').replace('/', '_')
    samples.append(('BOBseq', n, f'{bdir}/{n}.bam'))
for m, n, b in samples:
    assert os.path.exists(f'{PB}/{b}'), f'missing BAM {b}'
print(TAG, 'samples:', collections.Counter(m for m, _, _ in samples))


# ---- Picard, cached ----
def run_one(item):
    m, n, b = item
    out = f'{MET}/{n}.RNA_Metrics.txt'
    if os.path.exists(out) and '## HISTOGRAM' in open(out).read():
        return n, 'cached'
    cmd = [
        PICARD,
        'CollectRnaSeqMetrics',
        '-I',
        f'{PB}/{b}',
        '-O',
        f'{MET}/{n}.RNA_Metrics.txt',
        '--REF_FLAT',
        REFFLAT,
        '--STRAND_SPECIFICITY',
        'NONE',
        '--MINIMUM_LENGTH',
        '1000',
        '--VALIDATION_STRINGENCY',
        'LENIENT',
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    open(f'{MET}/{n}.log', 'w').write(r.stdout + r.stderr)
    return n, ('done' if r.returncode == 0 else 'FAILED')


from concurrent.futures import ThreadPoolExecutor

with ThreadPoolExecutor(NPAR) as ex_:
    for n, st in ex_.map(run_one, samples):
        print(' ', n, st, flush=True)


def parse(path):
    rows = [l.rstrip('\n').split('\t') for l in open(path)]
    hist = []
    inh = False
    metr = {}
    for i, r in enumerate(rows):
        if r and r[0] == 'PF_BASES':
            metr = dict(zip(r, rows[i + 1]))
        if r and r[0].startswith('## HISTOGRAM'):
            inh = True
            continue
        if inh and len(r) >= 2 and r[0].isdigit():
            hist.append(float(r[1]))
    return np.array(hist), metr


prof = collections.defaultdict(list)
bal = collections.defaultdict(list)
names = collections.defaultdict(list)
p = np.arange(101) / 100
for m, n, b in samples:
    h, metr = parse(f'{MET}/{n}.RNA_Metrics.txt')
    assert len(h) == 101 and h.sum() > 0, n
    prof[m].append(h)
    bal[m].append(float(2 * (p * h).sum() / h.sum()))
    names[m].append(n)
plt.rcParams.update(
    {
        'svg.fonttype': 'none',
        'font.family': ['Arial', 'Nimbus Sans', 'DejaVu Sans'],
        'font.size': 7.5,
        'axes.labelsize': 8,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'legend.fontsize': 7,
        'axes.linewidth': 0.7,
        'xtick.major.width': 0.7,
        'ytick.major.width': 0.7,
        'xtick.major.size': 3,
        'ytick.major.size': 3,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'axes.spines.top': False,
        'axes.spines.right': False,
        'text.color': '#000000',
        'axes.edgecolor': '#000000',
        'axes.labelcolor': '#000000',
        'xtick.color': '#000000',
        'ytick.color': '#000000',
        'axes.grid': False,
        'legend.frameon': False,
        'legend.handlelength': 1.6,
    }
)


def psave(fig, name):
    fig.tight_layout(pad=0.4)
    fig.savefig(f'{OUT}/{name}.svg', bbox_inches='tight', pad_inches=0.03)
    fig.savefig(f'{OUT}/{name}.png', dpi=200, bbox_inches='tight', pad_inches=0.03)
    plt.close(fig)
    s = open(f'{OUT}/{name}.svg').read()
    s = re.sub(r"font-family:\s*[^;\"]+", 'font-family: Arial', s)
    open(f'{OUT}/{name}.svg', 'w').write(s)
    print('wrote', name)


x = np.arange(101)
fig, ax = plt.subplots(figsize=(3.0, 2.5))
numbers = {
    'set': SET,
    'settings': 'Picard 3.5.0 CollectRnaSeqMetrics, STRAND_SPECIFICITY NONE, MINIMUM_LENGTH 1000, refFlat protein-coding Ensembl 113 transcripts without ribosomal-protein and MT genes; primary mapped records, raw reads',
    'per_method': {},
}
for m in METH:
    Y = np.array(prof[m])
    mu = Y.mean(0)
    n = len(Y)
    ci = stats.t.ppf(0.975, n - 1) * Y.std(0, ddof=1) / np.sqrt(n)
    ax.fill_between(x, mu - ci, mu + ci, color=COL[m], alpha=0.25, lw=0)
    ax.plot(x, mu, color=COL[m], lw=2.0, label=f'{m} (n={n})')
    numbers['per_method'][m] = {
        'n': n,
        'samples': names[m],
        'profile_mean': mu.tolist(),
        'profile_ci95': ci.tolist(),
        'peak_percentile': int(np.argmax(mu)),
        'peak': float(mu.max()),
        'balance': bal[m],
        'balance_mean': float(np.mean(bal[m])),
        'balance_sd': float(np.std(bal[m], ddof=1)),
    }
ax.axhline(1, color='#7f7f7f', lw=0.8, ls=(0, (4, 3)))
ax.set_xlim(0, 100)
ax.set_ylim(0, None)
ax.set_xlabel("gene-body percentile (5'→3')")
ax.set_ylabel('normalized read coverage')
ax.legend(loc='upper left')
psave(fig, f'p6_picard_profile_{TAG}')
for style in ('plain', 'box'):
    fig, ax = plt.subplots(figsize=(2.3, 2.5))
    rng = np.random.default_rng(1)
    for i, m in enumerate(METH):
        v = np.array(bal[m])
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        if style == 'box':
            ax.add_patch(plt.Rectangle((i - 0.3, q1), 0.6, q3 - q1, color=COL[m], alpha=0.25, lw=0))
            ax.vlines(i, v.min(), v.max(), color=COL[m], lw=0.8)
        ax.scatter(i + rng.normal(0, 0.08, len(v)), v, s=9, color=COL[m], alpha=0.9, lw=0, zorder=3)
        ax.hlines(med, i - 0.3, i + 0.3, color=COL[m] if style == 'box' else '#000000', lw=1.4, zorder=4)
    ax.axhline(1, color='#7f7f7f', lw=0.8, ls=(0, (4, 3)))
    ax.set_xticks(range(len(METH)))
    ax.set_xticklabels(METH, rotation=30, ha='right')
    ax.set_ylim(0, 2)
    ax.set_ylabel("5'→3' balance (2 x coverage centroid)")
    psave(fig, f'p6b_picard_balance_{TAG}' + ('_box' if style == 'box' else ''))
json.dump(numbers, open(f'{OUT}/p6_picard_numbers_{TAG}.json', 'w'), indent=1)
VAL = f'{OUT}/values'
os.makedirs(VAL, exist_ok=True)  # raw values for replotting
with open(f'{VAL}/p6_picard_profile_{TAG}_per_sample.tsv', 'w') as f:
    f.write('method\tsample\t' + '\t'.join(f'pct{i:03d}' for i in range(101)) + '\n')
    for m in METH:
        for n, h in zip(names[m], prof[m]):
            f.write(f'{m}\t{n}\t' + '\t'.join(f'{v:.5f}' for v in h) + '\n')
with open(f'{VAL}/p6_picard_profile_{TAG}_mean_ci95.tsv', 'w') as f:
    f.write('gene_body_percentile\t' + '\t'.join(f'{m}_mean\t{m}_ci95_halfwidth' for m in METH) + '\n')
    for i in range(101):
        f.write(
            f'{i}\t'
            + '\t'.join(
                f'{numbers["per_method"][m]["profile_mean"][i]:.5f}\t{numbers["per_method"][m]["profile_ci95"][i]:.5f}'
                for m in METH
            )
            + '\n'
        )
with open(f'{VAL}/p6b_picard_balance_{TAG}.tsv', 'w') as f:
    f.write('method\tsample\tbalance_2x_centroid\n')
    for m in METH:
        for n, b in zip(names[m], bal[m]):
            f.write(f'{m}\t{n}\t{b:.5f}\n')
print('DONE')

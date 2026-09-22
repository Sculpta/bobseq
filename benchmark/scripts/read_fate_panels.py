#!/usr/bin/env python3
"""Starting statistics for the main figure, one SVG per panel:
  p0a_read_fate_<set>        what survives from alignment input to a filtered read, % of reads into STAR; mean over the benchmark wells with a 95% t-interval ribbon and an in-plot legend with n (per-well funnel from funnel_per_sample.py, per-well values in values/)
  p0b_reads_per_sample_<set> filtered reads per sample as sequenced (one dot per sample, bar = median; from stats_E_U/rarefied.tsv reads_native)
BM_COVSET=native (default; native read length) or 50nt (read-length matched)."""

from settings import WORK
from scipy import stats as _st
import os, csv, json, re, logging, numpy as np, matplotlib, matplotlib.ticker

matplotlib.use('Agg')
import matplotlib.pyplot as plt

logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
from matplotlib.lines import Line2D

P = WORK
U = f'{P}/benchmark_uniform'
SET = os.environ.get('BM_COVSET', 'native')
FIG = 'preprint_figures_native' if SET == 'native' else f'preprint_figures_{SET}'
OUT = f'{P}/{FIG}/main_figure_panels'
os.makedirs(OUT, exist_ok=True)
TAG = 'native' if SET == 'native' else 'matched'
plt.rcParams.update(
    {
        'svg.fonttype': 'none',
        'font.family': ['Arial', 'Nimbus Sans', 'DejaVu Sans'],
        'font.size': 7,
        'axes.labelsize': 7,
        'xtick.labelsize': 6.5,
        'ytick.labelsize': 6.5,
        'legend.fontsize': 6.5,
        'axes.linewidth': 0.6,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 2.5,
        'ytick.major.size': 2.5,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'axes.spines.top': False,
        'axes.spines.right': False,
        'text.color': '#231F20',
        'axes.edgecolor': '#231F20',
        'axes.labelcolor': '#231F20',
        'xtick.color': '#231F20',
        'ytick.color': '#231F20',
        'axes.grid': False,
    }
)
import os as _os, sys as _sys

_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from palette import COL, METH  # one palette shared by every panel


def save(fig, name):
    fig.tight_layout(pad=0.4)
    fig.savefig(f'{OUT}/{name}.svg', bbox_inches='tight', pad_inches=0.03)
    fig.savefig(f'{OUT}/{name}.png', dpi=200, bbox_inches='tight', pad_inches=0.03)
    plt.close(fig)
    s = open(f'{OUT}/{name}.svg').read()
    s = re.sub(r"font-family:\s*[^;\"]+", 'font-family: Arial', s)
    open(f'{OUT}/{name}.svg', 'w').write(s)
    print('wrote', name)


LEGEND = (
    os.environ.get('BM_PANEL_LEGEND', '0') == '1'
)  # legends off by default; one shared legend SVG (panel_legends.py) is drawn instead


def note(ax, *lines, loc='tl'):
    """small grey in-axes note stating the set and the matching rule (off by default: no titles, no notes)"""
    if os.environ.get('BM_PANEL_NOTES', '0') != '1':
        return  # notes off: the matching rule lives in the axis labels only
    x, y, ha, va = {
        'tl': (0.03, 0.97, 'left', 'top'),
        'tr': (0.97, 0.97, 'right', 'top'),
        'bl': (0.03, 0.04, 'left', 'bottom'),
        'br': (0.97, 0.04, 'right', 'bottom'),
        'tc': (0.5, 0.97, 'center', 'top'),
        'ml': (0.03, 0.42, 'left', 'center'),
        'mr': (0.97, 0.42, 'right', 'center'),
    }[loc]
    ax.text(
        x,
        y,
        '\n'.join(lines),
        transform=ax.transAxes,
        fontsize=5.5,
        color='#6B6B6B',
        ha=ha,
        va=va,
        linespacing=1.3,
        zorder=5,
    )


SETLBL = 'native read length (BOBseq 2 x 150 PE)' if SET == 'native' else '50-nt single reads, all methods'


def legend(fig):
    if not LEGEND:
        return
    fig.legend(
        handles=[Line2D([], [], color=COL[m], lw=1.6) for m in METH],
        labels=METH,
        loc='lower center',
        frameon=False,
        ncol=4,
        handlelength=1.4,
        columnspacing=1.0,
        bbox_to_anchor=(0.5, -0.06),
    )


SETS = {
    'matched': {'DRUG-seq': 'drugseq', 'BRB-seq': 'brbseq', 'prime-seq': 'primeseq', 'BOBseq': 'bob57_24plex_pe_hs'},
    'native': {
        'DRUG-seq': 'drugseq_native',
        'BRB-seq': 'brbseq_native',
        'prime-seq': 'primeseq_native',
        'BOBseq': 'bob57_24plex_pe_native',
    },
}[TAG]
# ---- p0a: read fate, % of reads into STAR, summed over the SAME benchmark wells at every step (the library-level STAR log counts every
#      aligned well of a plate while 'filtered' covers only the benchmark wells, so every step is taken from the per-well funnel of funnel_per_sample.py) ----
steps = [
    ('reads_in', 'reads\ninto STAR'),
    ('mapped', 'mapped\n(unique +\nmultimapping)'),
    ('unique', 'uniquely\nmapped'),
    ('filtered', 'unique and\nexonic protein-\ncoding (filtered)'),
]  # explicit step labels
arms = {}
fun_by_method = {}
for m in METH:
    rows = list(csv.DictReader(open(f'{U}/{SETS[m]}/stats_E_U/rarefied.tsv'), delimiter='\t'))
    usable = {r['well']: float(r['reads_native']) for r in rows if r.get('reads_native')}
    wm = {}
    if os.path.exists(f'{U}/{SETS[m]}/well_map.tsv'):  # raw barcode in the read name -> error-corrected well (DRUG-seq)
        for ln in open(f'{U}/{SETS[m]}/well_map.tsv'):
            raw, well = ln.rstrip('\n').split('\t')[:2]
            wm[raw] = well
    fun = {}
    for r in csv.DictReader(open(f'{U}/{SETS[m]}/funnel_per_sample.tsv'), delimiter='\t'):
        w = wm.get(r['sample'], r['sample'])
        d = fun.setdefault(w, {'reads_in': 0, 'mapped': 0, 'unique': 0})
        for k in d:
            d[k] += int(r[k])
    ic = f'{U}/{SETS[m]}/input_counts.tsv'
    pmf = f'{U}/{SETS[m]}/perwell_mapped.tsv'  # per-well mapped / unique counts of the benchmark wells
    if os.path.exists(ic) and os.path.exists(
        pmf
    ):  # BOBseq: the pooled BAM is mapped-only, so reads into STAR / mapped / unique come from the per-well tables (same wells)
        inp = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(ic) if l.strip()}
        pm = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open(pmf) if l.count('\t') >= 2}
        fun = {
            w: {'reads_in': inp[w], 'mapped': int(pm[w][1]), 'unique': int(pm[w][2])}
            for w in usable
            if w in inp and w in pm
        }
    miss = [w for w in usable if w not in fun]
    assert not miss, f'{m}: wells without funnel counts: {miss[:5]}'
    arms[m] = {
        'n_samples': len(usable),
        'reads_in': sum(fun[w]['reads_in'] for w in usable),
        'mapped': sum(fun[w]['mapped'] for w in usable),
        'unique': sum(fun[w]['unique'] for w in usable),
        'filtered': sum(usable.values()),
    }
    fun_by_method[m] = {w: {**fun[w], 'filtered': usable[w]} for w in usable}
fig, ax = plt.subplots(figsize=(2.8, 2.2))
out = {}
bw = 0.19
perwell = {}
for i, m in enumerate(METH):
    a = arms[m]
    y = [100 * a[k] / a['reads_in'] for k, _ in steps]
    out[m] = dict(zip([k for k, _ in steps], y))
    out[m]['reads_in_M'] = a['reads_in'] / 1e6
    out[m]['n_samples'] = a['n_samples']
    pw = {
        w: [100 * fun_by_method[m][w][k] / fun_by_method[m][w]['reads_in'] for k, _ in steps] for w in fun_by_method[m]
    }
    perwell[m] = pw  # each replicate's own funnel
    Y = np.array([pw[w] for w in pw])
    mu = Y.mean(axis=0)
    k = len(Y)
    ci = (
        _st.t.ppf(0.975, k - 1) * Y.std(axis=0, ddof=1) / np.sqrt(k) if k >= 3 else np.zeros(len(steps))
    )  # mean over wells, 95% t-interval ribbon
    ax.fill_between(range(len(steps)), mu - ci, mu + ci, color=COL[m], alpha=0.3, lw=0, zorder=2)
    ax.plot(range(len(steps)), mu, color=COL[m], lw=1.3, zorder=3)
ax.set_xticks(range(len(steps)))
ax.set_xticklabels([l for _, l in steps], fontsize=6, linespacing=0.95)
ax.set_ylim(0, 100)
ax.set_xlim(-0.15, len(steps) - 0.85)
ax.set_ylabel('% of reads into STAR (per sample)')
ax.legend(
    handles=[Line2D([], [], color=COL[m], lw=1.3) for m in METH],
    labels=[f'{m} (n={len(perwell[m])})' for m in METH],
    loc='lower left',
    frameon=False,
    fontsize=6.5,
    handlelength=1.6,
    borderaxespad=0.3,
    labelspacing=0.3,
)
save(fig, f'p0a_read_fate_{TAG}')
json.dump({**out, '_arms': arms, '_perwell': perwell}, open(f'{OUT}/p0a_read_fate_{TAG}.json', 'w'), indent=1)
print(json.dumps({m: {k: round(v, 1) for k, v in out[m].items()} for m in METH}))
# ---- p0b: filtered reads per sample as sequenced ----
SETS = {
    'matched': {'DRUG-seq': 'drugseq', 'BRB-seq': 'brbseq', 'prime-seq': 'primeseq', 'BOBseq': 'bob57_24plex_pe_hs'},
    'native': {
        'DRUG-seq': 'drugseq_native',
        'BRB-seq': 'brbseq_native',
        'prime-seq': 'primeseq_native',
        'BOBseq': 'bob57_24plex_pe_native',
    },
}[TAG]
fig, ax = plt.subplots(figsize=(2.2, 2.2))
nums = {}
for i, m in enumerate(METH):
    rows = list(csv.DictReader(open(f'{U}/{SETS[m]}/stats_E_U/rarefied.tsv'), delimiter='\t'))
    v = np.array(sorted({(r['well'], float(r['reads_native'])) for r in rows if r.get('reads_native')}))[:, 1].astype(
        float
    )
    jit = np.random.default_rng(1).normal(0, 0.07, len(v))
    ax.scatter(i + jit, v, s=7, color=COL[m], alpha=0.8, lw=0, zorder=2)
    ax.hlines(np.median(v), i - 0.28, i + 0.28, color='#231F20', lw=1.0, zorder=3)
    nums[m] = dict(n=len(v), median=float(np.median(v)), min=float(v.min()), max=float(v.max()))
ax.set_yscale('log')
ax.set_yticks([2e5, 5e5, 1e6, 2e6, 5e6])
ax.set_yticklabels(['200k', '500k', '1M', '2M', '5M'])
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks(range(len(METH)))
ax.set_xticklabels(METH, rotation=30, ha='right')
note(ax, 'no subsampling; (n) = samples', loc='tl')
ax.set_xlim(-0.6, len(METH) - 0.4)
ax.set_ylabel('filtered reads per sample\n(as sequenced)')
ax.set_ylim(1.2e5, 1.2e7)
save(fig, f'p0b_reads_per_sample_{TAG}')
json.dump(nums, open(f'{OUT}/p0b_reads_per_sample_{TAG}.json', 'w'), indent=1)
print(json.dumps(nums))
# ---- p0c / p0d: reads per sample after de-duplication, and % duplicate rate ----
# per well at native depth: filtered reads (stats_E_U/basic/reads_per_well.tsv) and molecules = UCI (perwell_B.tsv B_corr, the currency of p1);
# duplicate rate = 100 x (1 - UCI / filtered reads). p0d2 = the same rate at a matched depth of 250k filtered reads (rarefied.tsv), because the
# native-depth rate rises with how deep a well was sequenced (BOBseq wells carry ~10x the filtered reads of the competitor samples).
os.makedirs(f'{OUT}/values', exist_ok=True)
D250 = 250000
dd = {}
for m in METH:
    B = {r['well']: r for r in csv.DictReader(open(f'{U}/{SETS[m]}/stats_E_U/basic/perwell_B.tsv'), delimiter='\t')}
    R = {
        l.split('\t')[0]: float(l.split('\t')[1])
        for l in open(f'{U}/{SETS[m]}/stats_E_U/basic/reads_per_well.tsv')
        if l.strip()
    }
    rar = {
        r['well']: float(r['B_corr'])
        for r in csv.DictReader(open(f'{U}/{SETS[m]}/stats_E_U/rarefied.tsv'), delimiter='\t')
        if r.get('B_corr') and r['at_native'] == '0' and int(r['depth']) == D250
    }
    dd[m] = {
        w: dict(
            filtered_reads=R[w],
            molecules_UCI=float(B[w]['B_corr']),
            duplicate_rate_pct=100 * (1 - float(B[w]['B_corr']) / R[w]),
            duplicate_rate_pct_at_250k=(100 * (1 - rar[w] / D250) if w in rar else None),
        )
        for w in B
        if w in R
    }
with open(f'{OUT}/values/p0c_p0d_dedup_per_sample_{TAG}.tsv', 'w') as f:
    f.write('method\tsample\tfiltered_reads\tmolecules_UCI\tduplicate_rate_pct\tduplicate_rate_pct_at_250k\n')
    for m in METH:
        for w, v in dd[m].items():
            d250 = v['duplicate_rate_pct_at_250k']
            f.write(
                f"{m}\t{w}\t{v['filtered_reads']:.0f}\t{v['molecules_UCI']:.0f}\t{v['duplicate_rate_pct']:.2f}\t{'' if d250 is None else f'{d250:.2f}'}\n"
            )


def strip_panel(name, key, ylabel, log=False, ylim=None, yticks=None, note_lines=()):
    fig, ax = plt.subplots(figsize=(2.2, 2.2))
    nums = {}
    for i, m in enumerate(METH):
        v = np.array([x[key] for x in dd[m].values() if x[key] is not None], float)
        jit = np.random.default_rng(1).normal(0, 0.07, len(v))
        ax.scatter(i + jit, v, s=7, color=COL[m], alpha=0.8, lw=0, zorder=2)
        ax.hlines(np.median(v), i - 0.28, i + 0.28, color='#231F20', lw=1.0, zorder=3)
        nums[m] = dict(
            n=int(len(v)), median=float(np.median(v)), mean=float(v.mean()), min=float(v.min()), max=float(v.max())
        )
    if log:
        ax.set_yscale('log')
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    if yticks:
        ax.set_yticks(yticks[0])
        ax.set_yticklabels(yticks[1])
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_xticks(range(len(METH)))
    ax.set_xticklabels(METH, rotation=30, ha='right')
    ax.set_xlim(-0.6, len(METH) - 0.4)
    ax.set_ylabel(ylabel)
    note(ax, *note_lines, loc='tl')
    save(fig, name)
    json.dump(nums, open(f'{OUT}/{name}.json', 'w'), indent=1)
    print(name, json.dumps({m: round(nums[m]['median'], 1) for m in METH}))


strip_panel(
    f'p0c_reads_after_dedup_{TAG}',
    'molecules_UCI',
    'filtered reads per sample\nafter de-duplication (molecules, UMI)',
    log=True,
    yticks=([1e5, 2e5, 5e5, 1e6, 2e6], ['100k', '200k', '500k', '1M', '2M']),
    note_lines=('no subsampling; (n) = samples',),
)
strip_panel(
    f'p0d_duplicate_rate_{TAG}',
    'duplicate_rate_pct',
    'duplicate rate (%)\n1 - molecules / filtered reads, as sequenced',
    ylim=(0, 100),
    note_lines=('as sequenced; (n) = samples',),
)
strip_panel(
    f'p0d2_duplicate_rate_at_250k_{TAG}',
    'duplicate_rate_pct_at_250k',
    'duplicate rate (%) at 250k filtered reads\n1 - molecules / reads, matched depth',
    ylim=(0, 100),
    note_lines=('every sample subsampled to 250k', 'filtered reads; (n) = samples'),
)


# ---- combined panels: two dot columns per method, filled = as sequenced, open = after de-duplication / at matched depth ----
def paired_panel(
    name,
    keys,
    labels,
    ylabel,
    log=False,
    ylim=None,
    yticks=None,
    note_lines=(),
    legend_loc='upper left',
    ticklabels=None,
):
    fig, ax = plt.subplots(figsize=(2.6, 2.2))
    nums = {}
    off = (-0.19, 0.19)
    for i, m in enumerate(METH):
        nums[m] = {}
        for k, (key, o) in enumerate(zip(keys, off)):
            v = np.array([x[key] for x in dd[m].values() if x[key] is not None], float)
            jit = np.random.default_rng(1 + k).normal(0, 0.05, len(v))
            if k == 0:
                ax.scatter(i + o + jit, v, s=7, color=COL[m], alpha=0.8, lw=0, zorder=2)
            else:
                ax.scatter(i + o + jit, v, s=8, facecolors='white', edgecolors=COL[m], lw=0.6, alpha=0.9, zorder=2)
            ax.hlines(np.median(v), i + o - 0.15, i + o + 0.15, color='#231F20', lw=1.0, zorder=3)
            nums[m][key] = dict(
                n=int(len(v)), median=float(np.median(v)), mean=float(v.mean()), min=float(v.min()), max=float(v.max())
            )
    if log:
        ax.set_yscale('log')
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    if yticks:
        ax.set_yticks(yticks[0])
        ax.set_yticklabels(yticks[1])
    if ylim:
        ax.set_ylim(*ylim)
    if (
        ticklabels
    ):  # one tick per column saying what is plotted, the method name once under each pair
        ax.set_xticks([i + o for i in range(len(METH)) for o in off])
        ax.set_xticklabels([t for _ in METH for t in ticklabels], fontsize=5.6, linespacing=0.95)
        ax.tick_params(axis='x', length=2)
        for i, m in enumerate(METH):
            ax.text(
                i,
                -0.2,
                m,
                transform=ax.get_xaxis_transform(),
                ha='center',
                va='top',
                fontsize=7,
                color=COL[m],
                fontweight='bold',
                clip_on=False,
            )
    else:
        ax.set_xticks(range(len(METH)))
        ax.set_xticklabels(METH, rotation=30, ha='right')
    ax.set_xlim(-0.6, len(METH) - 0.4)
    ax.set_ylabel(ylabel)
    ax.legend(
        handles=[
            Line2D([], [], marker='o', ls='', color='#555', markersize=3.5),
            Line2D([], [], marker='o', ls='', markerfacecolor='white', markeredgecolor='#555', markersize=3.8),
        ],
        labels=labels,
        loc=legend_loc,
        frameon=False,
        fontsize=6,
        handletextpad=0.3,
        borderaxespad=0.3,
        labelspacing=0.3,
    )
    note(ax, *note_lines, loc='tl')
    save(fig, name)
    json.dump(nums, open(f'{OUT}/{name}.json', 'w'), indent=1)
    print(name, json.dumps({m: {k: round(nums[m][k]['median'], 1) for k in keys} for m in METH}))


paired_panel(
    f'p0cd_reads_before_after_dedup_{TAG}',
    ('filtered_reads', 'molecules_UCI'),
    ('filtered reads (as sequenced)', 'after dedup (molecules, UMI)'),
    'reads per sample',
    log=True,
    yticks=([1e5, 2e5, 5e5, 1e6, 2e6, 5e6], ['100k', '200k', '500k', '1M', '2M', '5M']),
    note_lines=('no subsampling; bar = median',),
    legend_loc='upper left',
)  # legend top-left, away from the DRUG-seq dots
paired_panel(
    f'p0dd_duplicate_rate_{TAG}',
    ('duplicate_rate_pct', 'duplicate_rate_pct_at_250k'),
    ('as sequenced', 'at 250k filtered reads (matched depth)'),
    'duplicate rate (%)\n1 - molecules / filtered reads',
    ylim=(0, 100),
    note_lines=('bar = median',),
)

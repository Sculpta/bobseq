#!/usr/bin/env python3
"""p3 (read composition) and p3b (intronic reads per sample) from the rule-based per-sample composition (composition_rules.py:
the preprint's classification rules on the per-sample BAMs, multimappers included, native BOBseq mate 1). The read-table
versions from main_figure_panels.py are written to extra/ as *_readtable_*. Writes the panels, their json sidecars and values/ TSVs.
BM_COVSET=native|50nt."""

from settings import WORK
import os, re, sys, json, csv, numpy as np, matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

P = WORK
SET = os.environ.get('BM_COVSET', 'native')
TAG = 'native' if SET == 'native' else 'matched'
OUT = f'{P}/{"preprint_figures_native" if SET == "native" else f"preprint_figures_{SET}"}/main_figure_panels'
VAL = f'{OUT}/values'
os.makedirs(VAL, exist_ok=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from palette import COL, METH, CLS

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


def save(fig, name):
    fig.tight_layout(pad=0.4)
    fig.savefig(f'{OUT}/{name}.svg', bbox_inches='tight', pad_inches=0.03)
    fig.savefig(f'{OUT}/{name}.png', dpi=200, bbox_inches='tight', pad_inches=0.03)
    plt.close(fig)
    s = open(f'{OUT}/{name}.svg').read()
    s = re.sub(r"font-family:\s*[^;\"]+", 'font-family: Arial', s)
    open(f'{OUT}/{name}.svg', 'w').write(s)
    print('wrote', name)


rows = list(csv.DictReader(open(f'{VAL}/p3_composition_rules_{TAG}_per_sample.tsv'), delimiter='\t'))
keys = [k for k, _, _ in CLS]
per = {m: [r for r in rows if r['method'] == m] for m in METH}
comp = {
    m: {k: 100 * sum(int(r[k]) for r in per[m]) / sum(int(r['mapped_records']) for r in per[m]) for k in keys}
    for m in METH
}
fig, ax = plt.subplots(figsize=(3.0, 1.9))
left = np.zeros(len(METH))
for k, lab, c in CLS:
    v = np.array([comp[m][k] for m in METH])
    ax.barh(range(len(METH)), v, left=left, color=c, edgecolor='white', lw=0.4, label=lab, height=0.7)
    left += v
ax.set_yticks(range(len(METH)))
ax.set_yticklabels(METH)
ax.invert_yaxis()
ax.set_xlabel('% of mapped reads (multimappers included)')
ax.set_xlim(0, 100)
save(fig, f'p3_read_composition_{TAG}')
json.dump(
    {
        'classes': keys,
        'per_method_pct_of_mapped': comp,
        'n_samples': {m: len(per[m]) for m in METH},
        'source': f'values/p3_composition_rules_{TAG}_per_sample.tsv (composition_rules.py)',
    },
    open(f'{OUT}/p3_read_composition_{TAG}.json', 'w'),
    indent=1,
)
with open(f'{VAL}/p3_read_composition_{TAG}.tsv', 'w') as f:
    f.write('method\tn_samples\t' + '\t'.join(keys) + '\n')
    for m in METH:
        f.write(f'{m}\t{len(per[m])}\t' + '\t'.join(f'{comp[m][k]:.4f}' for k in keys) + '\n')
pw = {m: {r['sample']: 100 * int(r['intronic']) / int(r['mapped_records']) for r in per[m]} for m in METH}
fig, ax = plt.subplots(figsize=(2.2, 2.2))
for i, m in enumerate(METH):
    v = np.array(list(pw[m].values()))
    jit = np.random.default_rng(1).normal(0, 0.07, len(v))
    ax.scatter(i + jit, v, s=7, color=COL[m], alpha=0.8, lw=0, zorder=2)
    ax.hlines(np.median(v), i - 0.28, i + 0.28, color='#231F20', lw=1.0, zorder=3)
ax.set_xticks(range(len(METH)))
ax.set_xticklabels(METH, rotation=30, ha='right')
ax.set_ylabel('intronic reads (% of mapped reads)\nper sample')
ax.set_xlim(-0.6, len(METH) - 0.4)
ax.set_ylim(bottom=0)
save(fig, f'p3b_intronic_per_sample_{TAG}')
json.dump(pw, open(f'{OUT}/p3b_intronic_per_sample_{TAG}.json', 'w'), indent=1)
with open(f'{VAL}/p3b_intronic_per_sample_{TAG}.tsv', 'w') as f:
    f.write('method\tsample\tintronic_pct_of_mapped\n')
    for m in METH:
        for s_, v in pw[m].items():
            f.write(f'{m}\t{s_}\t{v:.4f}\n')
print(
    SET,
    {
        m: f'mRNA {comp[m]["mRNA"]:.1f} RP {comp[m]["ribosomal-protein"]:.1f} rRNA {comp[m]["rRNA"]:.1f} intronic {comp[m]["intronic"]:.1f}'
        for m in METH
    },
)
print('DONE')

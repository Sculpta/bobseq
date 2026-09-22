#!/usr/bin/env python3
"""p0_samples_<tag>: the benchmark's sample table as a panel: per method the number of replicates, the cell line and
the read layout of the set. Counts come from the data (BOBseq: wells_keep minus units_excluded; competitors: metadata in_benchmark).
BM_COVSET=native|50nt."""

from settings import WORK, RUN_JSON
import os, re, sys, json, matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

P = WORK
PB = f'{P}/benchmark_uniform/per_sample_bams'
SET = os.environ.get('BM_COVSET', 'native')
TAG = 'native' if SET == 'native' else 'matched'
OUT = f'{P}/{"preprint_figures_native" if SET == "native" else f"preprint_figures_{SET}"}/main_figure_panels'
VAL = f'{OUT}/values'
os.makedirs(VAL, exist_ok=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from palette import COL, METH

SETNAME = 'native' if SET == 'native' else 'matched_50nt'
n = {}
for l in open(f'{PB}/metadata.tsv'):
    f = l.rstrip('\n').split('\t')
    if (
        f[0] in ('drugseq', 'brbseq', 'primeseq')
        and {'drugseq': 'DRUG-seq', 'brbseq': 'BRB-seq', 'primeseq': 'prime-seq'}[f[0]] in METH
        and f[1] == SETNAME
        and f[8] == 'yes'
    ):
        m = {'drugseq': 'DRUG-seq', 'brbseq': 'BRB-seq', 'primeseq': 'prime-seq'}[f[0]]
        n[m] = n.get(m, 0) + 1
RUN = json.load(open(RUN_JSON))
ARM = f'{P}/benchmark_uniform/bob57_24plex_pe_native'
excl = {l.split('\t')[0] for i, l in enumerate(open(f'{ARM}/units_excluded.tsv')) if i > 0}
n['BOBseq'] = sum(
    1
    for w in (l.strip() for l in open(f'{ARM}/wells_keep.txt') if l.strip())
    if w not in excl and RUN['tso_species_map'][w] == 'human'
)
CELL = {'DRUG-seq': 'U-2 OS', 'BRB-seq': 'HEK293T', 'prime-seq': 'HEK293T', 'BOBseq': 'HEK293T'}
READS = (
    {
        'DRUG-seq': '52 nt, single',
        'BRB-seq': '80 nt, single',
        'prime-seq': '50 nt, single',
        'BOBseq': '2 x 150 nt, paired',
    }
    if SET == 'native'
    else {m: '50 nt, single' for m in METH}
)
UNIT = {'DRUG-seq': 'wells', 'BRB-seq': 'wells', 'prime-seq': 'samples', 'BOBseq': 'wells'}
plt.rcParams.update(
    {
        'svg.fonttype': 'none',
        'font.family': ['Arial', 'Nimbus Sans', 'DejaVu Sans'],
        'font.size': 7,
        'text.color': '#231F20',
    }
)
fig = plt.figure(figsize=(2.2, 2.2))
ax = fig.add_axes([0, 0, 1, 1])
ax.axis('off')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)  # square, one block per method
blk = 1.0 / len(METH)
x0 = 0.10
for i, m in enumerate(METH):
    top = 1.0 - i * blk
    y1, y2, y3 = top - 0.22 * blk, top - 0.52 * blk, top - 0.80 * blk
    ax.add_patch(plt.Rectangle((0.02, y3 - 0.09 * blk), 0.045, (y1 - y3) + 0.18 * blk, color=COL[m], lw=0))
    ax.text(x0, y1, m, va='center', ha='left', color=COL[m], fontweight='bold', fontsize=7.5)
    ax.text(x0, y2, f'{n[m]} {UNIT[m]}, {CELL[m]}', va='center', ha='left', fontsize=6.8)
    ax.text(x0, y3, READS[m], va='center', ha='left', fontsize=6.8, color='#444444')
    if i < len(METH) - 1:
        ax.plot([0.02, 0.98], [top - blk] * 2, color='#C8C8C8', lw=0.5)
name = f'p0_samples_{TAG}'
fig.savefig(f'{OUT}/{name}.svg', bbox_inches='tight', pad_inches=0.03)
fig.savefig(f'{OUT}/{name}.png', dpi=200, bbox_inches='tight', pad_inches=0.03)
plt.close(fig)
s = open(f'{OUT}/{name}.svg').read()
s = re.sub(r"font-family:\s*[^;\"]+", 'font-family: Arial', s)
open(f'{OUT}/{name}.svg', 'w').write(s)
with open(f'{VAL}/{name}.tsv', 'w') as f:
    f.write('method\treplicates\tunit\tcell_line\treads\n')
    for m in METH:
        f.write(f'{m}\t{n[m]}\t{UNIT[m]}\t{CELL[m]}\t{READS[m]}\n')
print('wrote', name, n)

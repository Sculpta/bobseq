#!/usr/bin/env python3
"""Standalone legends for the main figure (one legend is placed per figure in Illustrator): method colours (horizontal, vertical) and
the read-composition classes. Written to both figure sets' main_figure_panels/."""

from settings import WORK
import os, re, matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

P = WORK
plt.rcParams.update(
    {
        'svg.fonttype': 'none',
        'font.family': ['Arial', 'Nimbus Sans', 'DejaVu Sans'],
        'font.size': 7,
        'legend.fontsize': 7,
        'text.color': '#231F20',
    }
)
import os as _os, sys as _sys

_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from palette import COL, METH  # one palette shared by every panel
from palette import CLS as _CLS

CLS = [(lab, c) for _, lab, c in _CLS]  # the rule-based composition classes


def save(fig, out, name):
    fig.savefig(f'{out}/{name}.svg', bbox_inches='tight', pad_inches=0.02)
    fig.savefig(f'{out}/{name}.png', dpi=200, bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)
    s = open(f'{out}/{name}.svg').read()
    s = re.sub(r"font-family:\s*[^;\"]+", 'font-family: Arial', s)
    open(f'{out}/{name}.svg', 'w').write(s)
    print('wrote', out.split('/')[-2], name)


for fig_set in ('preprint_figures_native', 'preprint_figures_50nt'):
    out = f'{P}/{fig_set}/main_figure_panels'
    os.makedirs(out, exist_ok=True)
    for name, ncol in (('legend_methods_horizontal', 4), ('legend_methods_vertical', 1)):
        fig = plt.figure(figsize=(0.1, 0.1))
        fig.legend(
            handles=[Line2D([], [], color=COL[m], lw=1.6) for m in METH],
            labels=METH,
            loc='center',
            frameon=False,
            ncol=ncol,
            handlelength=1.4,
            columnspacing=1.0,
            labelspacing=0.5,
        )
        save(fig, out, name)
    for name, ncol in (('legend_methods_dots_horizontal', 4), ('legend_methods_dots_vertical', 1)):
        fig = plt.figure(figsize=(0.1, 0.1))
        fig.legend(
            handles=[Line2D([], [], color=COL[m], marker='o', ms=4, lw=0) for m in METH],
            labels=METH,
            loc='center',
            frameon=False,
            ncol=ncol,
            handlelength=1.0,
            columnspacing=1.0,
            labelspacing=0.5,
        )
        save(fig, out, name)
    for name, ncol in (('legend_composition_horizontal', 3), ('legend_composition_vertical', 1)):
        fig = plt.figure(figsize=(0.1, 0.1))
        fig.legend(
            handles=[Patch(facecolor=c, edgecolor='white', lw=0.4) for _, c in CLS],
            labels=[l for l, _ in CLS],
            loc='center',
            frameon=False,
            ncol=ncol,
            handlelength=1.0,
            columnspacing=1.2,
            labelspacing=0.5,
        )
        save(fig, out, name)

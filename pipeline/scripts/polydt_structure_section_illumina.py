#!/usr/bin/env python3
"""Reusable figure for the redesigned BOB-polydT "Read Structure" section.

render_polydt_structure_figure(all_results, bc_list, bc_short, out_dir) builds a
single 3-panel PNG:
  (1) construct schematic of the canonical / "Retained" molecule
  (2) per-BC stacked bar of read-structure outcome (Retained vs drop reasons)
  (3) per-BC polyT-length line plot (2-nt bins, 0-50 nt)

Data source: each BC's results['bob_polydt_retention'] written by
override_bob_polydt_illumina.py. Returns the PNG path, or None if no retention data.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import FancyBboxPatch

# Drop-reason buckets in stacked-bar order (top stack = Retained), with colors.
BUCKETS = [
    ('retained',    'Full structure',      '#2ca02c'),
    ('no_insert',   'no insert (<20bp)',   '#d62728'),
    ('short_polyt', 'short polyT (<5)',    '#9467bd'),
    ('short_c28',   'short C28 (<5bp)',    '#8c564b'),
    ('chimera',     'chimera (2 codes)',   '#7f7f7f'),
    ('c28_contam',  'C28 contam',          '#e377c2'),
    ('no_code',     'no code',             '#c7c7c7'),
]


def render_polydt_structure_figure(all_results, bc_list, bc_short, out_dir):
    rets = {bc: all_results[bc].get('bob_polydt_retention') for bc in bc_list}
    if not any(rets.values()):
        return None

    rcParams['font.family'] = 'Arial'
    rcParams['font.size'] = 9

    fig = plt.figure(figsize=(13, 9.5))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 2.4, 1.7], hspace=0.45)

    # ---- Panel 1: construct schematic --------------------------------------
    ax0 = fig.add_subplot(gs[0])
    ax0.set_xlim(0, 100); ax0.set_ylim(0, 10); ax0.axis('off')
    ax0.set_title('Full structure  (5′ → 3′)',
                  fontsize=11, fontweight='bold', loc='left')
    segs = [  # (label, width, color, sublabel)
        ('C28 primer', 16, '#1f77b4', '20 bp (≥5 req.)'),
        ('bobcode', 9, '#ff7f0e', '6mer (req.)'),
        ('polyT', 12, '#2ca02c', '≥5 (req.)'),
        ('cDNA insert', 34, '#dddddd', '≥20 bp (req.)'),
        ('TSO', 22, '#d62728', '25 bp (optional)'),
    ]
    x = 2
    for label, w, color, sub in segs:
        dashed = (label == 'cDNA insert')
        ax0.add_patch(FancyBboxPatch((x, 3.4), w, 3.0,
                      boxstyle='round,pad=0.1,rounding_size=0.3',
                      facecolor=color, edgecolor='black',
                      linewidth=0.8, linestyle='--' if dashed else '-'))
        tcol = 'white' if color in ('#1f77b4', '#2ca02c', '#d62728') else 'black'
        ax0.text(x + w / 2, 4.9, label, ha='center', va='center',
                 fontsize=9, fontweight='bold', color=tcol)
        ax0.text(x + w / 2, 2.4, sub, ha='center', va='top', fontsize=7, color='#333')
        x += w + 1
    ax0.text(2, 7.7, 'Full structure = code + polyT≥5 + insert ≥20bp + C28≥5bp  (TSO optional)',
             fontsize=8.5, style='italic', color='#222')

    # ---- Panel 2: per-BC stacked bar ---------------------------------------
    ax1 = fig.add_subplot(gs[1])
    xs = np.arange(len(bc_list))
    bottoms = np.zeros(len(bc_list))
    for key, label, color in BUCKETS:
        fracs = []
        for bc in bc_list:
            r = rets.get(bc) or {}
            tot = max(r.get('n_total', 0), 1)
            fracs.append((r.get('drop_reason', {}).get(key, 0)) / tot * 100)
        fracs = np.array(fracs)
        ax1.bar(xs, fracs, bottom=bottoms, color=color, label=label,
                edgecolor='white', linewidth=0.4, width=0.8)
        bottoms += fracs
    ax1.set_xticks(xs)
    ax1.set_xticklabels([bc_short[bc] for bc in bc_list])
    ax1.set_ylabel('% of reads')
    ax1.set_ylim(0, 100)
    ax1.set_title('Read-structure outcome per BC (stacked, % of total reads)',
                  fontsize=11, fontweight='bold', loc='left')
    ax1.legend(loc='center left', bbox_to_anchor=(1.005, 0.5),
               fontsize=8, frameon=False)
    # annotate Retained % at the base of each bar
    for i, bc in enumerate(bc_list):
        r = rets.get(bc) or {}
        tot = max(r.get('n_total', 0), 1)
        ret_pct = (r.get('funnel', {}).get('retained', 0)) / tot * 100
        ax1.text(i, ret_pct / 2, f'{ret_pct:.0f}%', ha='center', va='center',
                 fontsize=7.5, fontweight='bold', color='white')
    for s in ('top', 'right'):
        ax1.spines[s].set_visible(False)

    # ---- Panel 3: per-BC polyT-length lines (2-nt bins, 0-50 nt) -----------
    ax2 = fig.add_subplot(gs[2])
    BINW, XMAX = 2, 50
    centers = np.arange(0, XMAX, BINW) + BINW / 2.0   # 1, 3, ..., 49
    cmap = plt.get_cmap('tab10')
    for i, bc in enumerate(bc_list):
        hist = (rets.get(bc) or {}).get('polyt_hist', {}) or {}
        tot = sum(hist.values())
        if not tot:
            continue
        binned = np.zeros(len(centers))
        for k, v in hist.items():
            L = int(k)
            if 0 <= L < XMAX:
                binned[L // BINW] += v
        ax2.plot(centers, binned / tot * 100, '-', lw=1.3,
                 color=cmap(i % 10), label=bc_short[bc])
    ax2.axvline(5, color='#888', lw=1.0, ls='--')
    ax2.text(5.4, ax2.get_ylim()[1] * 0.95, 'polyT ≥5', fontsize=7,
             color='#888', va='top')
    ax2.set_xlim(0, XMAX)
    ax2.set_xlabel('polyT run length (nt, 2-nt bins)')
    ax2.set_ylabel('% of code-bearing reads')
    ax2.set_title('Observed polyT lengths per BC', fontsize=11,
                  fontweight='bold', loc='left')
    ax2.legend(loc='upper right', fontsize=7, frameon=False, ncol=2)
    for s in ('top', 'right'):
        ax2.spines[s].set_visible(False)

    out_png = os.path.join(out_dir, 'polydt_read_structure.png')
    try:
        fig.savefig(out_png, dpi=150, bbox_inches='tight', facecolor='white')
    finally:
        plt.close(fig)
    return out_png

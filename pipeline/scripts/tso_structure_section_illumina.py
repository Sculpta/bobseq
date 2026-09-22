#!/usr/bin/env python3
"""Reusable figure for the TSO-bob "Read Structure" section.

render_tso_structure_figure(all_results, bc_list, bc_short, out_dir) builds a
3-panel PNG:
  (1) construct schematic of the canonical TSO-bob molecule
  (2) per-BC stacked bar of read-structure outcome (QC-filtered vs drop reasons)
  (3) per-BC polyT/A-length line plot (2-nt bins, 0-50 nt)

Data source: each BC's results['tso7_retention'] written by
override_bob_tso_7mer_illumina.py. Returns the PNG path, or None if no retention data.

This is the TSO analog of polydt_structure_section_illumina.py — the species barcode is a
7-mer in the TSO (Nextera-R1 backbone), NOT a 6-mer in the polyT primer, so the
canonical molecule and the QC funnel are different.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import FancyBboxPatch

# Drop-reason buckets in stacked-bar order (top stack = QC-filtered), with colors.
# The RT-primer-end labels carry a {rt} placeholder so the figure can name the
# actual RT primer (C28, or R1 for the R1_polydT_VN groups).
BUCKETS = [
    ('retained',     '{rt}-pDt-insert-TSO',    '#2ca02c'),
    ('chimera',      '≥2 TSO (chimera)',       '#8c564b'),
    ('no_c28',       'no {rt} primer end',     '#ff7f0e'),
    ('short_polyt',  'polyT/A <10bp',          '#17becf'),
    ('no_insert',    'no STAR alignment',      '#d62728'),
    ('no_ggg',       'no GGG switch',          '#9467bd'),
    ('no_barcode',   'no 7-mer barcode',       '#c7c7c7'),
]


def render_tso_structure_figure(all_results, bc_list, bc_short, out_dir,
                                rets_key='tso7_retention',
                                title_suffix='',
                                out_name='tso_read_structure.png',
                                rt_label='C28'):
    """Render the §1 TSO structure figure. `rets_key` selects which retention
    dict to read (e.g. 'tso7_retention' for all reads, 'tso7_retention_no_rrna'
    for the rRNA-excluded variant). `title_suffix` is appended to the panel-2
    and panel-3 titles. `out_name` is the PNG filename. `rt_label` names the RT
    primer at the 3' end ('C28' or 'R1') and is substituted into all labels."""
    rets = {bc: all_results[bc].get(rets_key) for bc in bc_list}
    if not any(rets.values()):
        return None
    # RT-primer-end labels resolved for this group's primer.
    buckets = [(k, lbl.format(rt=rt_label), c) for k, lbl, c in BUCKETS]

    rcParams['font.family'] = 'Arial'
    rcParams['font.size'] = 9

    fig = plt.figure(figsize=(13, 9.5))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 2.4, 1.7], hspace=0.45)

    # ---- Panel 1: construct schematic --------------------------------------
    ax0 = fig.add_subplot(gs[0])
    ax0.set_xlim(0, 100); ax0.set_ylim(0, 10); ax0.axis('off')
    ax0.set_title(f'{rt_label}-pDt-insert-TSO read structure — TSO-bob  (5′ → 3′)',
                  fontsize=11, fontweight='bold', loc='left')
    segs = [  # (label, width, color, sublabel)
        ('Nextera-R1', 16, '#1f77b4', '≥5 bp (req.)'),
        ('7mer barcode', 10, '#ff7f0e', '7mer (req.)'),
        ('GGG', 6, '#9467bd', 'switch (req.)'),
        ('cDNA insert', 30, '#dddddd', '≥50 bp (req.)'),
        ('polyA/polyT', 8, '#2ca02c', '≥10 bp (req.)'),
        (rt_label, 20, '#d62728', '≥5 bp (req.)'),
    ]
    x = 2
    for label, w, color, sub in segs:
        dashed = (label == 'cDNA insert')
        ax0.add_patch(FancyBboxPatch((x, 3.4), w, 3.0,
                      boxstyle='round,pad=0.1,rounding_size=0.3',
                      facecolor=color, edgecolor='black',
                      linewidth=0.8, linestyle='--' if dashed else '-'))
        tcol = 'white' if color in ('#1f77b4', '#9467bd', '#2ca02c', '#d62728') else 'black'
        ax0.text(x + w / 2, 4.9, label, ha='center', va='center',
                 fontsize=8.5, fontweight='bold', color=tcol)
        ax0.text(x + w / 2, 2.4, sub, ha='center', va='top', fontsize=7, color='#333')
        x += w + 1
    ax0.text(2, 7.7, f'{rt_label}-pDt-insert-TSO = single 7-mer (≥5bp Nextera-R1) + GGG + insert '
             f'≥50bp + polyT/A≥10 + ≥5bp {rt_label}  (≥2-TSO chimeras excluded)',
             fontsize=8.0, style='italic', color='#222')

    # ---- Panel 2: per-BC stacked bar ---------------------------------------
    ax1 = fig.add_subplot(gs[1])
    xs = np.arange(len(bc_list))
    bottoms = np.zeros(len(bc_list))
    for key, label, color in buckets:
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
    ax1.set_xticklabels([bc_short.get(bc, bc) for bc in bc_list])
    ax1.set_ylabel('% of reads')
    ax1.set_ylim(0, 100)
    ax1.set_title('Read-structure outcome per BC (stacked, % of total reads)' + title_suffix,
                  fontsize=11, fontweight='bold', loc='left')
    ax1.legend(loc='center left', bbox_to_anchor=(1.005, 0.5),
               fontsize=8, frameon=False)
    for i, bc in enumerate(bc_list):
        r = rets.get(bc) or {}
        tot = max(r.get('n_total', 0), 1)
        ret_pct = (r.get('funnel', {}).get('retained', 0)) / tot * 100
        if ret_pct > 4:
            ax1.text(i, ret_pct / 2, f'{ret_pct:.0f}%', ha='center', va='center',
                     fontsize=7.5, fontweight='bold', color='white')
    for s in ('top', 'right'):
        ax1.spines[s].set_visible(False)

    # ---- Panel 3: per-BC polyT/A-length lines (2-nt bins, 0-50 nt) ----------
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
                 color=cmap(i % 10), label=bc_short.get(bc, bc))
    ax2.axvline(10, color='#888', lw=1.0, ls='--')
    ax2.text(10.4, ax2.get_ylim()[1] * 0.95, 'polyT ≥10', fontsize=7,
             color='#888', va='top')
    ax2.set_xlim(0, XMAX)
    ax2.set_xlabel('polyT/A run length (nt, 2-nt bins)')
    ax2.set_ylabel('% of barcode-bearing reads')
    ax2.set_title('Observed polyT/A lengths per BC' + title_suffix, fontsize=11,
                  fontweight='bold', loc='left')
    ax2.legend(loc='upper right', fontsize=7, frameon=False, ncol=2)
    for s in ('top', 'right'):
        ax2.spines[s].set_visible(False)

    out_png = os.path.join(out_dir, out_name)
    try:
        fig.savefig(out_png, dpi=150, bbox_inches='tight', facecolor='white')
    finally:
        plt.close(fig)
    return out_png

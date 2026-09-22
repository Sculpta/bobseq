#!/usr/bin/env python3
"""Figure for the "TSO Structure & Multiplicity" section (TSO-bob, pure-STAR).

render_tso_multiplicity_figure(all_results, bc_list, bc_short, out_dir) builds a
2-panel PNG from each BC's results['tso_multiplicity'] (written by
star_taxonomy._tso_multiplicity_stats):
  (A) per-BC TSO-unit count (no-TSO / 1 / 2 / >=3), % of all reads
  (B) per-BC barcode concordance among 2-TSO reads (different- vs same-species),
      with a black tick at the random-pairing expectation.

A TSO unit = Nextera-R1 backbone (>=5 bp) + a canonical 7-mer (<=1 mm). One TSO is
the canonical molecule; >=2 TSO are TSO-TSO chimeras, and a cross-species pair
(human + mouse) is inherently barcode-ambiguous.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.lines import Line2D


def render_tso_multiplicity_figure(all_results, bc_list, bc_short, out_dir):
    M = {bc: (all_results[bc].get('tso_multiplicity') or {}) for bc in bc_list}
    if not any(m.get('n_total') for m in M.values()):
        return None

    rcParams['font.family'] = 'Arial'
    rcParams['font.size'] = 9
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 5.2))
    xs = np.arange(len(bc_list))

    # ---- Panel A: TSO-unit count distribution (% of all reads) -------------
    def frac(bc, k):
        m = M[bc]; tot = max(m.get('n_total', 0), 1)
        return 100.0 * (m.get('count_dist', {}).get(k, 0)) / tot
    layers = [('1', '1 TSO (canonical)', '#2ca02c'),
              ('2', '2 TSO (chimera)', '#d62728'),
              ('3+', '≥3 TSO', '#7f1d1d'),
              ('0', 'no TSO (no barcode)', '#dddddd')]
    bottoms = np.zeros(len(bc_list))
    for k, lbl, col in layers:
        vals = np.array([frac(bc, k) for bc in bc_list])
        axA.bar(xs, vals, bottom=bottoms, color=col, label=lbl,
                edgecolor='white', linewidth=0.4, width=0.8)
        bottoms += vals
    for i, bc in enumerate(bc_list):
        v2 = frac(bc, '2') + frac(bc, '3+')
        if v2 > 3:
            axA.text(i, frac(bc, '1') + v2 / 2, f'{v2:.0f}%', ha='center', va='center',
                     fontsize=7.5, fontweight='bold', color='white')
    axA.set_xticks(xs); axA.set_xticklabels([bc_short.get(bc, bc) for bc in bc_list])
    axA.set_ylim(0, 100); axA.set_ylabel('% of reads')
    axA.set_title('TSO units per read (% of all reads)', fontsize=10,
                  fontweight='bold', loc='left')
    axA.legend(loc='center left', bbox_to_anchor=(1.005, 0.5), fontsize=8, frameon=False)
    for s in ('top', 'right'):
        axA.spines[s].set_visible(False)

    # ---- Panel B: 2-TSO barcode concordance --------------------------------
    # Pairing by CODE: any two samples, not only human vs mouse.
    diff = np.array([M[bc].get('diff_code2_pct') or 0 for bc in bc_list])
    same = np.array([M[bc].get('same_code2_pct') or 0 for bc in bc_list])
    axB.bar(xs, diff, color='#9467bd', label='two different barcodes',
            edgecolor='white', linewidth=0.4, width=0.8)
    axB.bar(xs, same, bottom=diff, color='#1f77b4', label='the same barcode twice',
            edgecolor='white', linewidth=0.4, width=0.8)
    for i, bc in enumerate(bc_list):
        e = M[bc].get('exp_diff_code2_pct')
        if e is not None:
            axB.plot([i - 0.4, i + 0.4], [e, e], color='black', lw=1.8, zorder=5)
        if M[bc].get('diff_code2_pct') is not None and diff[i] > 6:
            axB.text(i, diff[i] / 2, f'{diff[i]:.0f}%', ha='center', va='center',
                     fontsize=7.5, fontweight='bold', color='white')
    axB.set_xticks(xs); axB.set_xticklabels([bc_short.get(bc, bc) for bc in bc_list])
    axB.set_ylim(0, 100); axB.set_ylabel('% of 2-TSO reads')
    axB.set_title('2-TSO reads: are the two barcodes different?', fontsize=10, fontweight='bold', loc='left')
    handles = [
        Line2D([0], [0], color='#9467bd', lw=8, label='two different barcodes'),
        Line2D([0], [0], color='#1f77b4', lw=8, label='the same barcode twice'),
        Line2D([0], [0], color='black', lw=1.8, label='expected if molecules pair at random'),
    ]
    axB.legend(handles=handles, loc='center left', bbox_to_anchor=(1.005, 0.5),
               fontsize=8, frameon=False)
    for s in ('top', 'right'):
        axB.spines[s].set_visible(False)

    fig.tight_layout()
    out_png = os.path.join(out_dir, 'tso_multiplicity.png')
    try:
        fig.savefig(out_png, dpi=150, bbox_inches='tight', facecolor='white')
    finally:
        plt.close(fig)
    return out_png

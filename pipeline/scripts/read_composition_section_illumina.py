#!/usr/bin/env python3
"""Read start/end composition figure for the cross-barcode summary PDF.

Renders the 2x2 grid (true FIRST | true LAST ; AFTER ONT | BEFORE ONT) from the
per-barcode `read_composition` result dicts produced by read_end_composition.
Categories, order, and colours are imported from that module so the PDF section
never drifts from the classifier. Provenance is rendered separately as a native
reportlab table in analyze_speciesmix (crisper than baking text into the PNG)."""
import os

import read_end_composition_illumina as R

PANELS = [('true FIRST element', 'true_start'), ('true LAST element', 'true_end'),
          ('element AFTER the ONT adapter', 'after_ont'),
          ('element BEFORE the ONT adapter', 'before_ont')]
_INK, _INK2, _MUTED = '#0b0b0b', '#52514e', '#898781'


def render_grid(all_results, bc_list, bc_short, output_dir):
    """2x2 stacked-bar grid over all barcodes -> PNG path, or None if no data.
    Landscape aspect for a cross-barcode (landscape-letter) page."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    samples = [(bc_short[bc], all_results[bc].get('read_composition'))
               for bc in bc_list if all_results[bc].get('read_composition')]
    if not samples:
        return None
    ncol = len(samples)

    fig, axes = plt.subplots(2, 2, figsize=(11, 6.0))
    fig.patch.set_facecolor('white')
    for ax, (title, key) in zip(axes.ravel(), PANELS):
        ax.set_facecolor('white')
        for x, (lab, r) in enumerate(samples):
            n = r.get('n') or 1
            bottom = 0.0
            for g in R.ORDER:
                v = (r.get(key) or {}).get(g, 0)
                if not v:
                    continue
                ax.bar(x, 100 * v / n, bottom=bottom, width=0.86,
                       color=R.HUES[g], zorder=3)
                bottom += 100 * v / n
        ax.set_xticks(range(ncol))
        ax.set_xticklabels([lab for lab, _ in samples], fontsize=8, color=_INK2)
        ax.set_ylim(0, 100); ax.set_xlim(-0.6, ncol - 0.4)
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_yticklabels(['0', '25', '50', '75', '100%'], fontsize=7.5, color=_MUTED)
        ax.tick_params(length=0, pad=3)
        for s in ('top', 'right', 'left'):
            ax.spines[s].set_visible(False)
        ax.spines['bottom'].set_color('#c3c2b7')
        ax.set_title(title, fontsize=9.5, fontweight='bold', color=_INK, pad=6)

    from matplotlib.patches import Patch
    handles = [Patch(facecolor=R.HUES[g], label=g) for g in R.ORDER]
    fig.legend(handles=handles, loc='lower center', ncol=4, frameon=False,
               fontsize=7.6, labelcolor=_INK2, handlelength=1.0, handleheight=1.0,
               columnspacing=1.4, bbox_to_anchor=(0.5, -0.02))
    fig.subplots_adjust(left=0.055, right=0.99, top=0.93, bottom=0.13,
                        hspace=0.36, wspace=0.11)
    out = os.path.join(output_dir, 'read_composition_grid.png')
    fig.savefig(out, dpi=170, facecolor='white', bbox_inches='tight')
    plt.close(fig)
    return out

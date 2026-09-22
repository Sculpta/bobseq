#!/usr/bin/env python3
"""Keep the main-figure panel directory to the panels of the figure; everything else moves to extra/.
Kept: the documentation, the legends, p0a, p1/p2 depth curves of the set, p3, p3b, the threshold heatmap and profile,
the Picard profile and balance, and their number files. Idempotent; run after every panel generation.
BM_COVSET=native|50nt."""

from settings import WORK
import os, sys, shutil

P = WORK
SET = os.environ.get('BM_COVSET', 'native')
TAG = 'native' if SET == 'native' else 'matched'
D = f'{P}/{"preprint_figures_native" if SET == "native" else f"preprint_figures_{SET}"}/main_figure_panels'
X = f'{D}/extra'
os.makedirs(X, exist_ok=True)
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from panel_keep import keep

KEEP = keep(TAG)  # one list (panel_keep.py) shared with the generators
moved = 0
for f in sorted(os.listdir(D)):
    p = f'{D}/{f}'
    if os.path.isdir(p) or f.endswith('.md') or f.startswith('legend_'):
        continue
    stem = f.rsplit('.', 1)[0]
    if stem in KEEP:
        continue
    shutil.move(p, f'{X}/{f}')
    moved += 1
print(
    f'{SET}: {moved} files moved to extra/; main dir now holds',
    sum(1 for f in os.listdir(D) if not os.path.isdir(f'{D}/{f}')),
    'files',
)

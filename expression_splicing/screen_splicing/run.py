#!/usr/bin/env python3
"""BOBseq 24-plex differential splicing, five contrasts.   run.py [--only STAGE ...] [--list]"""
import argparse, os, subprocess, sys
P = os.path.dirname(os.path.abspath(__file__))
STAGES = [
    ("ds",        "five contrasts (CHX dose / CHX 50 / Ris 25 / ASO dose / ASO 100 nM): prefilter, Welch t or OLS dose slope, BH; positive hits; top 10 per contrast by dPSI within the padj gate", ["src/ds.py"]),
    ("volcano",   "per contrast: volcano without labels + volcano with the top 10 labelled (dodged labels, leader lines)", ["src/volcano.py"]),
    ("barcharts", "per contrast: top-10 PSI bar charts over every condition but Ris 500 mM (pooled fragment PSI, 95 % CI, well dots)", ["src/barcharts.py"]),
    ("loci",      "genome coordinates (junction, LSV extent, padded window) + BED of the top-10 junctions -> results/top10_loci.tsv", ["src/loci.py"]),
]
ap = argparse.ArgumentParser(description=__doc__); ap.add_argument("--only", nargs="+"); ap.add_argument("--list", action="store_true"); a = ap.parse_args()
if a.list:
    for n, d, _ in STAGES: print(f"  {n:<10} {d}")
    sys.exit(0)
for n, d, cmd in [s for s in STAGES if not a.only or s[0] in a.only]:
    print(f"== {n}: {d}", flush=True)
    if subprocess.run([sys.executable] + cmd, cwd=P).returncode: sys.exit(f"stage {n} failed")

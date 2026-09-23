#!/usr/bin/env python3
"""Browser-style coverage + junction-arc views of the CHX splicing loci and the positive controls.
   run.py [--only STAGE ...] [--list]"""
import argparse, os, subprocess, sys
P = os.path.dirname(os.path.abspath(__file__))
STAGES = [
    ("loci",    "config/loci.tsv: the CHX top-10 loci of ../screen_splicing (chx_dose, chx_high) + the CHX/NMD and risdiplam positive controls", ["src/loci.py"]),
    ("bams",    "depth-matched track BAMs (per-sample BAMs -> MAPQ 255 -> exact fragment subsampling; replicates + pseudobulk sets) -> .cache/bams/, results/track_depths.tsv", ["src/prepare_bams.py"]),
    ("figures", "per locus x track set: gene_with_introns, gene_exons_only, zoom, zoom_arcs (coverage tracks; arcs with fragment counts on zoom_arcs only) -> figures/<locus>_<tag>/, results/values/", ["src/gene_coverage_tracks.py"]),
    ("summary", "the LSV junction counts and shares per track behind every figure -> results/summary.tsv", ["src/summary.py"]),
]
ap = argparse.ArgumentParser(description=__doc__); ap.add_argument("--only", nargs="+"); ap.add_argument("--list", action="store_true"); a = ap.parse_args()
if a.list:
    for n, d, _ in STAGES: print(f"  {n:<8} {d}")
    sys.exit(0)
for n, d, cmd in [s for s in STAGES if not a.only or s[0] in a.only]:
    print(f"== {n}: {d}", flush=True)
    if subprocess.run([sys.executable] + cmd, cwd=P).returncode: sys.exit(f"stage {n} failed")

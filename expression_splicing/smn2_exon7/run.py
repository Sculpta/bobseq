#!/usr/bin/env python3
"""SMN2 exon 7: paralog-aware quantification.   run.py [--only STAGE ...] [--list]"""
import argparse, os, subprocess, sys
P = os.path.dirname(os.path.abspath(__file__))
STAGES = [
    ("blast_loci",  "BLAST the SMN1 vs SMN2 loci (30 kb each, from the reference) -> every diagnostic difference with coordinates in both", ["src/blast_loci.py"]),
    ("region_bams", "slice both loci (ALL MAPQ) out of every per-sample BAM, raw and UMI-deduplicated -> .cache/region_bams/", ["src/region_bams.py"]),
    ("fragments",   "per fragment: exon-7 splicing class + paralog call from the diagnostic sites (+ c.840 base)", ["src/fragments.py"]),
    ("blast_reads", "cross-check: blastn every read against both loci, compare best-hit calls with the site-based calls", ["src/blast_reads.py"]),
    ("psi",         "pooled / SMN2 / SMN1 PSI (shared skip pool) + c.840 T-fraction, per sample and per condition, with CIs", ["src/psi.py"]),
    ("simple",      "BLAST-native exon-7 inclusion index from reads at c.840 / c.*239 (no junctions, no skip attribution)", ["src/simple.py"]),
    ("figures",     "SMN2:SMN1 inclusion ratio (ex6>7) + junction reads per event (dedup, BOBseq PE)", ["src/figures.py"]),
]
ap = argparse.ArgumentParser(description=__doc__); ap.add_argument("--only", nargs="+"); ap.add_argument("--list", action="store_true"); a = ap.parse_args()
if a.list:
    for n, d, _ in STAGES: print(f"  {n:<12} {d}")
    sys.exit(0)
for n, d, cmd in [s for s in STAGES if not a.only or s[0] in a.only]:
    print(f"== {n}: {d}", flush=True)
    r = subprocess.run([sys.executable] + cmd, cwd=P)
    if r.returncode: sys.exit(f"stage {n} failed ({r.returncode})")

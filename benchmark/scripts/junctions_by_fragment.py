#!/usr/bin/env python3
"""Junction BED6 counted per fragment (distinct read name).

    junctions_by_fragment.py --bam q255.bam --sample <acc> --out <acc>_junctions.bed [--min-anchor 4]

A junction is the intron spanned by an N operation of a primary alignment. It is reported when the longest left and the longest
right anchor over all records spanning it are both >= min_anchor, where an anchor is the contiguous aligned run (M/=/X) flanking
the intron and any other operation (N, D, I, S) ends it: the same rule as regtools junctions extract -a, which filters by the
maximum anchor rather than per read. Support is the number of distinct read names, so a pair whose two mates both span the
junction counts once (per-record counting would count it twice). Output: chrom, intron start (0-based), intron end, sample,
count, strand ('?': the BAMs carry no XS tag), sorted by chrom, start, end.
"""

import argparse
from collections import defaultdict

import pysam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-anchor", type=int, default=4)
    a = ap.parse_args()
    frags = defaultdict(set)
    maxl = defaultdict(int)
    maxr = defaultdict(int)
    for r in pysam.AlignmentFile(a.bam, "rb").fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or "N" not in r.cigarstring:
            continue
        pos = r.reference_start
        blocks = []  # ('M'|'N'|'X', len, start, end) in CIGAR order
        for op, ln in r.cigartuples:
            if op in (0, 7, 8):
                blocks.append(("M", ln, pos, pos + ln))
                pos += ln
            elif op == 3:
                blocks.append(("N", ln, pos, pos + ln))
                pos += ln
            elif op == 2:
                blocks.append(("X", ln, pos, pos + ln))
                pos += ln  # deletion: breaks the anchor
            elif op == 1:
                blocks.append(("X", ln, pos, pos))  # insertion: breaks the anchor
            elif op == 4:
                blocks.append(("X", ln, pos, pos))  # soft clip: breaks the anchor
        for i, b in enumerate(blocks):
            if b[0] != "N":
                continue
            # anchor = the contiguous aligned (M/=/X) run immediately flanking the intron; any other op
            # (N, D, I, S) ends it -- this is regtools' definition (verified identical on SE and PE)
            left = right = 0
            for x in reversed(blocks[:i]):
                if x[0] != "M":
                    break
                left += x[1]
            for x in blocks[i + 1 :]:
                if x[0] != "M":
                    break
                right += x[1]
            k = (r.reference_name, b[2], b[3])
            frags[k].add(r.query_name)
            if left > maxl[k]:
                maxl[k] = left
            if right > maxr[k]:
                maxr[k] = right
    keep = [k for k in sorted(frags) if maxl[k] >= a.min_anchor and maxr[k] >= a.min_anchor]
    with open(a.out, "w") as out:
        for c, s, e in keep:
            out.write(f"{c}\t{s}\t{e}\t{a.sample}\t{len(frags[(c, s, e)])}\t?\n")
        if not keep:
            out.write(f".\t0\t0\t{a.sample}\t0\t?\n")


if __name__ == "__main__":
    main()

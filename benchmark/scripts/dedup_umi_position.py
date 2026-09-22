#!/usr/bin/env python3
"""UMI + position deduplication of one benchmark BAM (the `mol` count unit).

    dedup_umi_position.py --set <set> --in in.bam --out dedup.bam [--stats stats.tsv]

Keeps ONE read (or one read pair) per molecule, umi_tools-style: the molecule key is the UMI carried
in the read name (umi_of below) plus the alignment position, because BOBseq's 6 nt UMI
(4,096 values) cannot identify a molecule genome-wide on its own.

  single-end:  (contig, 5' end, strand, UMI)
  paired-end:  (contig, fragment start = min(pos, mate pos), |template length|, read1 strand, UMI)
               -- both mates of a pair produce the same key and share the read name, so the pair is
               kept or dropped as a unit in ONE pass: the first read of a new key registers its
               qname; a later record is written iff its qname was registered.
  a mate that is unmapped / on another contig falls back to the single-end key.

Input is coordinate-sorted (all benchmark BAMs are), so the seen/kept sets are flushed at every
contig change and memory stays bounded. Output order == input order (still sorted).
No UMI error-correction: with a 6 nt UMI, ed<=1 would merge 18 neighbours per UMI and over-collapse.
"""

import argparse
import os
import sys

import pysam

# Per-set facts about the per-sample BAMs (the sets of build_per_sample_bams.sh): method, read layout, paired, read length,
# UMI length, STAR --outFilter*OverLread. bobseq_pe_native alone was aligned at 0.33 (the pipeline's paired-end setting), the
# other sets at 0.66; documented in the sample sheet (star_filter), not corrected.
SETS = {
    "drugseq": ("drugseq", "50 nt SE (truncated)", False, 50, 10, 0.66),
    "drugseq_native": ("drugseq", "52 nt SE", False, 52, 10, 0.66),
    "primeseq": ("primeseq", "50 nt SE", False, 50, 10, 0.66),
    "primeseq_native": ("primeseq", "50 nt SE", False, 50, 10, 0.66),
    "bobseq_50nt": ("bobseq", "50 nt SE (R2 truncated)", False, 50, 6, 0.66),
    "bobseq_pe_native": ("bobseq", "2 x 150 nt PE", True, 150, 6, 0.33),
}


def umi_of(qname, set_name):
    """The UMI string carried in a read name, per set. None if it cannot be parsed."""
    if set_name == "bobseq_pe_native":
        # the pipeline writes <rid>__<umi>_<bobcode>_<glen>_<rt>: the umi field is UMI(6) + "HH"(2) spacer for the +8
        # codes, so the first 6 nt are the UMI (positions 7-8 are G-depleted, positions 1-6 uniform)
        parts = qname.split("__", 1)
        if len(parts) != 2:
            return None
        return parts[1].split("_", 1)[0][:6]
    # <rid>_<well>_<umi> (DRUG-seq, prime-seq, the 50-nt BOBseq set): the last field, cut to the set's UMI length. The
    # 50-nt BOBseq arm is prepped with the whole degenerate stretch after the bobcode in the read name (--umi full, up to
    # 14 nt); its UMI is the first 6 nt of that stretch, as in the paired-end set above and in the read table.
    return qname.rsplit("_", 1)[-1][: SETS[set_name][4]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True, choices=sorted(SETS))
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--stats", help="append one TSV line: set sample n_in n_out n_umi_unparsed")
    ap.add_argument("--sample", default="")
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.inp, "rb")
    tmp = args.out + ".tmp.bam"
    out = pysam.AlignmentFile(tmp, "wb", template=bam)

    seen, kept = set(), set()
    cur_tid = None
    n_in = n_out = n_unparsed = 0
    for r in bam.fetch(until_eof=True):
        n_in += 1
        if r.reference_id != cur_tid:
            seen.clear()
            kept.clear()
            cur_tid = r.reference_id
        umi = umi_of(r.query_name, args.set)
        if umi is None:
            n_unparsed += 1
            umi = ""
        paired = (
            r.is_paired and not r.mate_is_unmapped and r.next_reference_id == r.reference_id and r.template_length != 0
        )
        if paired:
            # read1's strand, computed identically from either mate (mate_is_reverse on read2 == read1's flag)
            strand = r.is_reverse if r.is_read1 else r.mate_is_reverse
            key = (min(r.reference_start, r.next_reference_start), abs(r.template_length), strand, umi)
        else:
            five = r.reference_end if r.is_reverse else r.reference_start
            key = (five, r.is_reverse, umi)
        if key not in seen:
            seen.add(key)
            kept.add(r.query_name)
        if r.query_name in kept:
            out.write(r)
            n_out += 1
    out.close()
    bam.close()
    os.replace(tmp, args.out)
    pysam.index(args.out)

    line = f"{args.set}\t{args.sample}\t{n_in}\t{n_out}\t{n_unparsed}\n"
    sys.stderr.write("dedup\t" + line)
    if args.stats:
        with open(args.stats, "a") as fh:
            fh.write(line)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build a Picard refFlat from the combined GTF.

EXPLORATORY TIER input. Nothing here reaches the Master Summary.

WHY NOT REUSE THE EXISTING HUMAN refFlat
A human-only refFlat built from the Ensembl GTF carries bare chromosome names (`1`),
while our combined BAMs use `HUMAN_1` / `MOUSE_17`. Picard matches refFlat contigs against
the BAM's sequence dictionary and reports ZERO coverage on a mismatch rather than failing,
so reusing it would produce a complete, plausible, entirely empty metrics file. Building
from combined_genome.gtf instead means the prefixes come along for free.

FORMAT (11 columns, no header; UCSC refFlat):
    geneName, transcriptName, chrom, strand,
    txStart, txEnd, cdsStart, cdsEnd,
    exonCount, exonStarts (CSV), exonEnds (CSV)

Coordinates are 0-based half-open, as UCSC expects; the GTF is 1-based inclusive, so every
start loses one. Getting that wrong shifts every transcript by a base and quietly changes
the coverage metrics rather than erroring.
"""

import argparse
import gzip
import re
import sys
from collections import defaultdict

_ATTR = re.compile(r'(\S+)\s+"([^"]*)"')


def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gtf', required=True)
    ap.add_argument('-o', '--out', required=True)
    ap.add_argument('--expect-prefixes', default='HUMAN_,MOUSE_',
                    help='comma-separated contig prefixes that must be present, or "" to skip')
    args = ap.parse_args()

    # transcript -> (gene, chrom, strand), exon list, cds list
    meta = {}
    exons = defaultdict(list)
    cds = defaultdict(list)

    with _open(args.gtf) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            feat = f[2]
            # stop_codon is folded into the CDS on purpose. Ensembl keeps it as its own
            # feature and EXCLUDES it from CDS; the UCSC genePred/refFlat convention
            # INCLUDES it. Without this the CDS is 3 bp short at the 3' end of every
            # coding transcript: compared with a UCSC-convention refFlat of the same
            # annotation, 46k transcripts come out off by 3, in both directions
            # depending on strand.
            if feat not in ('exon', 'CDS', 'stop_codon'):
                continue
            attrs = dict(_ATTR.findall(f[8]))
            tx = attrs.get('transcript_id')
            if not tx:
                continue
            start0 = int(f[3]) - 1          # GTF 1-based inclusive -> 0-based half-open
            end = int(f[4])
            if feat == 'exon':
                if tx not in meta:
                    # gene_id fallback, not the literal 'none' the existing human refFlat
                    # carries for 111,179 transcripts: Picard groups by geneName, so
                    # every unnamed transcript would otherwise collapse into a single
                    # pseudo-gene.
                    gene = (attrs.get('gene_name') or attrs.get('gene_id') or tx)
                    meta[tx] = (gene, f[0], f[6])
                exons[tx].append((start0, end))
            else:
                cds[tx].append((start0, end))     # CDS and stop_codon both land here

    if not meta:
        sys.exit(f'make_refflat: FATAL no transcripts parsed from {args.gtf}')

    per_prefix = defaultdict(int)
    n = 0
    with open(args.out, 'w') as out:
        for tx, (gene, chrom, strand) in meta.items():
            ex = sorted(exons[tx])
            tx_start = ex[0][0]
            tx_end = max(e for _s, e in ex)
            if cds[tx]:
                c = sorted(cds[tx])
                cds_start, cds_end = c[0][0], max(e for _s, e in c)
            else:
                # UCSC convention for a non-coding transcript: an empty CDS interval.
                # Picard treats cdsStart == cdsEnd as "no CDS"; using txStart/txEnd here
                # instead would make every lncRNA look like a coding gene.
                cds_start = cds_end = tx_end
            out.write('\t'.join([
                gene, tx, chrom, strand,
                str(tx_start), str(tx_end), str(cds_start), str(cds_end),
                str(len(ex)),
                ','.join(str(s) for s, _e in ex) + ',',
                ','.join(str(e) for _s, e in ex) + ',',
            ]) + '\n')
            n += 1
            for p in ('HUMAN_', 'MOUSE_'):
                if chrom.startswith(p):
                    per_prefix[p] += 1

    print(f'make_refflat: {n} transcripts -> {args.out}')
    for p, c in sorted(per_prefix.items()):
        print(f'make_refflat:   {p}* {c}')

    # A refFlat whose contigs do not match the BAM produces zero coverage silently, so
    # assert the prefixes are actually there rather than discovering it in a metrics file
    # full of zeros.
    want = [p for p in args.expect_prefixes.split(',') if p]
    missing = [p for p in want if per_prefix.get(p, 0) == 0]
    if missing:
        sys.exit(f'make_refflat: FATAL no transcripts on contigs prefixed '
                 f'{", ".join(missing)}. The GTF is not the combined build, or the '
                 f'prefixes changed; Picard would silently report zero coverage.')


if __name__ == '__main__':
    main()

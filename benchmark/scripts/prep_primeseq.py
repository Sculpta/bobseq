#!/usr/bin/env python3
"""prime-seq (E-MTAB-10142) -> single-end cDNA FASTQ with sample+UMI in the QNAME.

Same QNAME convention as the DRUG-seq prep, so every downstream script is shared
across all arms:  @<readid>_<BARCODE>_<UMI>

STRUCTURE, verified from the reads before writing this. The 16 nt read is BC6 + UMI10:
positions 1-6 carry ~0 bits of entropy and one 6-mer takes 91.7% of reads (these runs are
already demultiplexed per sample), while positions 7-16 sit at 2.0 bits, the maximum, with
169,014 distinct UMIs in 200,000 reads.

NOTE ON THE BARCODE VERSION. Janjic 2022 describes a 12 nt barcode + 16 nt UMI matched to
10x v3. E-MTAB-10142 is the Fig. 4 RNA-extraction comparison and predates that switch, so
it carries the older 6 + 10 layout. Measured, not assumed.
"""

import argparse, collections, gzip, sys


def op(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)


def main():
    a = argparse.ArgumentParser()
    a.add_argument('--bcumi', required=True)  # 16 nt read
    a.add_argument('--cdna', required=True)  # 50 nt read
    a.add_argument('--sample', required=True)
    a.add_argument('--out', required=True)
    a.add_argument('--bc-len', type=int, default=6)
    a.add_argument('--umi-len', type=int, default=10)
    a = a.parse_args()
    c = collections.Counter()
    bcs = collections.Counter()
    out = gzip.open(a.out, 'wt', compresslevel=1)
    with op(a.bcumi) as f1, op(a.cdna) as f2, out:
        while True:
            h1 = f1.readline()
            if not h1:
                break
            s1 = f1.readline().rstrip('\n')
            f1.readline()
            f1.readline()
            f2.readline()
            s2 = f2.readline().rstrip('\n')
            f2.readline()
            q2 = f2.readline().rstrip('\n')
            c['reads'] += 1
            if len(s1) < a.bc_len + a.umi_len:
                c['short'] += 1
                continue
            bc = s1[: a.bc_len]
            umi = s1[a.bc_len : a.bc_len + a.umi_len]
            if 'N' in bc or 'N' in umi:
                c['has_N'] += 1
                continue
            bcs[bc] += 1
            rid = h1[1:].split()[0]
            # SAMPLE NAME, not the observed BC6. These runs are already demultiplexed,
            # so the in-read barcode is redundant AND noisy: on HEK_2 the true barcode is
            # 86.1% of reads and a further 13.2% sit at Hamming distance 1 from it, i.e.
            # single sequencing errors (95 distinct 6-mers seen, of 4,096 possible).
            # Keying on the observed barcode would scatter ~13% of every sample's reads
            # into spurious wells and split their molecules apart.
            # Underscores stripped from the sample name: downstream shares drugseq's
            # rsplit('_', 2) parser, so 'HEK_2' would be read as well='2'. Harmless here
            # by luck (the numbers happen to be unique) but wrong, and it would break
            # silently on any other naming.
            out.write(f'@{rid}_{a.sample.replace("_","")}_{umi}\n{s2}\n+\n{q2}\n')
            c['written'] += 1
    dom, dn = bcs.most_common(1)[0]
    sys.stderr.write(
        f'  {a.sample}: reads {c["reads"]:,} | written {c["written"]:,} | '
        f'dominant BC {dom} = {100*dn/max(c["written"],1):.1f}% '
        f'| N-drop {c["has_N"]:,}\n'
    )


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Supplement (21 human wells): inputs for the three risdiplam 500 nM wells through the benchmark arm path.
  native: mate-1 records of the raw pair BAM (what the native arm BAM holds) -> <arm>/supp21/<well>/Aligned.out.bam
  50 nt:  primary mate-1 reads as FASTQ -> <arm 50nt>/supp21/reads.150.fq.gz, for the matched-depth tables
          (the per-sample supplement itself is native only)
Then (shell): build_read_table.sh per native BAM; align.sh + build_read_table.sh for the 50-nt reads; slices."""

from settings import WORK, PROJECT
import os, sys, gzip, pysam

P = WORK
PR = PROJECT
A = f'{P}/benchmark_uniform/bob57_24plex_pe_native/supp21'
A50 = f'{P}/benchmark_uniform/bob57_24plex_pe_hs/supp21'
WELLS = ['Ris_dose_500mM-1', 'Ris_dose_500mM-2', 'Ris_dose_500mM-3']
os.makedirs(A50, exist_ok=True)
fq = gzip.open(f'{A50}/reads.150.fq.gz', 'wt')
for w in WELLS:
    src = f'{PR}/data/bams/pairs_raw/{w}_pairs_raw.bam'
    os.makedirs(f'{A}/{w}', exist_ok=True)
    b = pysam.AlignmentFile(src, 'rb')
    out = pysam.AlignmentFile(f'{A}/{w}/Aligned.out.bam', 'wb', template=b)
    n1 = n2 = nf = 0
    for a in b.fetch(until_eof=True):
        if a.is_read2:
            n2 += 1
            continue
        # raw pair BAM names are '<id>__<umi>_<bobcode>_<x>_<y>'; the arm's are '<id>_<bobcode>_<umi>' (read_table.py: well, umi = rsplit('_', 2)[1:])
        rid, rest = a.query_name.split('__', 1)
        umi, bobcode = rest.split('_')[:2]
        a.query_name = f'{rid}_{bobcode}_{umi}'
        n1 += 1
        out.write(a)
        if a.is_secondary or a.is_supplementary:
            continue
        seq = a.get_forward_sequence()
        q = a.get_forward_qualities()
        if seq is None:
            continue
        fq.write(f'@{a.query_name}\n{seq}\n+\n{"".join(chr(x + 33) for x in q)}\n')
        nf += 1
    out.close()
    print(w, f'mate-1 records {n1:,} (mate-2 skipped {n2:,}); FASTQ reads {nf:,}', flush=True)
fq.close()
print('DONE')

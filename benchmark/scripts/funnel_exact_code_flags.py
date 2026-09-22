#!/usr/bin/env python3
"""exact_code_flags.py <R2 fastq.gz URI or path> <out.npz>
One pass over the raw read-2 FASTQ of the 24-plex run (read 2 begins with the bobcode at position 0). For every read the
first 7 nt (11 nt for the 11-mer BOB25D) are compared with the declared bobcodes: a read is EXACT when they equal one code with
zero mismatches. Writes the sorted uint64 keys of the exact reads (key = lane, tile, x, y of the Illumina read id packed into
one integer, so it joins the BAM read names without hashing) plus the tallies. The pipeline's prep step accepts one mismatch
and discards the count, which is why this pass exists (the preprint scores exact codes only, like the nanopore analysis)."""

from settings import RUN_JSON
import sys, json, subprocess, array, numpy as np

src, out = sys.argv[1:3]
RUN = json.load(open(RUN_JSON))
CODES7 = {c for c in RUN['tso_species_map'] if len(c) == 7}
CODES11 = {c for c in RUN['tso_species_map'] if len(c) == 11}


def key_of(rid):  # 'INSTRUMENT:RUN:FLOWCELL:7:1115:12911:4025' -> lane<<56 | tile<<40 | x<<20 | y
    p = rid.split(':')
    return (int(p[3]) << 56) | (int(p[4]) << 40) | (int(p[5]) << 20) | int(p[6])


def ham(a, b):
    return sum(x != y for x, y in zip(a, b))


cmd = ['bash', '-c', f"zcat '{src}'"]
p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
keys = array.array('Q')
n = n_exact = n_1mm = n_none = 0
it = p.stdout
for hdr in it:
    seq = next(it)
    next(it)
    next(it)
    n += 1
    s7 = seq[:7]
    s11 = seq[:11]
    if s11 in CODES11 or s7 in CODES7:
        n_exact += 1
        keys.append(key_of(hdr[1:].split(None, 1)[0]))
    else:
        d = min([ham(s7, c) for c in CODES7] + [ham(s11, c) for c in CODES11]) if CODES7 else 9
        if d == 1:
            n_1mm += 1
        else:
            n_none += 1
    if n == 200000:
        print(
            f'first {n:,} reads: exact {100*n_exact/n:.1f}%  1-mismatch {100*n_1mm/n:.1f}%  no code {100*n_none/n:.1f}%',
            flush=True,
        )
p.wait()
k = np.frombuffer(keys, dtype=np.uint64)
k.sort()
np.savez(out, keys=k, n=n, n_exact=n_exact, n_1mm=n_1mm, n_none=n_none)
print(
    f'reads {n:,}  exact {n_exact:,} ({100*n_exact/n:.2f}%)  one mismatch {n_1mm:,} ({100*n_1mm/n:.2f}%)  no code within 1 mismatch {n_none:,} ({100*n_none/n:.2f}%)  unique keys {len(np.unique(k)):,}',
    flush=True,
)

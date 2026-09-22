#!/usr/bin/env python3
"""The canonical per-read table. One row per mapped primary read, for every benchmark arm.

Every downstream figure is a filter on this table rather than its own pass over a BAM, so an
arm cannot silently acquire its own gene-assignment or multimapper policy. The region call is
recorded rather than applied: exon-only and exon+intron are then both available as a slice,
which matters because the arms differ enormously in intronic content (prime-seq ~47%,
BRB-seq ~25%) and exon-only quietly discards half of one arm's data.

Columns:  well  umi  contig  pos  nh  region  gene  strand

  pos     the read's 5'-end anchor (POS for +, last aligned base for -), NOT raw SAM POS
  strand  + or -

  nh      U unique (NH:i:1) | M multimapping
  region  E exonic protein-coding (gene named)   X exonic, other biotype or unnamed
          I intronic (inside a gene, no exon overlap; gene name still recorded)
          G intergenic     R rRNA/rDNA     T mitochondrial
  gene    gene name for E and I when the containing gene is protein-coding, else '.'

--mod N --rem K processes only records where index %% N == K, so N workers can split one BAM.
"""

from settings import PIPELINE_DIR
import os, subprocess, sys, bisect, argparse

sys.path.insert(0, PIPELINE_DIR)
import star_taxonomy_illumina as st  # noqa: E402
import umi_dedup_illumina as ud  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('bam')
ap.add_argument('--mod', type=int, default=1)
ap.add_argument('--rem', type=int, default=0)
ap.add_argument('--threads', type=int, default=2)
a = ap.parse_args()

genes, gst, exons, est = st.parse_gtf(os.environ['BOBSEQ_GTF'])
ridna, rdst = st.load_rdna(os.environ['BOBSEQ_RDNA_BED'])


def in_rdna(ch, s, e):
    v = ridna.get(ch)
    if not v:
        return False
    i = bisect.bisect_right(rdst[ch], e)
    for j in range(max(0, i - 1), -1, -1):
        lo, hi = v[j]
        if hi < s:
            break
        if lo < e and hi > s:
            return True
    return False


MOD, REM = a.mod, a.rem
n_rec = n_out = 0
out = sys.stdout
p = subprocess.Popen(
    ['samtools', 'view', '-@', str(a.threads), '-F', '0x904', a.bam], stdout=subprocess.PIPE, text=True
)
for line in p.stdout:
    n_rec += 1
    if MOD > 1 and (n_rec - 1) % MOD != REM:
        continue
    f = line.split('\t', 6)
    qname, ch, pos, cig = f[0], f[2], int(f[3]), f[5]
    parts = qname.rsplit('_', 2)
    if len(parts) != 3:
        continue
    well, umi = parts[1], parts[2]
    nh = 'U' if ('NH:i:1\t' in line or line.rstrip('\n').endswith('NH:i:1')) else 'M'
    s, e = pos - 1, ud._refend(pos, cig)
    # the molecule anchor is the read's 5' end: POS for a forward read, the last aligned
    # base for a reverse read, so a soft-clipped 3' end cannot move it
    rev = int(f[1]) & 16
    anchor, strand = (e, '-') if rev else (pos, '+')

    if ch.endswith('_MT'):
        region, gene = 'T', '.'
    elif in_rdna(ch, s, e):
        region, gene = 'R', '.'
    else:
        gs = st._overlap(genes, gst, ch, s, e, collect=True)
        if not gs:
            region, gene = 'G', '.'
        else:
            best = max(gs, key=lambda g: min(g[1], e) - max(g[0], s))
            named = best[3] if best[2] == 'protein_coding' and best[3] else None
            if not st._overlap(exons, est, ch, s, e):
                # intronic reads carry their containing gene too: the exon+intron slice
                # needs it, and prime-seq is ~47% intronic so dropping it is not neutral
                region, gene = 'I', (named or '.')
            elif named:
                region, gene = 'E', named
            else:
                region, gene = 'X', '.'
    n_out += 1
    out.write(f'{well}\t{umi}\t{ch}\t{anchor}\t{nh}\t{region}\t{gene}\t{strand}\n')
p.wait()
sys.stderr.write(f'  worker {REM}/{MOD}: {n_rec:,} scanned, {n_out:,} written\n')

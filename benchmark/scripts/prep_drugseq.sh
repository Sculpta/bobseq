#!/bin/bash
# Attach the DRUG-seq well barcode and UMI to the R2 read name, so a molecule survives
# alignment as (well, UMI) and can be deduplicated afterwards.
#
# DRUG-seq R1 is 20 nt: 10 nt well barcode then 10 nt UMI. R2 is 52 nt of cDNA and
# carries NO barcode, which is why this data cannot go through the BOBseq read
# preparation (prep_bobseq_illumina.py): that step calls a species bobcode positionally
# at R2[0:7], and here that position is cDNA (12,629 distinct 7-mers per 100k reads
# against ~2 in a real bobcode run). It would complete and report meaningless species
# metrics rather than fail.
#
# Output QNAME: <origid>_<well>_<umi>
# The '_' separator matches the BOBseq read names, so the same downstream
# splitting works and STAR preserves QNAME up to the first whitespace.
#
# Reads are dropped when the barcode or UMI contains N: an ambiguous base there is an
# unusable molecule identity, and keeping it would invent a distinct molecule per error.
#
# The benchmark runs this with --keep-dups: every read is kept and molecules are counted
# after alignment on (well, UMI, position) like every other arm (molecule_defs.sh), so a
# UMI shared by two molecules at different positions is not collapsed. Without the flag the
# script deduplicates in the same pass on (well, UMI) alone, keeping the first read seen,
# which aligns about 3.6x fewer reads on a full plate but merges such collisions.
#
# USAGE: prep_drugseq.sh <R1.fastq.gz> <R2.fastq.gz> <out_R2.fastq.gz> [--keep-dups]
set -euo pipefail
R1="$1"; R2="$2"; OUT="$3"
KEEPDUPS=0; [[ "${4:-}" == "--keep-dups" ]] && KEEPDUPS=1
BC=10; UMI=10

paste <(zcat "$R1" | paste - - - -) <(zcat "$R2" | paste - - - -) \
| "${AWK:-awk}" -v BC="$BC" -v UMI="$UMI" -v KEEPDUPS="$KEEPDUPS" -F'\t' '
    {
      n++
      r1seq = $2; r2id = $5; r2seq = $6; r2qual = $8
      if (length(r1seq) < BC + UMI) { short++; next }
      well = substr(r1seq, 1, BC)
      umi  = substr(r1seq, BC + 1, UMI)
      if (well ~ /N/ || umi ~ /N/) { ambig++; next }
      if (!KEEPDUPS) {
        key = well umi
        if (key in seen) { dup++; next }
        seen[key] = 1
      }
      split(r2id, a, " ")
      id = substr(a[1], 2)                       # drop the leading @
      print "@" id "_" well "_" umi "\n" r2seq "\n+\n" r2qual
      kept++
    }
    END {
      printf("PREP  total=%d kept=%d (%.2f%%) dropped_short=%d dropped_N=%d\n",
             n, kept, 100.0*kept/n, short+0, ambig+0) > "/dev/stderr"
      if (!KEEPDUPS) printf("PREP  deduplicated: %d duplicate reads dropped, %d molecules kept\n",
                            dup+0, kept) > "/dev/stderr"
    }' \
| ${GZIP_CMD:-gzip -1} > "$OUT"

echo "PREP_DONE $OUT"

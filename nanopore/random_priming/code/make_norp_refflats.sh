#!/bin/bash
# Build ribosomal-protein-EXCLUDED refFlats used by the gene-body coverage metric.
# Input: the combined human+mouse protein-coding refFlats (rRNA + MT already excluded),
# with HUMAN_/MOUSE_ contig prefixes. Ribosomal-protein genes are short and extremely
# abundant, so with short-insert selection they create a shoulder in the metagene; removing
# them isolates the coverage shape. The pattern is ANCHORED so it does NOT catch the RPS6K*
# ribosomal-S6 kinases (which are not ribosomal proteins).
#
# Usage: make_norp_refflats.sh <ref_annotation_dir> <out_dir>
#   <ref_annotation_dir> must contain:
#     combined_genome.HUMAN.proteincoding.nordna_nomt.refFlat
#     combined_genome.MOUSE.proteincoding.nordna_nomt.refFlat
set -euo pipefail
REF="${1:?ref annotation dir}"; OUT="${2:?out dir}"; mkdir -p "$OUT"
HREF="$REF/combined_genome.HUMAN.proteincoding.nordna_nomt.refFlat"
MREF="$REF/combined_genome.MOUSE.proteincoding.nordna_nomt.refFlat"
HPAT='^(RP[LS][0-9]+[A-Z]?|RPLP[0-9]|RPSA|MRP[LS][0-9]+[A-Z]?|FAU)$'
MPAT='^(Rp[ls][0-9]+[A-Za-z]?|Rplp[0-9]|Rpsa|Mrp[ls][0-9]+[A-Za-z]?|Fau)$'
awk -v p="$HPAT" '$1 !~ p' "$HREF" > "$OUT/HUMAN.noRP.refFlat"
awk -v p="$MPAT" '$1 !~ p' "$MREF" > "$OUT/MOUSE.noRP.refFlat"
echo "human RP removed: $(( $(wc -l < "$HREF") - $(wc -l < "$OUT/HUMAN.noRP.refFlat") ))  -> $OUT/HUMAN.noRP.refFlat"
echo "mouse RP removed: $(( $(wc -l < "$MREF") - $(wc -l < "$OUT/MOUSE.noRP.refFlat") ))  -> $OUT/MOUSE.noRP.refFlat"
echo "sanity — RPS6K* kinases kept: $(grep -cE '^RPS6K' "$OUT/HUMAN.noRP.refFlat")"

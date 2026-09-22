#!/bin/bash
# The competitor arms of the benchmark (DRUG-seq and prime-seq), from the public data fetched by fetch_competitors.sh.
#   build_arms.sh
#
# Every arm goes through the identical frozen path: prep (barcode and UMI into the read name) -> truncate to BM_READLEN ->
# STAR (settings.sh) -> read table -> molecule definitions -> composition -> rarefaction -> per-sample funnel.
#   drugseq          the whole DRUG-seq plate at 50 nt; the 24 DMSO wells are selected at the slice (well_map.tsv, wells_keep.txt)
#   drugseq_native   the reads of the 24 DMSO wells at their native 52 nt
#   primeseq         the eight prime-seq samples (50 nt is their native length)
#   primeseq_native  the same arm under the native name, so every method has a native set
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/scripts/settings.sh"
U=$BM_UNIFORM; DATA=$BM_WORK/data; mkdir -p "$U" "$BM_WORK/prepped"
log() { echo "[$(date '+%F %T')] $*"; }
composition() {
  local D=$1
  mawk -F'\t' '{r[$6]++; if($5=="U")u++} END{ t=0; for(k in r)t+=r[k]; printf "      unique(NH=1) %.1f%%\n", 100*u/t; split("E X I G R T", o, " "); for(i=1;i<=6;i++){ k=o[i]; printf "      %s %12d  %5.1f%%\n", k, r[k], 100*r[k]/t } }' "$D/reads.tsv" > "$D/composition.txt"
  bash "$HERE/scripts/composition_by_species.sh" "$D"
  mawk -F'\t' '$3 ~ /^HUMAN_/ {a[$6]++; if($5=="U") u[$6]++} END{for(k in a) print k"\t"a[k]"\t"(k in u ? u[k] : 0)}' "$D/reads.tsv" > "$D/composition_human.txt"
  bash "$HERE/scripts/composition_benchmark_wells.sh" "$D"
  bash "$HERE/scripts/composition_per_well.sh" "$D"
}
arm() {   # <name> <prepped reads.fq.gz> <truncate: yes|no>
  local NAME=$1 SRC=$2 TRUNC=$3 D=$U/$1
  mkdir -p "$D"
  log "[$NAME]"
  if [ "$TRUNC" = yes ]; then
    [ -s "$D/reads.$BM_READLEN.fq.gz" ] || bash "$HERE/scripts/truncate_fastq.sh" "$SRC" "$D/reads.$BM_READLEN.fq.gz"
    FQ=$D/reads.$BM_READLEN.fq.gz
  else
    FQ=$SRC
  fi
  [ -s "$D/Aligned.out.bam" ] || bash "$HERE/scripts/align.sh" "$FQ" "$D" ""
  grep -E "Number of input reads|Average input read length|Uniquely mapped reads %|multiple loci \|" "$D/Log.final.out" | sed 's/^ */      /'
  [ -s "$D/reads.tsv" ] || bash "$HERE/scripts/build_read_table.sh" "$D/Aligned.out.bam" "$D" "$BM_THREADS"
  bash "$HERE/scripts/molecule_defs.sh" "$D" E U 2>&1 | grep -E "mol_|SANITY|rror|wells" | head -14
  composition "$D"
  python3 "$HERE/scripts/rarefy.py" "$D" E_U 2>&1 | tail -1
  python3 "$HERE/scripts/funnel_per_sample.py" "$D" 2>&1 | grep -v "index file" | tail -1
}

# ---- DRUG-seq: well barcode (10 nt) and UMI (10 nt) of read 1 into the read name, duplicates kept ----
DS=$BM_WORK/prepped/drugseq_prepped.fq.gz
if [ ! -s "$DS" ]; then
  log "DRUG-seq prep"
  AWK=mawk GZIP_CMD="pigz -1 -p 6" bash "$HERE/scripts/prep_drugseq.sh" "$DATA/drugseq/SRR14730306_1.fastq.gz" "$DATA/drugseq/SRR14730306_2.fastq.gz" "$DS" --keep-dups
fi
python3 "$HERE/scripts/drugseq_well_map.py" "$U/drugseq"
arm drugseq "$DS" yes
# native: the reads of the 24 DMSO wells (observed barcodes within one substitution of a DMSO barcode) at 52 nt
DN=$U/drugseq_native; mkdir -p "$DN"; cp "$U/drugseq/well_map.tsv" "$U/drugseq/wells_keep.txt" "$DN/"
if [ ! -s "$DN/reads.native.fq.gz" ]; then
  mawk -F'\t' 'NR==FNR{k[$1]=1; next} ($2 in k){print $1}' "$DN/wells_keep.txt" "$DN/well_map.tsv" | sort -u > "$DN/keep_barcodes.txt"
  pigz -dc "$DS" | mawk -v kf="$DN/keep_barcodes.txt" 'BEGIN{while((getline l < kf)>0) keep[l]=1}
    NR%4==1{ n=split($1,f,"_"); k=(f[n-1] in keep) } k' | pigz -1 -p 4 > "$DN/reads.native.fq.gz"
fi
arm drugseq_native "$DN/reads.native.fq.gz" no

# ---- prime-seq: sample name and UMI (10 nt of read 1) into the read name, one file per sample, concatenated ----
PS=$BM_WORK/prepped/primeseq_prepped.fq.gz
if [ ! -s "$PS" ]; then
  log "prime-seq prep"
  mkdir -p "$BM_WORK/prepped/primeseq"
  for s in $(tail -n +2 "$BM_CONFIG/primeseq_samples.tsv" | cut -f1 | sort -u); do
    python3 "$HERE/scripts/prep_primeseq.py" --bcumi "$DATA/primeseq/${s}_R16.fastq.gz" --cdna "$DATA/primeseq/${s}_R50.fastq.gz" --sample "$s" --out "$BM_WORK/prepped/primeseq/$s.fq.gz"
  done
  cat $(ls "$BM_WORK/prepped/primeseq"/*.fq.gz | sort) > "$PS"
fi
arm primeseq "$PS" yes
[ -e "$U/primeseq_native" ] || ln -s primeseq "$U/primeseq_native"
log "ARMS DONE"

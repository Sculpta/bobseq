#!/bin/bash
# 5'->3' gene-body coverage for ONE coordinate-sorted combined-genome BAM, with the exact
# conventions used in the manuscript:
#   * inserts <= 500 bp only          (aligned cDNA length = sum of CIGAR M/I/=/X)
#   * ribosomal-protein mRNAs removed (HUMAN.noRP / MOUSE.noRP refFlats)
#   * MINIMUM_LENGTH = 1000           (only transcripts >= 1000 bp enter the metagene)
#   * STRAND_SPECIFICITY = NONE
# Two Picard runs per BAM (human refFlat, mouse refFlat) split species without splitting the BAM.
#
# Usage: picard_genebody.sh <sorted.bam> <HUMAN.noRP.refFlat> <MOUSE.noRP.refFlat> <out_prefix>
# Requires: samtools, picard (Picard 3.4.0 used).
set -euo pipefail
SB="${1:?sorted bam}"; HREF="${2:?human noRP refFlat}"; MREF="${3:?mouse noRP refFlat}"; PRE="${4:?out prefix}"
MAX=500; ML=1000
FB="${PRE}.ins${MAX}.bam"
samtools view -h -F 0x4 "$SB" \
 | awk -v M=$MAX 'BEGIN{OFS="\t"} /^@/{print;next}{c=$6;ql=0;n="";for(i=1;i<=length(c);i++){ch=substr(c,i,1);if(ch>="0"&&ch<="9")n=n ch;else{if(ch=="M"||ch=="I"||ch=="="||ch=="X")ql+=n;n=""}} if(ql<=M)print}' \
 | samtools view -b -o "$FB" -
samtools index "$FB"
picard CollectRnaSeqMetrics -I "$FB" -O "${PRE}.human.ins${MAX}_noRP_ml${ML}.rnaseq_metrics.txt" -REF_FLAT "$HREF" -STRAND_SPECIFICITY NONE -MINIMUM_LENGTH $ML
picard CollectRnaSeqMetrics -I "$FB" -O "${PRE}.mouse.ins${MAX}_noRP_ml${ML}.rnaseq_metrics.txt" -REF_FLAT "$MREF" -STRAND_SPECIFICITY NONE -MINIMUM_LENGTH $ML
# per-species primary-read counts that passed the <=500 filter (for figure legends / balance)
read h m < <(samtools view -F 0x904 "$FB" | awk '{if($3~/^HUMAN_/)h++;else if($3~/^MOUSE_/)m++}END{print (h+0)" "(m+0)}')
echo -e "human_reads\t$h\nmouse_reads\t$m" > "${PRE}.readcounts.tsv"
rm -f "$FB" "$FB.bai"
echo "wrote ${PRE}.{human,mouse}.ins${MAX}_noRP_ml${ML}.rnaseq_metrics.txt and .readcounts.tsv"

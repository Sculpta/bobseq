#!/bin/bash
# BOBseq Illumina processing: one multiplexed paired-end library (2 x 150 nt) to per-sample deduplicated alignments.
#
#   bash pipeline/run_illumina.sh <R1.fastq.gz> <R2.fastq.gz> <run.json> <reference_dir> <out_dir> [threads]
#
# Steps, in order:
#   1 guide       render the run configuration (bobcodes, species map, read layout) into the analysis guide the modules read
#   2 prep        call the bobcode at the start of read 2, move bobcode and UMI into the read name, trim the G-run and poly(A);
#                 every read pair of the library, one pass (prep_reads_illumina.py)
#   3 demux       split the prepped pairs by the bobcode in the read name into one pair per sample (demux_qname_illumina.py)
#   4 align       STAR per sample against the combined human + mouse reference, both mates; the mate-1 file holds read 2
#   5 dedup       remove PCR duplicates per sample: same bobcode, contig, strand, 5' position, fragment length and UMI
#                 (dedup_reads_illumina.py); originals are kept beside the deduplicated files as *.all
#   6 tables      per-sample read accounting, composition and the per-sample report (per_sample_*.py)
# Outputs under <out_dir>: prepped reads and prep_stats.json, demux/ (per-sample prepped pairs and demux stats),
# samples/<sample>/ (STAR logs, <sample>_star_pairs.bam and <sample>_star_combined.bam before and after deduplication,
# <sample>_dup_prefilter.json), <run>_per_sample.tsv, <run>_per_sample_composition.tsv, <run>_per_sample_report.pdf.
set -euo pipefail
R1=$1; R2=$2; RUN_JSON=$3; REF=$4; OUT=$5; THREADS=${6:-8}
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(dirname "$HERE"); SCRIPTS="$HERE/scripts"; RUN=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['run'])" "$RUN_JSON")
export PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg MPLCONFIGDIR="$OUT/.mpl"
export BOBSEQ_REF_DIR="$REF" BOBSEQ_GTF="$REF/combined_genome.gtf" BOBSEQ_RDNA_BED="$ROOT/config/rdna_loci.bed"
export BOBCODE_REF_DIR="$ROOT/config" BOBSEQ_NEB_INDEX_FILE="$ROOT/config/nebnext_e7600_dual_index_primers.tsv"
mkdir -p "$OUT/.mpl" "$OUT/samples"
log() { echo "[$(date '+%F %T')] $*"; }

log "1 guide"
python3 "$SCRIPTS/make_guide.py" "$RUN_JSON" -o "$OUT" > "$OUT/make_guide.log"
GUIDE="$OUT/machine_readable_guide/${RUN}_analysisguide.txt"; [ -s "$GUIDE" ] || { echo "guide not written"; exit 1; }

log "2 prep ($THREADS chunks in parallel)"
if [ ! -s "$OUT/prep_stats.json" ]; then
  # prep has no cross-read state, so the library is split into chunks of 4 M pairs, prepped in parallel and concatenated;
  # the per-chunk statistics are merged by prep_reads_illumina.py merge-stats (identical to a single pass)
  C="$OUT/prep_chunks"; mkdir -p "$C"; NL=$((4 * 4000000))
  pigz -dc "$R1" | split -l $NL -d -a 3 --filter='pigz -1 > $FILE.fq.gz' - "$C/R1."
  pigz -dc "$R2" | split -l $NL -d -a 3 --filter='pigz -1 > $FILE.fq.gz' - "$C/R2."
  ls "$C"/R1.*.fq.gz | sed -E 's/.*R1\.([0-9]+)\.fq\.gz/\1/' | xargs -P "$THREADS" -I{} bash -c "python3 '$SCRIPTS/prep_reads_illumina.py' '$C/R1.{}.fq.gz' '$C/R2.{}.fq.gz' '$C/prepped.{}.fastq.gz' --r1-out '$C/prepped_R1.{}.fastq.gz' --guide '$GUIDE' --stats '$C/prep_stats.{}.json' > '$C/prep.{}.log' 2>&1"
  n_chunks=$(ls "$C"/R1.*.fq.gz | wc -l); n_stats=$(ls "$C"/prep_stats.*.json | wc -l)
  [ "$n_chunks" = "$n_stats" ] || { echo "prep failed on $((n_chunks - n_stats)) chunk(s), see $C/prep.*.log"; exit 1; }
  cat $(ls "$C"/prepped.[0-9]*.fastq.gz | sort) > "$OUT/prepped.fastq.gz"
  cat $(ls "$C"/prepped_R1.[0-9]*.fastq.gz | sort) > "$OUT/prepped_R1.fastq.gz"
  python3 "$SCRIPTS/prep_reads_illumina.py" merge-stats --guide "$GUIDE" --out "$OUT/prep_stats.json" $(ls "$C"/prep_stats.[0-9]*.json | sort) > "$OUT/prep.log"
  cat "$C"/prep.[0-9]*.log >> "$OUT/prep.log"
  rm -f "$C"/R1.*.fq.gz "$C"/R2.*.fq.gz "$C"/prepped.[0-9]*.fastq.gz "$C"/prepped_R1.[0-9]*.fastq.gz
fi
python3 -c "import json,sys; s=json.load(open(sys.argv[1])); print(f\"  pairs {s['n']:,}, bobcode called {s['bobcode_called']:,}, written {s['written']:,}\")" "$OUT/prep_stats.json"

log "3 demux by bobcode"
if [ ! -s "$OUT/demux/library_demux_stats.json" ]; then
  python3 "$SCRIPTS/demux_qname_illumina.py" "$OUT/prepped.fastq.gz" - "$OUT/demux" --run-json "$RUN_JSON" --r1 "$OUT/prepped_R1.fastq.gz" > "$OUT/demux.log" 2>&1
fi
SAMPLES=$(python3 -c "import json,sys; print(' '.join(sorted(json.load(open(sys.argv[1]))['per_bobcode'][c]['sample'].replace(' ','_').replace('/','_') for c in json.load(open(sys.argv[1]))['per_bobcode'])))" "$OUT/demux/library_demux_stats.json")

log "4 align + 5 dedup, per sample"
for S in $SAMPLES; do
  D="$OUT/samples/$S"; mkdir -p "$D"
  PREP="$OUT/demux/samples/$S/prepped_${S}.fastq"; MATE="$OUT/demux/samples/$S/mateR1_${S}.fastq"
  [ -s "$PREP" ] || { log "  $S: no reads, skipped"; continue; }
  if [ ! -s "$D/${S}_star_pairs.bam" ]; then
    STAR --runThreadN "$THREADS" --genomeDir "$REF/STAR_index" --readFilesIn "$PREP" "$MATE" --outFileNamePrefix "$D/${S}_" \
         --outSAMtype BAM Unsorted --outSAMattributes NH HI AS nM --outSAMunmapped Within \
         --outSAMmultNmax 18446744073709551615 --outFilterMultimapNmax 10000 \
         --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 \
         --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.04 --alignIntronMax 1000000 --alignSJDBoverhangMin 3 --alignEndsType Local > "$D/star.log" 2>&1
    mv "$D/${S}_Aligned.out.bam" "$D/${S}_star_pairs.bam"
    samtools view -b -f 0x40 -o "$D/${S}_star_combined.bam" "$D/${S}_star_pairs.bam"   # mate 1 = read 2, the bobcode side: one record per fragment
    samtools quickcheck "$D/${S}_star_pairs.bam" "$D/${S}_star_combined.bam"
  fi
  # dedup_reads_illumina.py works on a run directory of fixed layout
  W="$D/dedup"; mkdir -p "$W/machine_readable_guide" "$W/fastq/$S" "$W/individual-analyses/$S"
  ln -sfn "$SCRIPTS" "$W/scripts"; cp "$GUIDE" "$W/machine_readable_guide/"
  if [ ! -s "$W/individual-analyses/$S/${S}_dup_prefilter.json" ]; then
    cp "$D/${S}_star_pairs.bam" "$W/individual-analyses/$S/${S}_star_pairs.bam"; cp "$D/${S}_star_combined.bam" "$W/individual-analyses/$S/${S}_star_combined.bam"
    cp "$PREP" "$W/fastq/$S/prepped_${S}.fastq"; cp "$MATE" "$W/fastq/$S/prepped_R1_${S}.fastq"
    python3 "$SCRIPTS/dedup_reads_illumina.py" "$W" "$S" > "$D/dedup.log" 2>&1
    cp "$W/individual-analyses/$S/${S}_dup_prefilter.json" "$D/${S}_dup_prefilter.json"
  fi
  log "  $S: $(python3 -c "import json,sys; s=json.load(open(sys.argv[1])); print(f\"{s['n_umi']:,} aligned fragments, {s['unique']:,} molecules kept ({s['dup_pct']}% duplicates)\")" "$D/${S}_dup_prefilter.json")"
done

log "6 per-sample tables and report"
BAMS=$(ls "$OUT"/samples/*/dedup/individual-analyses/*/*_star_combined.bam); SIDECARS=$(ls "$OUT"/samples/*/*_dup_prefilter.json)
python3 "$SCRIPTS/per_sample_table.py" --run-json "$RUN_JSON" --demux-stats "$OUT/demux/library_demux_stats.json" --bams $BAMS --rdna-bed "$BOBSEQ_RDNA_BED" --out "$OUT/${RUN}_per_sample.tsv"
PYTHONPATH="$SCRIPTS" python3 "$SCRIPTS/per_sample_composition.py" --bams $BAMS --run-json "$RUN_JSON" --per-sample "$OUT/${RUN}_per_sample.tsv" --out "$OUT/${RUN}_per_sample_composition.tsv"
python3 "$SCRIPTS/per_sample_report.py" --per-sample "$OUT/${RUN}_per_sample.tsv" --composition "$OUT/${RUN}_per_sample_composition.tsv" --demux-stats "$OUT/demux/library_demux_stats.json" \
  --prep-stats "$OUT/prep_stats.json" --sidecars $SIDECARS --run-json "$RUN_JSON" --assets "$ROOT/config" --run-name "$RUN" --out "$OUT/${RUN}_per_sample_report.pdf"
log "done: $OUT"

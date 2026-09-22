#!/bin/bash
# The per-sample BAM store (benchmark_uniform/per_sample_bams/): one coordinate-sorted, indexed BAM per sample and set, human
# contigs only with Ensembl names (1..22, X, Y, MT, scaffolds), primary mapped records, duplicates not removed. The coverage,
# Picard, junction and pooled-track steps read these.
#   build_per_sample_bams.sh
#
#   drugseq / primeseq / bobseq_50nt   the arms' 50-nt alignments split by sample (split_pooled.py)
#   drugseq_native                     the native 52-nt DRUG-seq alignments split by sample
#   primeseq_native                    = primeseq (50 nt is native)
#   bobseq_pe_native                   both mates of the pipeline's raw paired-end alignments, per well (human_pair_bam.py);
#                                      every human well of the 24-plex, so the Ris 500 nM wells are available to the supplement
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/scripts/settings.sh"
U=$BM_UNIFORM; B=$U/per_sample_bams; DATA=$BM_WORK/data; mkdir -p "$B"; cd "$B"
log() { echo "[$(date '+%F %T')] $*"; }
# ---- sample maps: key in the read name -> sample name ----
python3 - "$BM_CONFIG" "$BM_RUN_JSON" "$U" <<'PY'
import csv, json, re, sys
CONFIG, RJ, U = sys.argv[1:4]
def safe(n): return re.sub(r'_+', '_', ''.join(c if (c.isalnum() or c in '-_.') else '_' for c in n)).strip('_')
# DRUG-seq: observed barcode -> corrected barcode; the 24 DMSO wells, named by plate well
wm = dict(l.split() for l in open(f'{U}/drugseq/well_map.tsv'))
dmso = set(open(f'{CONFIG}/drugseq_dmso_wells.txt').read().split())
well_of = {r['well_index']: r['plate_well'] for r in csv.DictReader(open(f'{CONFIG}/drugseq_plate_VH02001704.csv'))}
with open('samples_drugseq.tsv', 'w') as f:
    for obs, corr in wm.items():
        if corr in dmso: f.write(f'{obs}\tDMSO_{well_of[corr]}_{corr}\n')
# prime-seq: the sample name is the key
samples = sorted({r['sample'].replace('_', '') for r in csv.DictReader(open(f'{CONFIG}/primeseq_samples.tsv'), delimiter='\t')})
with open('samples_primeseq.tsv', 'w') as f:
    for s in samples: f.write(f'{s}\t{s}\n')
# BOBseq: bobcode -> sample name, the 18 benchmark wells (the matched arm) and every human well (the native set)
cfg = json.load(open(RJ)); lab = {k.upper(): v for k, v in cfg['bobcode_labels'].items()}
keep = [l.strip() for l in open(f'{U}/bob57_24plex_pe_hs/wells_keep.txt') if l.strip()]
with open('samples_bobseq.tsv', 'w') as f:
    for c in keep: f.write(f'{c}\t{safe(lab[c])}\n')
with open('samples_bobseq_human.tsv', 'w') as f:
    for c, sp in cfg['tso_species_map'].items():
        if sp == 'human': f.write(f'{c.upper()}\t{safe(lab[c.upper()])}\n')
print('sample maps:', {k: sum(1 for _ in open(f'samples_{k}.tsv')) for k in ('drugseq', 'primeseq', 'bobseq', 'bobseq_human')})
PY
# ---- split the pooled arm BAMs ----
split_set() {   # <set dir> <pooled bam> <sample map>
  local m=$1 src=$2 map=$3
  [ -d "$m" ] && [ "$(ls "$m"/*.bam.bai 2>/dev/null | wc -l)" -gt 0 ] && { log "$m: present"; return 0; }
  log "split $m"; mkdir -p "$m"
  python3 "$HERE/scripts/split_pooled.py" "$src" "$B/$m" "$m" "$B/$map" >> records.tsv
  for f in "$m"/*.unsorted.bam; do s=$(basename "$f" .unsorted.bam); samtools sort -@ 4 -m 2G -o "$m/$s.bam" "$f" && samtools index "$m/$s.bam" && rm -f "$f"; done
  log "$m: $(ls "$m"/*.bam | wc -l) BAMs"
}
: > records.tsv
split_set primeseq "$U/primeseq/Aligned.out.bam" samples_primeseq.tsv
split_set bobseq_50nt "$U/bob57_24plex_pe_hs/Aligned.out.bam" samples_bobseq.tsv
split_set drugseq "$U/drugseq/Aligned.out.bam" samples_drugseq.tsv
split_set drugseq_native "$U/drugseq_native/Aligned.out.bam" samples_drugseq.tsv
[ -e primeseq_native ] || ln -s primeseq primeseq_native
# ---- BOBseq paired-end native: both mates, human contigs, every human well ----
mkdir -p bobseq_pe_native
cut -f2 samples_bobseq_human.tsv | xargs -P "$BM_THREADS" -I{} bash -c "[ -s '$B/bobseq_pe_native/{}.bam.bai' ] || python3 '$HERE/scripts/human_pair_bam.py' '$DATA/bams/pairs_raw/{}_pairs_raw.bam' '$B/bobseq_pe_native/{}.bam'"
log "bobseq_pe_native: $(ls bobseq_pe_native/*.bam | wc -l) BAMs"
python3 "$HERE/scripts/per_sample_metadata.py"
# the coverage sample sheets (label, BAM relative to this store, method) ship with the repository
mkdir -p "$BM_COVERAGE"; cp "$BM_CONFIG/coverage_samples_native.tsv" "$BM_COVERAGE/samples.tsv"; cp "$BM_CONFIG/coverage_samples_50nt.tsv" "$BM_COVERAGE/samples_50nt.tsv"
for f in "$BM_COVERAGE/samples.tsv" "$BM_COVERAGE/samples_50nt.tsv"; do
  while IFS=$'\t' read -r lab bam m; do [ -s "$B/$bam" ] || echo "missing BAM for $lab: $bam" >&2; done < "$f"
done
log "PER-SAMPLE BAMS DONE"

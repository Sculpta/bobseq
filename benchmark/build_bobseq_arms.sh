#!/bin/bash
# The two BOBseq arms of the benchmark, from the 24-plex lane and the pipeline's per-sample outputs.
#   build_bobseq_arms.sh <pipeline_out_dir> <R1.fastq.gz> <R2.fastq.gz>
#
#   <pipeline_out_dir>  the output directory of pipeline/run_illumina.sh for the 24-plex library (config/run_24plex.json)
#
# Arms (benchmark_uniform/):
#   bob57_24plex_pe_hs      read-length matched: the benchmark's own prep of the lane (prep_bobseq_illumina.py, the whole degenerate
#                           stretch after the bobcode as UMI), the 18 benchmark wells, read 2 truncated to 50 nt, STAR with the frozen
#                           parameters (settings.sh), read table with the UMI cut to 8 nt (6N + HH), then the common path
#   bob57_24plex_pe_native  native read length: mate-1 primary records of the pipeline's raw paired-end alignments (2 x 150 nt,
#                           STAR --outFilter*OverLread 0.33) of the same 18 wells, read names in the benchmark convention, then the
#                           same common path; reads with an N in the UMI are dropped from the read table
# Common path per arm: molecule_defs.sh E U (slice + molecule definitions), composition tables, rarefy.py, per-well input and
# mapping counts, the STAR log of the benchmark wells (make_selected_log.py), funnel_per_sample.py.
# Also stages the pipeline's per-sample BAMs and STAR logs under data/ in the layout the downstream scripts read.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/scripts/settings.sh"
PIPE=$1; R1=$2; R2=$3
U=$BM_UNIFORM; HS=$U/bob57_24plex_pe_hs; NAT=$U/bob57_24plex_pe_native; DATA=$BM_WORK/data
mkdir -p "$HS" "$NAT" "$DATA/bams/pairs_raw" "$DATA/bams/pipeline_dedup_pe/dedup_pairs_6NHH" "$DATA/bams/pipeline_dedup_pe/dedup_analysis_6NHH" "$DATA/pipeline_stats/star_logs"
log() { echo "[$(date '+%F %T')] $*"; }
composition() {   # <arm dir>: the composition tables every arm carries
  local D=$1
  mawk -F'\t' '{r[$6]++; if($5=="U")u++} END{ t=0; for(k in r)t+=r[k]; printf "      unique(NH=1) %.1f%%\n", 100*u/t; split("E X I G R T", o, " "); for(i=1;i<=6;i++){ k=o[i]; printf "      %s %12d  %5.1f%%\n", k, r[k], 100*r[k]/t } }' "$D/reads.tsv" > "$D/composition.txt"
  bash "$HERE/scripts/composition_by_species.sh" "$D"
  mawk -F'\t' '$3 ~ /^HUMAN_/ {a[$6]++; if($5=="U") u[$6]++} END{for(k in a) print k"\t"a[k]"\t"(k in u ? u[k] : 0)}' "$D/reads.tsv" > "$D/composition_human.txt"
  bash "$HERE/scripts/composition_benchmark_wells.sh" "$D"
  bash "$HERE/scripts/composition_per_well.sh" "$D"
}
common_path() {   # <arm dir>: definitions, composition, rarefaction, per-well tables, selected STAR log, per-sample funnel
  local D=$1
  log "[$(basename "$D")] molecule definitions"; bash "$HERE/scripts/molecule_defs.sh" "$D" E U 2>&1 | grep -E "mol_|SANITY|rror|wells" | head -14
  log "[$(basename "$D")] composition"; composition "$D"
  log "[$(basename "$D")] rarefaction"; python3 "$HERE/scripts/rarefy.py" "$D" E_U 2>&1 | tail -1
  python3 "$HERE/scripts/make_selected_log.py" "$D" "$D/input_counts.tsv" "$D/perwell_mapped.tsv"
  log "[$(basename "$D")] per-sample funnel"; python3 "$HERE/scripts/funnel_per_sample.py" "$D" 2>&1 | grep -v "index file" | tail -1
}

# ---- sample sheet of the run: species per bobcode, the excluded wells, the 18 benchmark wells ----
SAMPLES=$(python3 -c "import json,sys; r=json.load(open(sys.argv[1])); print(' '.join(sorted(r['bobcode_labels'][c].replace(' ','_').replace('/','_') for c in r['bobcode_labels'])))" "$BM_RUN_JSON")
python3 - "$BM_RUN_JSON" "$BM_CONFIG/units_excluded_24plex.tsv" "$HS" "$NAT" <<'PY'
import json, sys
rj, ex_table, hs, nat = sys.argv[1:5]
cfg = json.load(open(rj))
excluded = [l.split('\t')[0] for l in open(ex_table).read().splitlines()[1:] if l.strip()]
for D in (hs, nat):
    open(f'{D}/well_species.tsv', 'w').writelines(f"{c.upper()}\t^{sp.upper()}_\n" for c, sp in cfg['tso_species_map'].items())
    open(f'{D}/wells_keep.txt', 'w').writelines(f'{c.upper()}\n' for c in cfg['tso_species_map'] if c.upper() not in excluded)
    open(f'{D}/units_excluded.tsv', 'w').write(open(ex_table).read())
    open(f'{D}/arm.env', 'w').write("SPECIES='^(HUMAN|MOUSE)_'\nUMI_MAXED=0\nUMI_TAG_SPACE=\"8:36864,7:9216\"\n")
print('benchmark wells:', sum(1 for c in cfg['tso_species_map'] if c.upper() not in excluded), '| excluded:', excluded)
PY
EX=$(awk -F'\t' 'NR>1{printf "%s ", $1}' "$HS/units_excluded.tsv")

# ---- stage the pipeline outputs: raw pair BAMs, deduplicated BAMs, STAR logs ----
log "staging pipeline outputs from $PIPE"
for S in $SAMPLES; do
  D=$PIPE/samples/$S
  [ -s "$D/${S}_star_pairs.bam" ] || { echo "missing $D/${S}_star_pairs.bam" >&2; exit 1; }
  ln -sfn "$D/${S}_star_pairs.bam" "$DATA/bams/pairs_raw/${S}_pairs_raw.bam"
  ln -sfn "$D/dedup/individual-analyses/$S/${S}_star_pairs.bam" "$DATA/bams/pipeline_dedup_pe/dedup_pairs_6NHH/${S}_star_pairs.bam"
  ln -sfn "$D/dedup/individual-analyses/$S/${S}_star_combined.bam" "$DATA/bams/pipeline_dedup_pe/dedup_analysis_6NHH/${S}_star_combined.bam"
  cp "$D/${S}_Log.final.out" "$DATA/pipeline_stats/star_logs/${S}_Log.final.out"
done
cp "$PIPE/prep_stats.json" "$DATA/pipeline_stats/prep_stats.json"
# the exclusion rule, checked on this run: uniquely mapped reads per human well below 25% of the median over the human wells
python3 - "$BM_RUN_JSON" "$DATA/pipeline_stats/star_logs" "$HS/units_excluded.tsv" <<'PY'
import json, sys, statistics
rj, logs, ex_table = sys.argv[1:4]
cfg = json.load(open(rj)); lab = {c.upper(): l.replace(' ', '_').replace('/', '_') for c, l in cfg['bobcode_labels'].items()}
uniq = {}
for c, sp in cfg['tso_species_map'].items():
    if sp != 'human': continue
    kv = {x.split('|')[0].strip(): x.split('|')[1].strip() for x in open(f'{logs}/{lab[c.upper()]}_Log.final.out') if '|' in x}
    uniq[c.upper()] = int(kv['Uniquely mapped reads number'])
med = statistics.median(uniq.values()); low = {c for c, n in uniq.items() if n < 0.25 * med}
listed = {l.split('\t')[0] for l in open(ex_table).read().splitlines()[1:] if 'input failure' in l}
print('input-failure rule on this run (unique reads < 25% of the human-well median):', sorted(lab[c] for c in low))
assert low == listed, f'wells below the rule {sorted(low)} differ from the recorded exclusion table {sorted(listed)}'
PY

# ---- matched arm: the benchmark's own prep of the lane, 18 wells, 50 nt ----
if [ ! -s "$HS/prepped.fq.gz" ]; then
  log "[bob57_24plex_pe_hs] prep of the lane ($BM_THREADS chunks in parallel)"
  C=$HS/chunks; mkdir -p "$C"; NL=$(( 4 * 12000000 ))
  pigz -dc "$R1" | split -l $NL -d -a 3 --filter='pigz -1 > $FILE.fq.gz' - "$C/R1."
  pigz -dc "$R2" | split -l $NL -d -a 3 --filter='pigz -1 > $FILE.fq.gz' - "$C/R2."
  ls "$C"/R1.*.fq.gz | sed -E 's/.*R1\.([0-9]+)\.fq\.gz/\1/' | xargs -P "$BM_THREADS" -I{} bash -c "python3 '$HERE/scripts/prep_bobseq_illumina.py' '$C/R1.{}.fq.gz' '$C/R2.{}.fq.gz' '$C/prepped.{}.fq.gz' --run-json '$BM_RUN_JSON' --umi full --stats '$C/prep_stats.{}.json' 2> '$C/prep.{}.log'"
  python3 - "$C" "$HS" "$EX" <<'PY'
import json, glob, sys, collections
C, HS, ex = sys.argv[1], sys.argv[2], sys.argv[3].split()
tot = collections.Counter(); mode = None
for f in sorted(glob.glob(f'{C}/prep_stats.[0-9]*.json')):
    for k, v in json.load(open(f)).items():
        if isinstance(v, int): tot[k] += v
        else: mode = v
out = dict(tot); out['umi_mode'] = mode
out['written'] = sum(v for k, v in tot.items() if k.startswith('well:') and k[5:] not in ex)
for e in ex: out.pop(f'well:{e}', None)
out['units_excluded'] = ex
json.dump(out, open(f'{HS}/prep_stats.json', 'w'), indent=1)
print('  prep: ' + ' | '.join(f'{k}={v:,}' for k, v in out.items() if isinstance(v, int) and not k.startswith('well:')))
PY
  # the 18 benchmark wells, in chunk order; the excluded wells are dropped by the bobcode in the read name
  cat $(ls "$C"/prepped.[0-9]*.fq.gz | sort) | pigz -dc | mawk -v ex="$EX" 'BEGIN{n=split(ex,a," "); for(i=1;i<=n;i++) drop[a[i]]=1}
    NR%4==1{ split($1,f,"_"); keep=!(f[length(f)-1] in drop) } keep' | pigz -1 -p 4 > "$HS/prepped.fq.gz"
  rm -rf "$C"/R1.*.fq.gz "$C"/R2.*.fq.gz "$C"/prepped.[0-9]*.fq.gz
fi
log "[bob57_24plex_pe_hs] truncate to $BM_READLEN nt"; bash "$HERE/scripts/truncate_fastq.sh" "$HS/prepped.fq.gz" "$HS/reads.$BM_READLEN.fq.gz"
log "[bob57_24plex_pe_hs] STAR"; bash "$HERE/scripts/align.sh" "$HS/reads.$BM_READLEN.fq.gz" "$HS" ""
grep -E "Number of input reads|Uniquely mapped reads %|multiple loci \|" "$HS/Log.final.out" | sed 's/^ */      /'
log "[bob57_24plex_pe_hs] read table"; bash "$HERE/scripts/build_read_table.sh" "$HS/Aligned.out.bam" "$HS" "$BM_THREADS"
# the UMI of the benchmark is the first 8 nt after the bobcode (6N + HH); the prep kept the whole degenerate stretch
mawk -F'\t' -v OFS='\t' '{ $2 = substr($2, 1, 8); print }' "$HS/reads.tsv" > "$HS/reads.umi8.tsv" && mv "$HS/reads.umi8.tsv" "$HS/reads.tsv"
python3 -c "import json,sys; st=json.load(open(sys.argv[1])); open(sys.argv[2],'w').write(''.join(f'{k[5:]}\t{v}\n' for k,v in st.items() if k.startswith('well:')))" "$HS/prep_stats.json" "$HS/input_counts.tsv"
mawk -F'\t' -v OFS='\t' '{m[$1]++; if($5=="U") u[$1]++} END{for(w in m) print w"\t"m[w]"\t"(w in u?u[w]:0)}' "$HS/reads.tsv" > "$HS/perwell_mapped.tsv"
common_path "$HS"

# ---- native arm: mate-1 records of the pipeline's raw paired-end alignments of the 18 wells ----
log "[bob57_24plex_pe_native] mate-1 records of the raw pair BAMs, names in the benchmark convention"
mkdir -p "$NAT/parts"
for S in $SAMPLES; do
  code=$(python3 -c "import json,sys; r=json.load(open(sys.argv[1])); print(next(c.upper() for c,l in r['bobcode_labels'].items() if l.replace(' ','_').replace('/','_')==sys.argv[2]))" "$BM_RUN_JSON" "$S")
  grep -qx "$code" "$NAT/wells_keep.txt" || continue
  echo "$S"
done | xargs -P "$BM_THREADS" -I{} bash -c "[ -s '$NAT/parts/{}.bam' ] || python3 '$HERE/scripts/native_arm_mate1.py' '$DATA/bams/pairs_raw/{}_pairs_raw.bam' '$NAT/parts/{}.bam'"
samtools cat -o "$NAT/Aligned.out.bam" "$NAT"/parts/*.bam
echo "native: paired-end 2x150, mate-1 (read 2) primary records only in Aligned.out.bam (one record per fragment), raw pre-deduplication alignments of the pipeline" > "$NAT/README.native.txt"
python3 - "$NAT" "$DATA/pipeline_stats/star_logs" "$NAT/parts" <<'PY'
# a library-level STAR log for the arm: the sums over the per-sample logs of the pipeline (reads = read pairs)
import sys, os, glob
D, logs, parts = sys.argv[1:4]
keep = [os.path.basename(p)[:-4] for p in sorted(glob.glob(f'{parts}/*.bam'))]
KEYS = ['Number of input reads', 'Uniquely mapped reads number', 'Number of reads mapped to multiple loci', 'Number of reads mapped to too many loci', 'Number of reads unmapped: too short', 'Number of reads unmapped: other']
tot = {k: 0 for k in KEYS}; lens = []
for s in keep:
    for l in open(f'{logs}/{s}_Log.final.out'):
        if '|' not in l: continue
        k, v = [x.strip() for x in l.split('|', 1)]
        if k in KEYS: tot[k] += int(v)
        if k == 'Average input read length': lens.append(float(v))
n = tot['Number of input reads']
with open(f'{D}/Log.final.out', 'w') as f:
    f.write(f"                          Number of input reads |\t{n}\n                      Average input read length |\t{int(round(sum(lens) / len(lens)))}\n")
    for k, pk in (('Uniquely mapped reads number', 'Uniquely mapped reads %'), ('Number of reads mapped to multiple loci', '% of reads mapped to multiple loci'), ('Number of reads mapped to too many loci', '% of reads mapped to too many loci'), ('Number of reads unmapped: too short', '% of reads unmapped: too short'), ('Number of reads unmapped: other', '% of reads unmapped: other')):
        f.write(f"{k:>47s} |\t{tot[k]}\n{pk:>47s} |\t{100 * tot[k] / n:.2f}%\n")
    f.write(f"# summed over the {len(keep)} per-sample STAR logs of the paired-end pipeline run; 'reads' = read pairs\n")
print('Log.final.out assembled from', len(keep), 'samples; input pairs', n)
PY
log "[bob57_24plex_pe_native] read table"; bash "$HERE/scripts/build_read_table.sh" "$NAT/Aligned.out.bam" "$NAT" "$BM_THREADS"
# reads with an N in the UMI leave the table; counted per well so the input counts can be corrected
mawk -F'\t' -v OFS='\t' -v D="$NAT" '$2 ~ /N/ {d[$1]++; next} {print; m[$1]++; if($5=="U") u[$1]++} END{for(w in m) print w"\t"m[w]"\t"(w in u?u[w]:0) > (D "/perwell_mapped.tsv"); for(w in d) print w"\t"d[w] > (D "/n_umi_dropped_per_well.tsv")}' "$NAT/reads.tsv" > "$NAT/reads.noN.tsv"
mv "$NAT/reads.noN.tsv" "$NAT/reads.tsv"; touch "$NAT/n_umi_dropped_per_well.tsv"
python3 - "$BM_RUN_JSON" "$DATA/pipeline_stats/star_logs" "$NAT" <<'PY'
# input reads per well = the pipeline's STAR input pairs of the sample minus the pairs dropped for an N in the UMI. Every well of the
# run is listed: the 18 benchmark wells after the drop, the other wells (risdiplam 500 nM, mouse) as sequenced, for the supplements
import json, sys, os
rj, logs, D = sys.argv[1:4]
cfg = json.load(open(rj)); lab = {c.upper(): l.replace(' ', '_').replace('/', '_') for c, l in cfg['bobcode_labels'].items()}
drop = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(f'{D}/n_umi_dropped_per_well.tsv') if l.strip()}
keep = [l.strip() for l in open(f'{D}/wells_keep.txt') if l.strip()]
with open(f'{D}/input_counts.tsv', 'w') as out:
    for code in keep + [c for c in lab if c not in keep]:
        kv = {x.split('|')[0].strip(): x.split('|')[1].strip() for x in open(f'{logs}/{lab[code]}_Log.final.out') if '|' in x}
        out.write(f"{code}\t{int(kv['Number of input reads']) - drop.get(code, 0)}\n")
print('input counts for', len(lab), 'wells')
PY
common_path "$NAT"
log "BOBSEQ ARMS DONE"

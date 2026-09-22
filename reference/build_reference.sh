#!/bin/bash
# Build the combined human + mouse reference used throughout: Ensembl release 113 primary assemblies with contigs prefixed
# HUMAN_ / MOUSE_, the combined GTF, the STAR index (sjdbOverhang 99, genomeSAindexNbases 14) and the Picard refFlat.
# The combined genome is 5.8 Gb; STAR needs about 60 GB of RAM to index it. Idempotent: finished stages are skipped.
#   bash reference/build_reference.sh [WORKDIR]      (default: reference/combined_GRCh38.113_GRCm39.113, where benchmark/settings.sh expects it)
set -euo pipefail
WORK=${1:-$(dirname "$0")/combined_GRCh38.113_GRCm39.113}; THREADS=${THREADS:-16}; BASE=https://ftp.ensembl.org/pub/release-113; HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$WORK/dl"; cd "$WORK"
step() { echo "=== [$(date -u +%H:%M:%S)] $* ==="; }
command -v STAR >/dev/null || { echo "STAR not found (environment.yml provides 2.7.11b)"; exit 1; }
command -v samtools >/dev/null || { echo "samtools not found"; exit 1; }
STAR --version

step "downloading Ensembl 113 primary assemblies and GTFs (unmasked)"
cd dl
for u in $BASE/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
         $BASE/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz \
         $BASE/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz \
         $BASE/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz; do
  f=$(basename "$u"); [ -s "$f" ] || wget -q -O "$f" "$u"
done
for f in *.gz; do gzip -t "$f" || { echo "corrupt download: $f"; exit 1; }; done
cd ..

if [ ! -s combined_genome.fa.fai ]; then
  step "combined_genome.fa"
  zcat dl/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | sed 's/^>/>HUMAN_/' >  combined_genome.fa
  zcat dl/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz | sed 's/^>/>MOUSE_/' >> combined_genome.fa
  samtools faidx combined_genome.fa
fi
echo "contigs: $(wc -l < combined_genome.fa.fai) (HUMAN_ $(grep -c '^HUMAN_' combined_genome.fa.fai), MOUSE_ $(grep -c '^MOUSE_' combined_genome.fa.fai))"

# FS and OFS must both be tab: with the default FS, awk would rewrite the spaces inside the attribute column into tabs,
# STAR would still index the file, and every downstream attribute regex would silently match nothing.
if [ ! -s combined_genome.gtf ]; then
  step "combined_genome.gtf"
  zcat dl/Homo_sapiens.GRCh38.113.gtf.gz | awk 'BEGIN{FS=OFS="\t"} !/^#/{$1="HUMAN_"$1} 1' >  combined_genome.gtf
  zcat dl/Mus_musculus.GRCm39.113.gtf.gz | awk 'BEGIN{FS=OFS="\t"} !/^#/{$1="MOUSE_"$1} 1' >> combined_genome.gtf
fi
step "validating the GTF"
awk -F'\t' '!/^#/{print $1}' combined_genome.gtf | sort -u > gtf_chroms.txt; cut -f1 combined_genome.fa.fai | sort -u > fa_chroms.txt
[ "$(comm -23 gtf_chroms.txt fa_chroms.txt | wc -l)" -eq 0 ] || { echo "GTF contigs missing from the FASTA:"; comm -23 gtf_chroms.txt fa_chroms.txt | head; exit 1; }
[ "$(awk -F'\t' '!/^#/ && NF!=9 {n++} END{print n+0}' combined_genome.gtf)" -eq 0 ] || { echo "GTF lines without 9 fields: attribute column corrupt"; exit 1; }
for tag in 'gene_biotype "' 'gene_name "' 'transcript_id "'; do grep -q "$tag" combined_genome.gtf || { echo "GTF lacks attribute $tag"; exit 1; }; done

if [ ! -s STAR_index/SA ]; then
  step "STAR genomeGenerate (about 45 minutes)"
  mkdir -p STAR_index
  STAR --runMode genomeGenerate --genomeDir STAR_index --genomeFastaFiles combined_genome.fa --sjdbGTFfile combined_genome.gtf \
       --sjdbOverhang 99 --genomeSAindexNbases 14 --runThreadN "$THREADS" --limitGenomeGenerateRAM 120000000000 --outFileNamePrefix STAR_index/
fi
step "Picard refFlat"
[ -s combined_genome.refflat ] || python3 "$HERE/make_refflat.py" --gtf combined_genome.gtf -o combined_genome.refflat
cp "$HERE/../config/rdna_loci.bed" .
step "done: $WORK holds combined_genome.fa(.fai), combined_genome.gtf, combined_genome.refflat, rdna_loci.bed, STAR_index/"

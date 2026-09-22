#!/bin/bash
# Download the public competitor data used by the benchmark into the working directory.
#   fetch_competitors.sh
#
#   DRUG-seq   GSE176150 / SRR14730306 (Li et al. 2022): one 384-well U-2 OS plate, read 1 = 10 nt well barcode + 10 nt UMI,
#              read 2 = 52 nt cDNA. The 24 DMSO wells are the benchmark samples (config/drugseq_dmso_wells.txt).
#   prime-seq  E-MTAB-10142 (Janjic et al. 2022): the eight "Incubation + ProtK" 10,000-cell HEK293T samples
#              (config/primeseq_samples.tsv), read 1 (16 nt) = 6 nt barcode + 10 nt UMI, read 2 = 50 nt cDNA.
# Files are fetched from the ENA mirror and left under $BOBSEQ_WORK/data/{drugseq,primeseq}/; existing intact files are kept.
set -euo pipefail
source "$(dirname "$0")/settings.sh"
D=$BM_WORK/data; mkdir -p "$D/drugseq" "$D/primeseq"
get() {   # url  out
  local url=$1 out=$2
  if [ -s "$out" ] && gzip -t "$out" 2>/dev/null; then echo "  kept $out"; return 0; fi
  curl -sS --retry 5 --retry-delay 30 -o "$out.part" "$url" && gzip -t "$out.part" && mv "$out.part" "$out" && echo "  fetched $out" \
    || { rm -f "$out.part"; echo "FAILED $url" >&2; return 1; }
}
echo "DRUG-seq SRR14730306"
for m in 1 2; do
  get "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/006/SRR14730306/SRR14730306_${m}.fastq.gz" "$D/drugseq/SRR14730306_${m}.fastq.gz"
done
echo "prime-seq E-MTAB-10142, eight samples"
tail -n +2 "$BM_CONFIG/primeseq_samples.tsv" | while IFS=$'\t' read -r s L acc url; do
  get "https://$url" "$D/primeseq/${s}_R${L}.fastq.gz"
done
echo "FETCH DONE"

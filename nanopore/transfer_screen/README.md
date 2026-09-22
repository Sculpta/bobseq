# Species-mixing bobcode barcode accuracy (transfer mechanism screen)

Part of the bobseq repository that accompanies the preprint (citation, contact and licence in the top-level README).
Preprint panels: Fig. 3, S2, S3 and Table S1.

**Data:** NCBI BioProject PRJNA1532586, one SRA run per dataset; `fastq_manifest.tsv` lists the 186 datasets with read
count, MD5 and run accession (section 3).

Calculates barcode accuracy for every species-mixing dataset (splint, polydT, 5'TSO
and 3'TSO bobcodes) from raw Oxford Nanopore reads. Each read's species barcode is
compared with the species its cDNA aligns to.

## Files

| File | What it is |
|---|---|
| `speciesmix_bobcode_metrics.py` | the program (Python ≥ 3.8, standard library only) |
| `make_accuracy_figure.py` | draws the barcode accuracy figure from the program's output (needs matplotlib) |
| `make_supplementary_table.py` | writes Supplementary Table S1 from the program's output (needs openpyxl for .xlsx) |
| `tools/build_fastq_repository.py` | provenance: how the deposited FASTQ files were assembled from the MinKNOW output (not needed to reproduce the analysis) |
| `bobcode_general_rules.md` | the rules applied to every dataset; the program follows them exactly |
| `bobcode_run_rules.tsv` | everything that differs between datasets: chemistry, construct, codes, species, barcodes to score |
| `fastq_manifest.tsv` | the 186 deposited datasets: file, run, barcode, plot number, read count, MD5, source files and SRA run accession |
| `rdna_loci.bed` | rDNA regions on the combined genome (used to call rRNA reads); made for this project, so it cannot be downloaded |
| `reference_contigs.tsv` | order, length and MD5 of every contig in the combined genome, to check a rebuilt reference |

## 1. Software

- **Python** 3.8 or later. The main program needs no packages; the figure and table scripts need
  `pip install matplotlib openpyxl`.
- **samtools** on the `PATH` (1.20 was used).
- **STARlong 2.7.11b.** Use the official precompiled binary from
  <https://github.com/alexdobin/STAR/releases/tag/2.7.11b> (`bin/MacOSX_x86_64/STARlong` or
  `bin/Linux_x86_64_static/STARlong`). On Apple-silicon Macs it runs under Rosetta.
  Homebrew and bioconda builds of STAR were found to exit without error while reading
  0 reads on this machine, so use the official binary.
- **Memory:** STARlong loads the whole index (~57 GB), so the machine needs ~64 GB of RAM
  or more. The program runs one STARlong at a time and refuses to start a second.

## 2. Reference files

Everything comes from Ensembl. Chromosome names get a `HUMAN_` or `MOUSE_` prefix so that
one index holds both species.

### 2a. Download

```bash
mkdir -p refs && cd refs
# genome sequences (primary assembly, unmasked)
curl -O https://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
curl -O https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
# gene annotation, Ensembl release 113 for both species
curl -O https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
curl -O https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz
```

Human sequence: the original analysis used the Ensembl release 99 primary-assembly file,
so that is the file given here. The mouse file used carries no release number; any
Ensembl GRCm39 primary-assembly release is expected to give the same sequence. Step 2d
checks both, contig by contig.

### 2b. Combined genome (`combined_genome.fa`)

```bash
gunzip -c Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | sed -E 's/^>([^ ]+).*/>HUMAN_\1/' >  unordered.fa
gunzip -c Mus_musculus.GRCm39.dna.primary_assembly.fa.gz | sed -E 's/^>([^ ]+).*/>MOUSE_\1/' >> unordered.fa
samtools faidx unordered.fa
tail -n +2 ../reference_contigs.tsv | cut -f2 > contig_order.txt     # the index's contig order
samtools faidx -r contig_order.txt unordered.fa > combined_genome.fa
rm unordered.fa unordered.fa.fai
```

### 2c. Combined annotation (`combined_genome.gtf`)

Human lines first, then mouse, header (`#`) lines removed, chromosome column prefixed:

```bash
{ gunzip -c Homo_sapiens.GRCh38.113.gtf.gz | grep -v '^#' | sed 's/^/HUMAN_/'
  gunzip -c Mus_musculus.GRCm39.113.gtf.gz | grep -v '^#' | sed 's/^/MOUSE_/'; } > combined_genome.gtf
```

This is exactly the file used: 4,114,450 human + 2,475,479 mouse lines, identical in content
and order to the two Ensembl 113 GTFs.

### 2d. Check the genome

Every contig's sequence must match `reference_contigs.tsv` (255 contigs: 194 human, 61 mouse):

```bash
python3 - combined_genome.fa ../reference_contigs.tsv <<'EOF'
import sys, hashlib
exp = {l.split('\t')[1]: l.split('\t')[3].strip() for l in open(sys.argv[2]).readlines()[1:]}
got, name, h = {}, None, None
for line in open(sys.argv[1], 'rb'):
    if line[:1] == b'>':
        if name: got[name] = h.hexdigest()
        name, h = line[1:].split()[0].decode(), hashlib.md5()
    else:
        h.update(line.rstrip(b'\n').upper())
got[name] = h.hexdigest()
bad = [k for k in exp if got.get(k) != exp[k]]
print('OK: all', len(exp), 'contigs match' if not bad and len(got) == len(exp) else f'MISMATCH: {bad[:5]}')
EOF
```

### 2e. STAR index

```bash
STAR --runMode genomeGenerate --genomeDir STAR_index \
     --genomeFastaFiles combined_genome.fa --sjdbGTFfile combined_genome.gtf \
     --sjdbOverhang 99 --runThreadN 8 --limitGenomeGenerateRAM 90000000000
```

This is the command used (STAR 2.7.11b, default dense suffix array). Building it needs about
90 GB of RAM, and the finished index is ~57 GB. Do not use a sparse index
(`--genomeSAsparseD` > 1): it drops ~1% of reads on noisy libraries, so counts would not match.

## 3. Sequencing reads

The FASTQ files are deposited as one file per dataset, arranged as

```
<chemistry>_<construct>/<YYMMDD>_BC<NN>.fastq.gz      e.g. polydT_bob-polydt/260507_BC12.fastq.gz
```

These are the original MinKNOW reads, untrimmed and unfiltered. A dataset sequenced in
several chunks or MinKNOW runs has them joined into one file; one dataset (260120) is
limited to its first 10,000 reads per barcode. `fastq_manifest.tsv` in this directory lists
every dataset with its plot number (the "figure dataset" number in the SRA run titles), read
count, MD5 and SRA run accession; the reads themselves are in the SRA under BioProject
PRJNA1532586 (`data/sra_runs.tsv` in the repository root). The SRA file names carry the
directory in the name, `<chemistry>_<construct>__<YYMMDD>_BC<NN>.fastq.gz`: split at the
double underscore to obtain the layout above (the apostrophe of 3'TSO and 5'TSO is dropped
in both).

## 4. Running

All 186 datasets (one STARlong at a time; expect many hours):

```bash
python3 speciesmix_bobcode_metrics.py \
    --fastq-dir /path/to/fastqs --out-dir results \
    --starlong /path/to/STARlong --star-index refs/STAR_index --gtf refs/combined_genome.gtf
```

One chemistry, run or barcode: add `--chemistry polydT`, `--run 260522_BC9-17_polydtbob` and/or
`--barcode BC12` (`--exclude-run` leaves a run out). Instead of `--fastq-dir`, `--data-dir` reads
the original run folders named in the run rules; a few runs there point into the MinKNOW output
folder, given as `$MINKNOW_DATA` (set that environment variable). Finished
barcodes are skipped on a rerun; `--force` redoes them. `rdna_loci.bed`, `bobcode_run_rules.tsv`
and `bobcode_general_rules.md` are found next to the program automatically.

## 5. Figure and supplementary table

```bash
python3 make_accuracy_figure.py --results results/bobcode_metrics.tsv --out-dir figures
python3 make_supplementary_table.py --results-dir results \
    --index figures/barcode_accuracy_index.tsv --out-dir figures
```

The figure script writes `barcode_accuracy_all_chemistries_log.pdf/.png` and
`barcode_accuracy_index.tsv`. The index holds every plotted point and its number: points
are ranked by accuracy and numbered from 1 within each chemistry. The table script reads
that index, the per-dataset `metrics.json` files and the run rules. Sample conditions come
from the run file's `sample_notes` column and per-barcode flags from `barcode_flags`. It
writes `Supplementary_Table_S1.xlsx` and `.tsv`, which contain the datasets, column
definitions, exclusions and method.

## 6. Output

- `results/bobcode_metrics.tsv`, one row per dataset:
  - read categories: mRNA, rRNA, MT, other;
  - human and mouse mRNA reads;
  - barcode accuracy on human mRNA, on mouse mRNA and on all mRNA, each with the number of
    reads scored and a 95% Wilson interval for the total.

  `bobcode_general_rules.md` §9 defines every column.
- `results/<run>/<barcode>/metrics.json`: the count at every filter step, the code × species
  table, and provenance (checksums of the program and both rules files, software versions).

## Licence

MIT (`LICENSE` in the repository root). Copyright (c) 2026 Sculpta, Inc.

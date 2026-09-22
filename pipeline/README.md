# The BOBseq Illumina pipeline

`run_illumina.sh` processes one multiplexed BOBseq library (paired-end, 2 x 150 nt) into per-sample de-duplicated alignments
and a per-sample report. The steps are the Python modules in `scripts/`; the driver runs them in order on one machine.

```
bash pipeline/run_illumina.sh <R1.fastq.gz> <R2.fastq.gz> config/run_24plex.json <reference_dir> <out_dir> [threads]
```

The driver takes one R1/R2 pair, the sequenced lane. The SRA holds the lane as 24 per-sample pairs (BioProject PRJNA1532586): concatenate the 24 R1 files into one R1 and the 24 R2 files into one R2, in the same sample order, and run the driver on that pair. The deposited files are untrimmed and unmodified; the 4% of pairs without an assignable bobcode were not deposited.

Inside the container: `docker run --rm -v $PWD:/work -v <reference_dir>:/ref -w /work bobseq:0.1 bash pipeline/run_illumina.sh ... /ref out/ 16`.

| Step | Module | What it does |
|---|---|---|
| 1 guide | `make_guide.py` | renders the run configuration (bobcodes, sample labels, species map, read layout, G-run offsets) into the analysis guide the modules read |
| 2 prep | `prep_reads_illumina.py` | calls the bobcode at the start of read 2 (one substitution allowed, ties dropped), moves bobcode and UMI (6N + spacer) into the read name, measures the non-templated G-run, trims poly(A) read-through, writes the prepped pair; runs on 4 M-pair chunks in parallel, statistics merged with `merge-stats` |
| 3 demux | `demux_qname_illumina.py` | splits the prepped pairs by the bobcode in the read name into one pair of FASTQ files per sample |
| 4 align | STAR 2.7.11b | per sample, both mates, combined GRCh38.113 + GRCm39.113 reference, all multimappers kept (`--outFilterMultimapNmax 10000`, `--outFilterScoreMinOverLread 0.33`, `--outFilterMatchNminOverLread 0.33`) |
| 5 dedup | `dedup_reads_illumina.py`, `umi_dedup_illumina.py` | one molecule per (bobcode, contig, strand, 5' position, fragment length, UMI); the minimum read id is kept so the choice does not depend on thread order; originals kept as `*.all` |
| 6 tables | `per_sample_table.py`, `per_sample_composition.py`, `per_sample_report.py` | per-sample read accounting (assigned, mapped, species, rDNA), read composition by the classification rules of `star_taxonomy_illumina.py`, and the PDF report |

Outputs under `<out_dir>`: `prepped.fastq.gz`, `prepped_R1.fastq.gz`, `prep_stats.json`; `demux/` (per-sample prepped pairs,
`library_demux_stats.json`); `samples/<sample>/` (STAR logs, `<sample>_star_pairs.bam` and `<sample>_star_combined.bam` raw
and de-duplicated, `<sample>_dup_prefilter.json`); `<run>_per_sample.tsv`, `<run>_per_sample_composition.tsv`,
`<run>_per_sample_report.pdf`. The benchmark (`benchmark/build_bobseq_arms.sh`) reads `samples/` directly.

The other modules in `scripts/` are the library-level analysis the report draws on (`analyze_speciesmix_illumina.py` and
its sections, `filter_funnel_illumina.py`, `full_mol_funnel_illumina.py`, the chemistry overrides, `validate_outputs.py`,
`detect_species_map.py`); they are imported by the steps above and by the benchmark's read table
(`benchmark/scripts/read_table.py` imports `star_taxonomy_illumina.py` and `umi_dedup_illumina.py`, so the benchmark classifies
reads with the pipeline's own rules).

Test: `bash tests/run_smoke_test.sh` (prep and demultiplexing of the bundled 200,000-pair example, counts checked);
`bash tests/run_smoke_test.sh --full <reference_dir>` runs every step on it.

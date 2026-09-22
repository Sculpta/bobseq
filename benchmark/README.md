# The benchmark

Every number and panel of the preprint's benchmark (BOBseq vs DRUG-seq vs prime-seq) comes out of this directory.
`run_benchmark.sh <pipeline_out_dir> <R1.fastq.gz> <R2.fastq.gz>` runs the stage drivers below in order; each driver can be
run alone and skips outputs that exist. The drivers sit here; every script they call is in `scripts/`, where
`settings.sh` / `settings.py` hold every location and every frozen parameter (STAR parameters, read length of the matched
set, molecule definition); nothing is hard-coded elsewhere.

| Stage | Script | What it produces |
|---|---|---|
| 1 | `fetch_competitors.sh` | the public DRUG-seq and prime-seq reads (`data/`) |
| 2 | `build_arms.sh` | the competitor arms: prep (`prep_drugseq.sh`, `prep_primeseq.py`, `drugseq_well_map.py`), `truncate_fastq.sh`, `align.sh`, `build_read_table.sh` (`read_table.py`), `molecule_defs.sh` (`slice.sh`, `umi_b.py`, `dedup.c`), composition tables, `rarefy.py`, `funnel_per_sample.py` |
| 3 | `build_bobseq_arms.sh` | the BOBseq arms from the lane and the pipeline output: `prep_bobseq_illumina.py` for the read-length matched arm, `native_arm_mate1.py` for the native arm, then the same common path; stages the pipeline's per-sample BAMs and STAR logs under `data/` |
| 4 | `build_per_sample_bams.sh` | one human BAM per sample and set (`split_pooled.py`, `human_pair_bam.py`, `per_sample_metadata.py`) and the coverage sample sheets |
| 5 | `extract_junctions.sh` | per-sample junction tables: `dedup_umi_position.py`, MAPQ 255, `junctions_by_fragment.py` |
| 6 | `run_coverage.sh` | read positions on canonical transcripts (`coverage_lib.py`, `build_positions.py`), `junctions_per_sample.py`, `library_metrics.py`, the Picard refFlats (`build_refflat_filtered.py`), the Ris 500 nM wells (`supp21_*.py`) and the mouse wells (`mouse3_*.py`) through the same path |
| 7 | `run_funnel.sh` | the strict barcode-accuracy funnel per well (`funnel_*.py`) |
| 8 | `run_panels.sh` | the panels and their values tables (`composition_rules.py`, `picard_profile.py`, `sample_table_panel.py`, `read_fate_panels.py`, `main_figure_panels.py`, `composition_panels.py`, `threshold_panels.py`, `library_panels.py`, `panel_legends.py`, `junction_rarefaction*.py`, `supplement_*.py`, `stacked_supplement_merged.py`, `ucsc_pooled_tracks.py`, `gene_coverage_tracks.py`, the prime-seq peak check), `tidy_panels.py`, `panel_values.py`, `panel_readme.py`, `panel_sources.py`, `organize_figures.py` (the curated figure tree of the working directory) and the validators (`validate_*.py`, `check_genes_noRP.py`) |

`BM_COVSET=native|50nt` selects the figure set in the panel scripts (native = each method at its own read length, the
preprint's figures; 50nt = every method reduced to one 50-nt read, only its poly(A)-distance panel is used).
`BM_SPECIES_PREFIX=MOUSE_` switches `coverage_lib.py`, `composition_rules.py` and `library_metrics.py` to the mouse half
of the reference (the mouse wells).

## The working directory (`BOBSEQ_WORK`, default `<repo>/work`)

```
data/                       fetched competitor reads; staged pipeline BAMs and STAR logs; junctions/ (per-sample BEDs and the sample sheet)
benchmark_uniform/<arm>/    one directory per arm: prepped reads, Aligned.out.bam, reads.tsv (the read table), stats_E_U/ (slice,
                            molecule definitions, rarefied.tsv), composition tables, per-well counts, funnel_per_sample.tsv
benchmark_uniform/per_sample_bams/   one BAM per sample and set, metadata.tsv
coverage_architecture/      samples.tsv, positions_canonical*/, junctions*/, read_length_*.json, aligned_bases_*.tsv, the refFlats
results/                    library_metrics/ (insert per well), barcode_accuracy_funnel/, supplement_21/ (the Ris 500 nM and mouse wells)
preprint_figures_<set>/     generator output per figure set: main_figure_panels/ (panels, values/, README, VALIDATION), supplement_*, gene_coverage_tracks/
figures/                    the curated figure tree written by organize_figures.py (main/, supplement/, legends/, and the panels dropped from the figures)
```

Arms: `drugseq`, `primeseq`, `bob57_24plex_pe_hs` (the read-length matched set: 50 nt, one read per fragment) and
`drugseq_native`, `primeseq_native`, `bob57_24plex_pe_native` (each method at its sequenced read length). The read table
has one row per mapped primary read: well, UMI, contig, 5' position, unique/multimapping, region class (E exonic
protein-coding, X exonic other biotype, I intronic, G intergenic, R rRNA, T mitochondrial), gene, strand. A filtered read is
a unique read of class E on a human contig; a molecule is a unique (well, UMI, gene) after the over-use filter and the
collision correction of `umi_b_lib.py` (the UCI of the panels). `figures/METHODS_panels.md` describes every panel.

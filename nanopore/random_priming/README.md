# Random priming (3′TSO bobcode)

Everything needed to recreate the random-priming (3′TSO bobcode) figures of the preprint (Fig. 4):

1. **RACK1 genome-browser view** of two samples over all isoforms
2. **Picard gene-body coverage plot** (poly(dT) vs random priming)
3. **5′→3′ coverage-balance bar/box plot**
4. **Barcode-accuracy plot**

Plus the barcode-accuracy output table. The raw reads of the nine random-priming samples are in the NCBI SRA (BioProject
PRJNA1532586; the `randomPriming_bob-v1` rows of `data/sra_runs.tsv` in the repository root).

The two samples in the RACK1 panel are:
- **260730 bc12** — poly(dT)-primed 3′TSO bobcode (run `260730_RT_conditions_BC11-13`, barcode12) — *comparison*
- **260827 BC15** — random-priming 3′TSO bobcode (run `260827_BC14_16`, barcode15; 1× poly(A)) — *this study*

## Layout

```
README.md              this file
METHODS.md             manuscript methods section
code/                  analysis and plotting code
  genome_browser_gene.py          RACK1/ACTB/GAPDH two-sample browser (names the two samples)
  make_norp_refflats.sh           build ribosomal-protein-removed refFlats
  picard_genebody.sh              per-BAM gene-body coverage (inserts ≤500, RP removed, MINIMUM_LENGTH=1000)
  plot_genebody_profile.py        the Picard coverage plot (poly(dT) vs random, mean + 95% CI)
  plot_balance_box.py             the 5′→3′ balance box plot
  barcode_accuracy/
    bobcode_run_rules_randompriming.tsv   run rules for the 9 samples (the only per-run input)
    run_accuracy.sh                       driver: align + score + draw the figure
    make_accuracy_figure_random.py        the barcode-accuracy figure (paper style)
raw_fastq/             not shipped: the raw ONT reads, downloaded from the SRA for step 4 as
  randomPriming_bob-v1/<YYMMDD>_BC<NN>.fastq.gz   (9 datasets; the SRA file name with the double underscore replaced by a slash)
tables/                barcode-accuracy output
  barcode_accuracy_metrics_random_priming.tsv     full per-sample metrics (pipeline output)
  barcode_accuracy_index_random_priming.tsv       plotted points + 95% Wilson CI
data_derived/         inputs derived from the reads, so the figures are recreatable here
  genebody_metrics/random_priming/*.rnaseq_metrics.txt   Picard histograms for the 9 samples
  genebody_metrics/random_priming_counts.tsv             per-sample read counts (≤500 bp)
  browser_bams/*.rack1_actb_gapdh.bam                     region-sliced BAMs for the two browser samples
figures_reference/    the rendered figures (PNG + vector PDF/SVG) for reference
```

## Prerequisites

- Python 3 with `numpy`, `matplotlib` (and `scipy` for t-based CIs); `samtools` on `PATH`.
- For barcode accuracy (aligning from raw): the species-mixing barcode-accuracy program
  `speciesmix_bobcode_metrics.py` (in `../transfer_screen/`), **STARlong 2.7.11b**, and the
  combined `HUMAN_`/`MOUSE_`-prefixed GRCh38+GRCm39 **STAR index** and **combined_genome.gtf**
  (Ensembl 113). See `../transfer_screen/README.md` for building the reference. On Apple-silicon run STARlong
  under Rosetta (`DYLD_FALLBACK_LIBRARY_PATH=…/star_x86/lib`).

## Recreate each figure

### 1) RACK1 genome-browser view (two samples, all isoforms)
```bash
python3 code/genome_browser_gene.py --gene RACK1 --gtf /path/to/combined_genome.gtf
# also: --gene ACTB  /  --gene GAPDH   -> figures_reference/<gene>_two_samples_isoforms.{png,pdf,svg}
```
Uses the region-sliced BAMs in `data_derived/browser_bams/` (no full BAMs or alignment needed).

### 2) Picard gene-body coverage plot  &  3) 5′→3′ balance box plot
The Picard histograms for the nine random-priming samples are in
`data_derived/genebody_metrics/random_priming/`. To reproduce the **comparison** against poly(dT),
drop the poly(dT) metric files into `data_derived/genebody_metrics/polydt/` (same naming) and a
`polydt_counts.tsv`; these come from the main species-mixing deposit. Then:
```bash
python3 code/plot_genebody_profile.py     # -> figures_reference/genebody_profile_polydt_vs_random.*
python3 code/plot_balance_box.py          # -> figures_reference/balance_box_random_vs_polydt.*
```
Without the poly(dT) folder both scripts still run and plot the random-priming side alone.

To regenerate the metric files from a sorted BAM (e.g. after re-alignment):
```bash
bash code/make_norp_refflats.sh <ref_annotation_dir> <refflat_out_dir>
bash code/picard_genebody.sh <sample>.sorted.bam <refflat_out_dir>/HUMAN.noRP.refFlat \
     <refflat_out_dir>/MOUSE.noRP.refFlat data_derived/genebody_metrics/random_priming/<sample>
```

### 4) Barcode-accuracy table + figure
Download the nine SRA runs into `raw_fastq/randomPriming_bob-v1/<YYMMDD>_BC<NN>.fastq.gz` (layout above), then:
```bash
export PIPE=$PWD/../transfer_screen
export STARLONG=/path/to/STARlong  STAR_INDEX=/path/to/STAR_index  GTF=/path/to/combined_genome.gtf
bash code/barcode_accuracy/run_accuracy.sh
# -> results_barcode_accuracy/bobcode_metrics.tsv  (== tables/barcode_accuracy_metrics_random_priming.tsv)
# -> figures_reference/barcode_accuracy_all_chemistries_log.{png,pdf}
```
The pre-computed table is in `tables/` and the rendered figure in `figures_reference/`.

## Key result

All nine random-priming samples reach **98.3–100 % barcode accuracy** (median ~99 %), matching the
poly(dT)/3′TSO libraries, while giving **even gene-body coverage** (5′→3′ balance ≈ 1) in contrast
to poly(dT)'s strong 3′ bias — illustrated at single-gene resolution for RACK1 (and ACTB, GAPDH).

Contact: neal@sculpta.bio

# Within-method replicate reproducibility: BOBseq vs prime-seq vs DRUG-seq (Fig. 6D, Fig. S6A)

Code to reproduce the two-panel figures comparing replicate-to-replicate reproducibility of
three bulk RNA-seq library methods on the human mRNA transcriptome. One script produces both the
**Pearson-r** figure and the **RMSE** figure via `--metric`.

- **P1** — native sequencing depth; genes with ≥ 5 counts in both samples of a pair.
- **P2** — all libraries down-sampled to the smallest library (mRNA-assigned counts); metric
  computed over each pair's top-N most-expressed co-detected genes.

The reproducibility metric is per-replicate-pair, over an mRNA gene set, on log2(CPM+1):

- `--metric pearson` — **Pearson correlation** (higher = more reproducible).
- `--metric rmse` — **RMSE from the identity line y = x** (lower = more reproducible).

For BOBseq, only within-condition replicate pairs are used (triplicates), pooled across conditions;
prime-seq and DRUG-seq are each a single condition. Significance is a two-sided **Mann–Whitney U**
on the **per-sample mean** of the metric (the ~independent unit, since replicate pairs are not
independent). Note that the two metrics agree at native depth but can diverge at matched shallow
depth (e.g. down-sampled, BOBseq's correlation is significantly higher than prime-seq's while its
RMSE is comparable) — report both.

## Datasets

BOBseq libraries were generated in this study. The comparator datasets are public
(see `datasets.csv` for per-run accessions):

| Method | Cell line | Repository | Accession | Samples used | Reference |
|---|---|---|---|---|---|
| **prime-seq** | HEK293T | ArrayExpress | **E-MTAB-10142** | 8 untreated samples, "Incubation + Proteinase K", 10k cells (runs ERR5375207/241/365/212/248/284/320/356) | Janjic et al. 2022, *Genome Biol* 23:88 |
| **DRUG-seq** | U-2 OS | GEO | **GSE176150** (GSM5357049 / run SRR14730306, plate `VH02001704_S4`) | 24 DMSO vehicle wells, 24 h | Li et al. 2022, *ACS Chem Biol* 17(6):1401 — assay: Ye et al. 2018, *Nat Commun* 9:4307 |
| **BOBseq** | HEK | — | this study | 24-plex control + treatment triplicates | this study |

> Note: DRUG-seq here is the **2022** dataset (Li et al., GSE176150), not the original 2018 method paper.

## Inputs

1. **Gene-count matrix** (`--counts`), CSV or CSV.GZ: rows = genes, first two columns `GENEID,SYMBOL`,
   remaining columns one per sample. Sample columns are prefixed by method set
   (`bobseq_pe_native__`, `primeseq_native__`, `drugseq_native__`; edit the prefixes at the top of
   the script if yours differ). Counts are integer gene-level counts.
   The matrix of the preprint is the supplementary data table `S11_gene_counts_featurecounts.csv.gz`
   (78,932 Ensembl 113 genes x 119 sample-arms), produced by uniformly re-processing all raw reads: **STAR** (GRCh38,
   Ensembl 113) → UMI + position deduplication (`benchmark/scripts/dedup_umi_position.py`) → MAPQ 255 →
   **featureCounts** 2.1.1 on the Ensembl 113 annotation (`-s 1`, primary alignments, pairs counted once).
2. **Gene annotation** (`--annot`): either the Ensembl **GTF** (`.gtf`/`.gtf.gz`, e.g.
   `Homo_sapiens.GRCh38.113.gtf.gz` from https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/)
   or a TSV with columns `gene_id, biotype, chromosome, gene_name`.

The **mRNA gene set** is derived from the annotation: `protein_coding` genes, excluding
ribosomal-protein genes (symbols `RPL*`,`RPS*`,`MRPL*`,`MRPS*`, keeping `RPS6K*` kinases) and
mitochondrial genes (chromosome `MT`). This yields 19,934 genes for GRCh38/Ensembl 113.

## Usage

```bash
# Pearson-r figure (Fig. 6D)
python reproduce_reproducibility_figure.py \
  --counts <tables>/S11_gene_counts_featurecounts.csv.gz \
  --annot  Homo_sapiens.GRCh38.113.gtf.gz \
  --metric pearson --out figure_pearson

# RMSE figure (Fig. S6A)
python reproduce_reproducibility_figure.py \
  --counts <tables>/S11_gene_counts_featurecounts.csv.gz \
  --annot  Homo_sapiens.GRCh38.113.gtf.gz \
  --metric rmse --out figure_rmse

# ...or both at once
./run_all.sh <tables>/S11_gene_counts_featurecounts.csv.gz Homo_sapiens.GRCh38.113.gtf.gz
```

Each run writes `<out>.svg` (Illustrator-editable text), `.pdf`, `.png`, and a `.tsv` with
per-group summary statistics and the pairwise Mann–Whitney U p-values. The two `.tsv` files of the preprint's run
(`figure_pearson.tsv`, `figure_rmse.tsv`) are kept here as the numbers behind Fig. 6D and Fig. S6A. Requires Python 3
with numpy, scipy and matplotlib.

## Notes on choices

- **Why exclude ribosomal-protein and mitochondrial genes.** These few, ultra-high-expression genes
  are captured very differently across chemistries (e.g. ~19% of prime-seq and ~29% of DRUG-seq
  gene-assigned reads, vs <1% MT / ~15% RP for BOBseq) and would otherwise dominate the correlation.
- **CPM is computed within the mRNA set** (re-normalized after gene filtering).
- **P2 top-N.** `P2_TOP_N` (default 4168) is the largest per-pair top-N for which every replicate
  pair has that many genes detected in both samples at the full-set depth floor (48,837 counts);
  above this, correlations would be inflated/deflated by dropout. Adjust in the script if your data
  differ.
- **QC.** A BOBseq risdiplam-500 triplicate that failed input QC is excluded; all 24 DMSO DRUG-seq
  wells (including the systematically shallower plate-column-24 wells) are retained.

## Citation

If you use this code, please cite the BOBseq manuscript and the source datasets above
(Janjic et al. 2022; Li et al. 2022; Ye et al. 2018), and the tools STAR (Dobin et al. 2013),
featureCounts (Liao et al. 2014), and SciPy (Virtanen et al. 2020).

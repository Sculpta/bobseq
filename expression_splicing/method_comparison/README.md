# Method comparison: BOBseq controls against prime-seq, BRB-seq and TruSeq (Fig. 6A, 6B, S5A-C)

**Question.** Where do BOBseq's untreated HEK293T controls sit relative to the same cells profiled by three other library
methods, and which genes differ between the methods? Twenty-nine samples on the mRNA gene set: BOBseq's three control wells,
every control sample of prime-seq (8), BRB-seq (8) and TruSeq (10, three laboratories), every sample thinned to the smallest
library, the two full-length libraries (BOBseq, TruSeq) length-normalised.

| method | n | samples | source | unit |
|---|---|---|---|---|
| BOBseq | 3 | `Hek_control_1/2/3` of set `bobseq_pe_native` (random-primed, 2 x 150 nt, 6-nt UMI) | this study | TPM |
| prime-seq | 8 | HEK2, 3, 9, 20, 31, 42, 53, 64 (`primeseq_native`) | E-MTAB-10142 | CPM |
| BRB-seq | 8 | 2 untransfected + 6 empty-vector wells (`brbseq_native`) | GSE334309, batch 5 | CPM |
| TruSeq | 10 | 3 + 3 + 4 wild-type runs (`results/truseq_samples.tsv`) | GSE296724, GSE298560, GSE181978 | TPM |

## Method

1. **Counts.** Salmon on the Ensembl 113 cDNA index, summarised to genes with the tximport rule (count = sum of the
   transcripts' NumReads, TPM = sum of the transcripts' TPM, effective length = TPM-weighted mean of the transcripts'
   effective lengths); the mRNA gene set (S10, 19,934 genes) with alternative-haplotype copies folded onto their gene. Stage
   `tables` takes these from S7 / S8 / S9; stages `counts` + `mrna` compute them from Salmon `quant.sf` files.
2. **Thinning.** Every sample subsampled without replacement to the smallest mRNA library (seed 0), so that the 100-fold depth
   range does not drive gene selection or correlations. Expression = count / effective length for BOBseq and TruSeq (TPM),
   count for prime-seq and BRB-seq (CPM), per million over the gene set; counts are thinned first, divided second.
3. **PCA** (`src/pca.py`; Fig. 6A = `figures/pca_thinned_labs.svg`). log2(x + 1), genes with mean expression >= 1, the 5,000
   most variable, centred and scaled to unit variance, PCA with full SVD; PC1/PC2 and PC3/PC4, one marker per TruSeq laboratory.
   A clustered Spearman heatmap and its dendrogram on the same genes. Native and thinned runs.
4. **DESeq2** (`src/de.py`, `src/deseq2_contrast.R`; R). Six pairwise contrasts on the thinned counts, design `~ condition`,
   Wald test, apeglm shrinkage, genes with >= 10 counts, DE = padj < 0.001. For a contrast with BOBseq or TruSeq the gene
   length enters as DESeq2 normalisation factors (one length per gene per method: the mean effective length over the method's
   samples), so the fold change compares TPM with CPM. Run twice: on the full gene set and on the capped common set, the 7,500
   genes with the highest mean log2(CPM + 1) in their least-expressing method (a gene qualifies through the method that
   expresses it least; DESeq2 re-fitted on the subset).
5. **Signed Venn** (`src/de_capped_views.py`; Fig. 6B = `figures/de_capped/venn_vs_truseq_small.svg`). The capped-set DE genes of
   the three contrasts against TruSeq as (gene, direction) pairs: a region shared by two or three methods holds only genes DE in
   the same direction in each; up (higher than TruSeq) over down in every region; the direction-discordant genes counted below.
6. **Scatter** (`src/scatter.py`; Fig. S5A-C = `figures/scatter_truseq_vs_methods.svg`). Method means of log2(x + 1) on the thinned
   matrix, TruSeq against each other method, genes with mean expression >= 5 in both; outlier = 4-fold beyond the diagonal shifted
   by the median offset; the 15 per side with the largest deviation x mean expression named.

## Results of the preprint's run (`results/`)

PCA: native PC1 to PC4 = 41 / 22 / 14 / 4 %, thinned 29 / 17 / 10 / 4 %; the four methods form four groups, the TruSeq laboratories
separate on PC4 only (`pca_scores_thinned.tsv`, `pca_variance_thinned.tsv`, `genes_thinned.tsv`).

DESeq2 against TruSeq, capped common set (`de_capped/de_summary.tsv`, `de_capped_gene_set.tsv`): BOBseq 3,089 DE genes (1,591
up, 1,498 down), prime-seq 3,311 (1,695 / 1,616), BRB-seq 2,894 (1,342 / 1,552). Signed Venn (`de_capped/venn_vs_truseq.tsv`):
all three methods up 154 / down 319; BOBseq and prime-seq 147 / 102, BOBseq and BRB-seq 240 / 444, prime-seq and BRB-seq
497 / 378; BOBseq only 1,050 / 633, prime-seq only 897 / 817, BRB-seq only 451 / 411. Of the genes DE in both BOBseq and
prime-seq (1,499), 777 have opposite directions; BOBseq and BRB-seq 312 of 1,469; prime-seq and BRB-seq 227 of 1,575.

Scatter, TruSeq as the focal method (`scatter_truseq_summary.tsv`, `_labelled.tsv`, `_outliers.tsv`, `_outlier_overlap.tsv`):
TruSeq vs BOBseq 7,943 genes plotted, r 0.697, RMSE 1.51; vs prime-seq 9,735, r 0.779, RMSE 1.15; vs BRB-seq 9,664, r 0.796,
RMSE 1.10; 14 genes over- and 18 under-represented in TruSeq against all three other methods.

`results/sample_table.tsv` lists the 29 samples with the Salmon version, library type, fragment-length distribution and mapped
fraction of their quantification; `results/truseq_samples.tsv` the public runs.

## Run

```
python3 run.py --list
python3 run.py                                   # tables, expression, pca, de (R), capped_views, scatter
python3 run.py --only tables expression pca      # Fig. 6A without R, about 1 min
python3 run.py --only de capped_views            # DESeq2 (6 contrasts x 2, about 4 min) then the Venn
python3 run.py --only scatter                    # Fig. S5A-C
```

Inputs: the supplementary data tables (`../locations.py`). The thinning is a random draw: from the deposited tables the
smallest library rounds to 325,110 reads where the preprint's run from `quant.sf` had 325,111, so the thinned-matrix numbers
(thinned PCA, DE, outlier counts) regenerate to within sampling noise of `results/` and the native-depth numbers exactly. The
mRNA-length columns of the tables need `results/mrna_length.tsv` (gene id, length in nt) or the `quant.sf` of run SRR33776933 with
the Ensembl 113 tx2gene table; without them they stay empty. No figure depends on them.

Caveats: three BOBseq pairs against 28 to 45 for the other methods; the BOBseq contrasts are 3 vs 8 or 3 vs 10; TPM
over-corrects long mRNAs against the 3' methods (the fold changes of every BOBseq-vs-3'-method contrast fall with mRNA
length); TruSeq's within-method spread is a cross-laboratory number.

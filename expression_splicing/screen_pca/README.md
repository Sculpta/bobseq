# PCA of the 24-plex wells (Fig. 7A)

`pca_tpm.py` runs the PCA of gene expression over the human wells of the 24-plex, twice: all 21 wells (`all21`) and the 18
wells without the three input-failed `Ris dose 500mM` wells (`no_ris500`, Fig. 7A = `figures/pca_no_ris500.svg`).

Input: gene TPM of set `bobseq_pe_native` (S8: Salmon on the UMI-deduplicated reads, tximport to genes, the 19,934-gene mRNA
set, TPM re-scaled over it) and the sample sheet. Per run: log2(TPM + 1); genes with mean TPM >= 1 (12,050 in `all21`, 11,924 in
`no_ris500`); the 5,000 most variable by variance of log2(TPM + 1), re-selected per run; centred and scaled to unit variance;
PCA with full SVD, 10 components. Figures: PC1/PC2 and PC3/PC4 scores coloured by condition; a clustered heatmap of the
well x well Spearman correlation of log2(TPM + 1) over the same genes (average linkage, correlation distance) and its
dendrogram alone.

Result (`results/`): with all 21 wells PC1 (33 %) is the three Ris 500 nM wells; without them PC1 (22 %) orders the CHX wells
by dose (control, Ris 25 nM and both ASO arms overlap; CHX 1 and CHX 50 separate on PC2, 8 %). `pca_scores_<run>.tsv`,
`pca_variance_<run>.tsv`, `correlation_<run>.tsv`, `genes_<run>.tsv` (the selected genes, ranked).

```
python3 pca_tpm.py        # about 30 s; needs pandas, scikit-learn, seaborn
```

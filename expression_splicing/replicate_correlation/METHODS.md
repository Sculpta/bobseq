# Sample to sample variance calculation

## Datasets

Replicate transcriptomes were compared for three bulk RNA-seq library-preparation methods. BOBseq libraries were
generated in this study (HEK cells; a 24-plex experiment comprising untreated controls and drug-treatment conditions,
each in triplicate; full-length, randomly-primed, paired-end libraries with a 6-nt unique molecular identifier [UMI]).
prime-seq data were obtained from ArrayExpress accession E-MTAB-10142 (Janjic et al., 2022)¹; we used the eight
untreated HEK293T samples of a single lysis condition ("Incubation + Proteinase K", 10,000 cells per lysate), which are
3′-end poly-dT libraries (50-nt single-end, 10-nt UMI). DRUG-seq data were obtained from Gene Expression Omnibus
accession GSE176150 (Li et al., 2022)²; we used the 24 DMSO vehicle-control wells (24 h) of one U-2 OS plate
(VH02001704_S4), a 3′-end tag method with UMIs (Ye et al., 2018)³ originally describes the assay. One BOBseq
risdiplam-500 triplicate that failed input QC was excluded; the BOBseq Hek_control_2 well (a T-depleted UMI oligo) was
retained as it was otherwise a normal control.

## Uniform re-processing and quantification

To remove processing as a confounder, all datasets were re-processed from raw reads with an identical pipeline. Reads
were aligned to the human genome (GRCh38, Ensembl release 113) with STAR v2.7.11b⁴. For the UMI-containing libraries
(BOBseq, prime-seq, DRUG-seq), aligned reads were deduplicated to one read (or read pair) per unique combination of
mapping position, strand, and UMI. Uniquely mapped reads (MAPQ 255) were assigned to genes with featureCounts⁵ against
the Ensembl 113 gene annotation, yielding an integer gene-by-sample count matrix.

## mRNA gene set

To compare methods on the messenger transcriptome and to exclude gene classes captured very differently across
chemistries, an "mRNA" gene set was defined as protein-coding genes (Ensembl biotype protein_coding) after removing (i)
cytoplasmic and mitochondrial ribosomal-protein genes (symbols matching RPL, RPS, MRPL, or MRPS, excluding RPS6K kinases)
and (ii) all mitochondrially-encoded genes (chromosome MT), giving 19,934 genes. This exclusion was motivated by strong,
method-specific composition differences in the raw libraries (e.g., ~19% of prime-seq and ~2% of DRUG-seq gene-assigned
reads were mitochondrial, versus <1% for BOBseq; ~29% of DRUG-seq reads were ribosomal-protein), which would otherwise
dominate correlations. Counts within the mRNA set were converted to counts per million (CPM) re-normalized within that
set, and log-transformed as log₂(CPM + 1).

## Depth normalization

Analyses were performed either at native sequencing depth (no downsampling) or after depth-matching. For depth-matching,
each sample's mRNA gene-count vector was thinned to a common target by multinomial subsampling without replacement
(variance-preserving, fixed random seed). Depth targets were chosen as the minimum mRNA-assigned depth of a defined
reference set — the smallest prime-seq library (296,992 counts), the smallest of the eight highest-depth DRUG-seq wells
(141,053 counts), or the smallest library across the full sample set (48,837 counts) — as specified per analysis.

## Gene inclusion and reproducibility metric

For each pair of replicate samples, genes were included either by a detection threshold (≥1, ≥3, or ≥5 counts in both
samples of the pair) or by a rank threshold (the N most highly expressed genes for that pair, ranked by the mean of the
two samples' CPM, with no count cutoff). When a common rank threshold across a shallow set was used, N was capped at the
largest value for which every pair had that many genes detected in both samples, to avoid inflating correlations with
dropout (genes present in only one replicate). Replicate-to-replicate reproducibility was quantified as the Pearson
correlation coefficient of log₂(CPM + 1) over the included genes. For BOBseq, correlations were computed only between
replicates within the same experimental condition (the triplicate of each of six conditions) and pooled across
conditions; for prime-seq and DRUG-seq, which are each a single condition, all within-method replicate pairs were used.

## Statistics

Because replicate pairs are not independent (each sample contributes to multiple pairwise comparisons), significance
was assessed on a per-sample summary rather than on individual pairs: for each sample we computed its mean correlation to
the other replicates of its group, and compared these per-sample distributions between methods with a two-sided
Mann–Whitney U test. Interquartile range of the per-pair correlations was reported as a measure of within-method
consistency. Analyses used Python 3 with NumPy, SciPy⁶ (scipy.stats.mannwhitneyu), and Matplotlib.

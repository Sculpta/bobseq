# Expression and splicing analyses

The analyses behind the expression and splicing panels of the preprint: the comparison of BOBseq with published HEK293T
datasets (Fig. 6A, 6B, 6D, S5A-C, S6A) and the analysis of the 24-plex drug screen (Fig. 7A, 7C, 7D, 7E, 7F). Every analysis
starts from the supplementary data tables deposited with the preprint (S1 to S11) and, where reads are needed, from the
per-sample BAMs of the benchmark; nothing here needs the raw sequencing data or a cloud account.

| Directory | Preprint panels | What it computes | Inputs |
|---|---|---|---|
| `method_comparison/` | Fig. 6A, 6B, S5A-C | BOBseq's three control wells against the prime-seq, BRB-seq and TruSeq controls: PCA, DESeq2 against TruSeq on a capped common gene set (signed Venn), gene-vs-gene scatters | S7, S8, S9; R with DESeq2 for the Venn |
| `replicate_correlation/` | Fig. 6D, S6A | replicate-to-replicate Pearson r and RMSE, BOBseq vs prime-seq vs DRUG-seq | S11, the Ensembl 113 GTF |
| `screen_pca/` | Fig. 7A | PCA of the 24-plex wells on gene TPM | S8 |
| `screen_splicing/` | Fig. 7C | differential splicing of the 24-plex, five contrasts, on junction fragment counts | S6 (S4 alongside) |
| `smn2_exon7/` | Fig. 7D | paralog-aware SMN2 exon 7 inclusion from the reads at both SMN loci | per-sample BAMs, genome FASTA |
| `splicing_tracks/` | Fig. 7E, 7F | depth-matched coverage tracks with junction arcs for the CHX hits and the positive-control exons | per-sample BAMs, GTF, `screen_splicing/results/` |
| `supplement_tables/` | | how the supplementary data tables were assembled (provenance only) | |

Each directory has its own README with the method, the numbers behind its panels and the commands. `sample_metadata.csv`
(= S1) is the sample sheet every analysis reads: one row per sample-arm, `acc_number` = `<set>__<sample>` is the column name
in every count matrix. Every analysis writes its figures to `<directory>/figures/` (not part of the repository) and its
tables to `<directory>/results/`, where the tables of the preprint's run are kept.

## Inputs

**Supplementary data tables.** Unpack the archive into `expression_splicing/data/supplement_tables/` (or point
`BOBSEQ_TABLES` at the directory; see `locations.py` for every location and its environment variable):

| file | rows x columns | content |
|---|---|---|
| `S1_sample_metadata.csv` | 119 sample-arms | the sample sheet (a copy is `sample_metadata.csv` here) |
| `S2_public_truseq_runs.tsv` | 10 runs | the public TruSeq HEK293T runs of the method comparison |
| `S3_psi_dedup.csv.gz`, `S4_junction_counts_dedup.csv.gz` | 8,701 junctions x 119 | PSI and counts on UMI-deduplicated molecules (empty = not detected) |
| `S5_psi_fragments.csv.gz`, `S6_junction_counts_fragments.csv.gz` | 36,433 junctions x 119 | the same per fragment, without deduplication: the unit of `screen_splicing` |
| `S7_gene_counts_salmon.csv.gz`, `S8_gene_tpm_salmon.csv.gz`, `S9_gene_effective_length_salmon.csv.gz` | 19,934 genes x 129 | Salmon / tximport gene counts, TPM and effective length on the mRNA gene set; the 119 arms plus the 10 TruSeq runs |
| `S10_mrna_gene_set.tsv` | 19,934 genes | the mRNA gene set: Ensembl 113 protein-coding genes without the 13 mitochondrially encoded genes and the 169 ribosomal-protein genes |
| `S11_gene_counts_featurecounts.csv.gz` | 78,932 genes x 119 | featureCounts per Ensembl gene, all biotypes |

Row ids of S3 to S6 are `chr<N>:GENE:s|t:anchor|<N>:start-end:strand`: the LSV (local splicing variation: the junctions that
share a source, `s`, or target, `t`, splice site) and the junction (contig without prefix, intron start 0-based, intron end,
strand). PSI = fragments (or molecules) of the junction over the sum of the detected junctions of its LSV side.

**Sets and samples.** `bobseq_pe_native` (the 21 human wells of the 24-plex, both mates, 2 x 150 nt), `bobseq_50nt` (read 2
cut to 50 nt, the 18 benchmark wells), `drugseq` / `drugseq_native` (24 DMSO wells, 50 / 52 nt), `primeseq` /
`primeseq_native` (8 samples, 50 nt; identical), `brbseq` / `brbseq_native` (8 control wells of GSE334309, 50 / 80 nt). The
sample-sheet labels of the risdiplam wells, `Ris dose 25mM` and `Ris dose 500mM`, carry the legacy unit of the run
configuration; the doses are 25 nM and 500 nM. The three `Ris dose 500mM` wells are input failures and enter no figure;
`Hek_control_2` (`in_benchmark = no`, a T-depleted UMI oligo) is a normal HEK293T well by every other measure and is the third
control of every analysis here.

**Reads.** `smn2_exon7` and `splicing_tracks` read the per-sample BAMs of the benchmark, `<set>/<sample>.bam`
(`benchmark/build_per_sample_bams.sh`: primary alignments on the human contigs, Ensembl contig names, not deduplicated), and make
UMI-deduplicated copies with `benchmark/scripts/dedup_umi_position.py` where they need them.

## How the tables were made

Alignment as in the benchmark (STAR 2.7.11b, combined GRCh38.113 + GRCm39.113 index, human contigs, primary alignments),
one per-sample BAM per set. UMI + position deduplication (`benchmark/scripts/dedup_umi_position.py`: one molecule per UMI,
contig, strand and fragment position), MAPQ 255. Junctions with `benchmark/scripts/junctions_by_fragment.py` (the intron of
every N in a CIGAR, anchors of at least 4 nt, counted once per read name), mapped onto the LSVs of the junction catalogue
the atlas is built on; PSI per LSV side over the detected members, with the gates LSV total >= 3 per sample, junction maximum
>= 3 in >= 3 samples of a set, >= 2 members (S3 / S4 on the deduplicated BAMs, S5 / S6 on the raw BAMs). Gene counts:
featureCounts 2.1.1 (`-s 1`, MAPQ 255, primary, pairs once) on the deduplicated BAMs (S11); Salmon 1.10.3 (`quant -l A
--validateMappings`, Ensembl 113 cDNA index) on the reads of the deduplicated BAMs, summarised to genes with tximport 1.40
(counts, TPM, effective length), alternative-haplotype copies folded onto their gene, restricted to the mRNA gene set (S7 to
S9); the 10 TruSeq runs were trimmed with Trim Galore and quantified with Salmon 1.11.4 the same way. BRB-seq (GSE334309) was
processed like the other public sets; it is not an arm of the repository's benchmark.

## Software

Python 3 with numpy, scipy and matplotlib for everything; pandas, scikit-learn and seaborn for `method_comparison` and
`screen_pca`; pysam and samtools for `smn2_exon7` and `splicing_tracks`; BLAST+ (blastn) for `smn2_exon7`; R with DESeq2 and
apeglm for the DESeq2 stage of `method_comparison`. Text in every SVG stays text (`svg.fonttype none`); Arial is requested
and, where it is absent, the metric-identical Liberation Sans is registered under that name.

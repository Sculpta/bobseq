# SMN2 exon 7: paralog-aware quantification of the risdiplam readout (Fig. 7D)

**Question.** Does the risdiplam arm show the expected rise in SMN2 exon 7 inclusion, and can a random-primed full-length
method read it when the reads that skip exon 7 are perfect SMN1/SMN2 multimappers? Fig. 7D = `figures/SMN_junction_reads_pct.svg`.

**The problem.** SMN1 and SMN2 are 99.86 % identical over 29 kb; inside exon 7 they differ at one base, c.840 C>T. A read that
includes exon 7 carries c.840 and can be assigned to a paralog; a read that skips it (ex6>ex8) contains no exon-7 sequence, so
STAR places it at one locus at random with MAPQ 3, and a MAPQ-255 filter drops it. This analysis goes back to the BAMs at both
loci with all MAPQ values. Exon numbering follows the SMN literature (exon 7 = the 54-nt alternatively spliced exon = Ensembl
exon 8 of the canonical transcripts).

## Stages (`run.py --list`)

1. `blast_loci`: the two loci (SMN2 chr5:70,049,000-70,078,600, SMN1 chr5:70,924,500-70,953,600) cut from the reference and
   aligned against each other with blastn; the collinear alignment gives 42 differences (24 mismatches, 18 indels), the
   paralog-diagnostic sites (`results/SMN_diagnostic_sites.tsv`).
2. `region_bams`: both loci, all MAPQ, sliced out of every per-sample BAM, raw and UMI-deduplicated (made on demand with
   `benchmark/scripts/dedup_umi_position.py`).
3. `fragments`: every fragment (both mates) typed by its junctions (ex6>ex7, ex7>ex8 = inclusion; ex6>ex8 = skip) and by the
   diagnostic bases it covers (SMN2 / SMN1 / ambiguous / conflict), plus the base at c.840 and at c.*239 in exon 8.
4. `blast_reads`: cross-check of the site-based paralog calls by blasting every read against both loci.
5. `psi`: per well and per condition (all wells pooled), Clopper-Pearson 95 % intervals. A skip fragment belongs to no paralog:
   the skip pool is SMN exon-7 skipping as a whole, and each paralog's PSI is its own inclusion against that shared pool:
   psi_SMN2 = incl_SMN2 / (incl_SMN2 + skip_all). Fisher's exact test on the pooled junction reads and a Welch t-test on the
   per-well PSIs against the controls (`results/SMN2_tests.tsv`).
6. `simple`: a read-only inclusion index from the bases at c.840 and c.*239, no junctions.
7. `figures`: from `results/SMN2_psi_conditions.tsv` and `SMN2_psi_samples.tsv` (UMI-deduplicated, BOBseq 2 x 150 nt, Ris 500 nM
   excluded): `SMN_junction_reads_pct.svg` (Fig. 7D: every event as a percentage of the condition's exon-7 junction reads, one
   100 % bar per condition), `SMN_junction_reads.svg` (the counts), `SMN_inclusion_ratio.svg` (SMN2 : SMN1 ex6>7 reads, Katz
   interval, dots = wells).

## Result (BOBseq 2 x 150 nt, deduplicated molecules, wells of a condition pooled; `results/SMN2_psi_conditions.tsv`)

| condition | SMN2 PSI [95 % CI] | SMN1 PSI | SMN2 ex6>7 / ex7>8 · SMN1 ex6>7 / ex7>8 · ambiguous ex6>7 / ex7>8 · shared skip |
|---|---|---|---|
| control | 0.44 [0.27-0.62] | 0.73 | 9 / 13 · 34 / 46 · 2 / 10 · 19 |
| Ris 25 nM | 0.78 [0.58-0.91] | 0.85 | 20 / 14 · 20 / 26 · 0 / 4 · 6 |
| CHX 1 / 50 µg/mL | 0.34 / 0.36 | 0.63 / 0.63 | control-like |
| ASO 1 / 100 nM | 0.30 / 0.41 | 0.64 / 0.56 | control-like |

Risdiplam raises SMN2 exon 7 inclusion from 0.44 to 0.78 (Fisher exact p = 0.0098 on the pooled molecules; the per-well Welch
test, 3 vs 3 wells of 5 to 25 molecules each, gives p = 0.28). Under Ris 25 nM the SMN2 ex6>7 share of the exon-7 junction reads
rises from 7 % to 22 % while skipping falls from 14 % to 7 %. The other 98 sample-arms are in the tables: BRB-seq and prime-seq
reach exons 7 and 8 with few informative molecules, DRUG-seq has no reads at the locus.

## Run

```
python3 run.py --only blast_loci region_bams fragments blast_reads    # needs the per-sample BAMs, the genome FASTA, samtools, blastn, pysam
python3 run.py --only psi simple figures                              # seconds, from results/
```
Caveats: n = 3 wells per condition with few molecules each, so every inference is on pooled reads; the skips are not attributed
to a paralog; the inclusion ratio and the c.840 T-fraction move with the SMN2/SMN1 expression ratio as well and are not PSI.

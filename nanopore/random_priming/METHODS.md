# Methods: random-priming bobcode (3′TSO) analyses

## Libraries and samples

Random-priming bobcode libraries use the same 3′TSO **bob-v1** construct as the poly(dT) 3′TSO
libraries (identical template-switch oligo, 12-nt anchor `CTCTTCCGATCT`, and 7-mer species bobcodes),
differing only in the reverse-transcription primer: a random N-mer (R1-6N or R1-9N) replaces the
oligo-d(T) primer. Nine human/mouse species-mixing samples across three sequencing runs were
analysed:

| run | date | RT primer | barcodes (condition) |
|---|---|---|---|
| 260819_BC9-11_6N | 2026-08-19 | R1-6N | BC09 (TM1 50 °C 30 min), BC10 (TM2 temp-cycling), BC11 (TM3 50 °C 60 min) |
| 260808_9N_exoi_full | 2026-08-08 | R1-9N | BC01 (ExoI), BC02 (ExoI 20 min), BC03 (ExoI 10 min), BC04 (ExoI + RNeasy) |
| 260827_BC14_16 | 2026-08-27 | R1-6N (after poly(A) selection) | BC15 (1× poly(A)), BC16 (2× poly(A)) |

The RACK1 genome-browser panel compares two representative samples: **260730 bc12** (poly(dT)
3′TSO; run 260730_RT_conditions_BC11-13, barcode12) and **260827 BC15** (random-priming 3′TSO; run
260827_BC14_16, barcode15).

## Alignment

Reads (Oxford Nanopore) were trimmed of the ONT/native-barcode flank and aligned with **STARlong
2.7.11b** to a combined human+mouse genome index built from *Homo sapiens* GRCh38 and *Mus musculus*
GRCm39 (Ensembl release 113), with contigs prefixed `HUMAN_`/`MOUSE_` so that read species is given
by the chromosome of alignment. Local alignment retains the adapter and bobcode as soft-clipped
sequence. A single dense (non-sparse) index was used; one STARlong process was run at a time.

## Barcode accuracy

Barcode accuracy was computed exactly as for the other bobcode chemistries: for each ONT barcode,
species-unique protein-coding (mRNA) reads passing low-complexity, single-barcode-unit and
exact-code filters were scored, and a read was counted correct when the species declared by its
7-mer bobcode matched the species (`HUMAN_`/`MOUSE_`) of its alignment. Ribosomal-RNA and
mitochondrial reads were excluded (random priming yields 60–70 % rRNA, which is expected and does
not enter scoring). Accuracy is reported on all scored mRNA reads with a 95 % Wilson interval; reads
were capped at 10,000 per sample. The random-priming construct required no code change — only a
run-rules entry (`code/barcode_accuracy/bobcode_run_rules_randompriming.tsv`).

## Gene-body (5′→3′) coverage

Per-transcript-normalized gene-body coverage was computed with **Picard CollectRnaSeqMetrics 3.4.0**
(`STRAND_SPECIFICITY=NONE`) run twice per library (a `HUMAN_`- and a `MOUSE_`-prefixed refFlat) so
that human- and mouse-gene metagenes are obtained from the same species-mixed alignment. Gene models
were protein-coding transcripts of the merged Ensembl-113 annotation with ribosomal-RNA and
mitochondrial genes excluded; cytoplasmic and mitochondrial ribosomal-protein genes were
additionally removed by anchored gene-symbol match (`RPL*`, `RPS*`, `RPLP*`, `RPSA`, `MRPL*`,
`MRPS*`, `FAU`; the `RPS6K*` kinases were retained). Prior to quantification, alignments were
restricted to primary records with an aligned cDNA length ≤ 500 bp (aligned length = sum of CIGAR
`M/I/=/X`). Picard `MINIMUM_LENGTH` was set to **1000**, so only transcripts ≥ 1000 bp contribute
to the metagene; this removes short-transcript artefacts under short-insert selection. Profiles are
the mean of per-sample normalized-coverage histograms (101 positions), with a 95 % confidence
interval (t-distribution). Samples entered the summary only if **both** the human and the mouse track
had ≥ 1000 primary reads passing the ≤ 500-bp filter; if either failed, the whole sample was excluded.

## 5′→3′ coverage balance

For each sample and species, a scalar balance was defined as **twice the coverage centroid** of the
gene-body metagene, `2 · Σ(pᵢ·cᵢ) / Σcᵢ`, where `pᵢ ∈ [0,1]` is the 5′→3′ position and `cᵢ` the
normalized coverage. This equals 1 for an evenly balanced profile, is < 1 for 5′-leaning and > 1 for
3′-leaning coverage, and is bounded on [0,2] (robust to the near-zero 5′ coverage of strongly
3′-piled poly(dT) libraries, for which a simple 3′/5′ ratio diverges).

## Single-gene coverage (genome browser)

Per-base coverage over RACK1 (and ACTB, GAPDH) was computed from primary alignments with aligned
cDNA length ≤ 500 bp and displayed against the full set of annotated transcript isoforms
(Ensembl 113), coloured by transcript biotype.

## Software

STARlong 2.7.11b; samtools 1.20; Picard 3.4.0; Python 3.13 with NumPy, SciPy and Matplotlib.
References: Ensembl release 113 (GRCh38, GRCm39).

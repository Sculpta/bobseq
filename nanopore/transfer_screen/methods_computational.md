# Methods: computational analysis

### Nanopore sequencing and basecalling

Barcoded cDNA libraries from mixed human and mouse input were barcoded with the Native Barcoding
Kit 24 V14 (SQK-NBD114.24) and sequenced on a MinION Mk1D with R10.4.1 flow cells. Reads were basecalled
and demultiplexed in MinKNOW with the high-accuracy model
`dna_r10.4.1_e8.2_400bps_hac@v5.2.0`. All passing reads of each barcode were used, with one
exception: for one 22-hour run (three barcodes, 0.94–1.43 million reads each), only the
first 10,000 reads per barcode were analysed. Where the same libraries had been sequenced in
two MinKNOW runs, the reads were pooled (two sets of libraries, 15 datasets). Read identifiers were checked to be
unique within every dataset.

### Species barcodes

Each library carries a species code: one sequence for human input and one for mouse input.
It sits at a fixed position immediately 3′ of a constant primer sequence (the "anchor"),
in one of four chemistries:

- **poly(dT) priming:** a 6-nt code (or 12-nt, in one 90-mer design) in the polyT RT primer;
- **splint-mediated ligation:** a 4-nt code (or 8-nt BOB-A code) on a primer ligated to the
  RT primer through a splint oligonucleotide;
- **TSO bobcode (5′DBCO):** a 7-nt code (or 4-nt) in a Nextera-R1 or SMARTer TSO;
- **TSO bobcode (3′DBCO):** a 7-nt code after a SMARTer TSO backbone or a TruSeq R1 handle.

The anchor, the two codes and their species assignment for every run are listed in the run
table of the code repository (`bobcode_run_rules.tsv`).

### Read processing and alignment

ONT adapter and native-barcode sequences were removed from both read ends. The 5′ trim ran up
to and including the post-barcode flank `CAGCACCT`, searched within the first 80 nt, and the
3′ trim started at its reverse complement, searched within the last 80 nt. Trimmed reads were
aligned with STARlong 2.7.11b to a combined human and mouse genome:
- **Genome:** GRCh38 and GRCm39 primary assemblies from Ensembl, chromosome names prefixed
  `HUMAN_` and `MOUSE_`.
- **Annotation:** Ensembl release 113 for both species.
- **Index:** built with `--sjdbOverhang 99` and a dense suffix array.

Alignment parameters were `--outFilterMultimapNmax 10000 --outFilterScoreMinOverLread 0.3
--outFilterMatchNminOverLread 0.3 --outFilterMismatchNmax 1000
--outFilterMismatchNoverLmax 0.3 --alignIntronMax 1000000 --alignSJDBoverhangMin 1
--alignEndsType Local --seedPerReadNmax 100000`, with all alignments of a read reported.
Each barcode was aligned separately with identical settings.

### Read classification

A read was assigned to **human** or **mouse** only if all of its alignments (primary,
secondary and supplementary) were on that species' chromosomes. Reads with alignments on
both species were classed as ambiguous and not scored.

Each mapped read was then placed in one RNA category. The first matching rule applied:
1. **rRNA:** any alignment overlapping one of 99 empirically defined rDNA loci (provided with
   the code), or a primary alignment on an exon of a gene with an rRNA biotype.
2. **Mitochondrial:** a primary alignment on the mitochondrial chromosome.
3. **Other:** a primary alignment on an exon of a ribosomal-protein gene (RPL*, RPS*, MRPL*,
   MRPS*; 342 genes).
4. **mRNA:** a primary alignment on an exon of any other protein-coding gene.
5. **Other:** everything else, including lncRNA and other non-coding, intronic and
   intergenic reads.

Exon overlap was assessed on aligned blocks, ignoring strand.

Reads whose aligned bases had a normalized dinucleotide Shannon entropy below 0.65 were
classed as low-complexity and excluded from scoring. These are alignments that rest on
homopolymers or simple repeats.

### Barcode reading and accuracy

The species code was read from the bases immediately 3′ of an exact 12-nt match to the
run's anchor, searched on both strands. A barcode unit was an anchor followed by a sequence
within one mismatch of either of the sample's two codes. Reads with no unit, or with two or
more units (chimeras), were not scored. In reads with a single unit, the code had to match
one of the two codes exactly, with no mismatches or indels. Codes of other libraries were
ignored.

A read was **scored** if it passed all of these filters:
- species-unique (human or mouse);
- mRNA;
- not low-complexity;
- a single exact code.

It was **correct** if its code's species matched its aligned species. Barcode accuracy is
correct / scored reads, with a 95% Wilson score interval. It was also reported separately
for reads aligned to human and to mouse.

For each dataset, the code × species table was inspected: a code whose declared species was
not the majority of its reads was flagged, but accuracy was always reported. Datasets with
fewer than 100 scored reads were excluded (11 datasets). So were designed negative controls,
single-species controls and single-code designs, where accuracy cannot be defined. The
excluded datasets and reasons are listed in Supplementary Table S1.

In total, 512,565 reads were scored across 186 datasets (median 1,639 per dataset). Median
barcode accuracy was:
- poly(dT) priming: 69.9% (59 datasets);
- splint-mediated ligation: 49.9% (34);
- TSO bobcode (5′DBCO): 80.0% (20);
- TSO bobcode (3′DBCO): 99.1% (73).

In the 3′DBCO chemistry, 45 of 73 datasets reached at least 99%.

### Known limitation

In two constructs, one species code is identical to the barcode position of the standard
(unbarcoded) SMARTer TSO. These are the human code AGTACAT of the SMARTer 7-nt TSO design,
and the mouse code TACA of one 4-nt design. Reads from an unbarcoded TSO would therefore be
counted as carrying that code. The accuracies of these constructs were in line with those of
constructs without the collision.

### Software

Analysis: Python 3.12 (standard library only), samtools 1.20, STARlong 2.7.11b.
Figure and table: matplotlib 3.10 and openpyxl 3.1.

The complete pipeline is available in the code repository accompanying the preprint:
- the processing program;
- the general rules and the per-run table;
- the rDNA loci;
- the reference-building instructions with contig checksums;
- the figure and supplementary-table scripts.

It reproduces the figure and Supplementary Table S1 from the deposited reads.

### Data availability

Raw nanopore reads (186 datasets, one FASTQ file per dataset) are deposited in the NCBI
Sequence Read Archive under BioProject PRJNA1532586. `fastq_manifest.tsv` lists every
dataset with its read count, checksum, source files and run accession, and
`tools/build_fastq_repository.py` documents how MinKNOW output files were joined and which
datasets were pooled or subsampled.

---

### Figure legend

**Barcode accuracy of species bobcodes across four chemistries.** Each point is one library
(one ONT barcode) from a human + mouse mixing experiment. Panels, left to right: poly(dT)
priming, splint-mediated ligation, TSO bobcode (5′DBCO) and TSO bobcode (3′DBCO). Within
each panel, libraries are ranked by accuracy and numbered from 1 (Supplementary Table S1
lists every library). Barcode accuracy is the fraction of scored mRNA reads whose species
code matches the species they align to. Error bars are 95% Wilson confidence intervals from
the number of reads scored (counting error only). The y axis is log-scaled in inaccuracy
(100 − accuracy) and runs from 40% to 99.5%. Points beyond these limits are drawn as
arrowheads at the axis edge. Shading marks accuracy above 99% (bobcode benchmark), 96–99%
(DRUG-seq benchmark), 90–96% and below 90%. The red line marks 50%, the accuracy expected if
codes were assigned at random.

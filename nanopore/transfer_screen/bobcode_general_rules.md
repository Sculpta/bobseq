# Bobcode barcode accuracy (splint, polydT, 5'TSO, 3'TSO): general rules

These rules apply to every splint, polydT, 5'TSO and 3'TSO dataset. Anything that differs
between datasets (chemistry, construct/anchor, the two codes and their species,
input files, barcodes to score) lives in `bobcode_run_rules.tsv`, and nowhere else.
The code (`speciesmix_bobcode_metrics.py`) must follow this file; if the code and
this file disagree, the file wins and the code is fixed. Every cell of the run file is
filled: `none` where there is nothing to record, `all` in `max_reads` when every read is
used, and `not scored (see barcodes_excluded)` for runs listed only to record why they
are not scored.

In polydT chemistry the species code is in the RT primer (3' end of the cDNA, next
to the polyT); in splint chemistry it is in a BOB primer ligated onto the RT primer
through a splint (also the 3' end); in 5'TSO and 3'TSO it is in the template-switch
oligo (5' end). The rules below treat them identically: find the anchor, read the
code after it.

Key choices: exact match everywhere; only the sample's two codes count; protein-coding
exons only (ribosomal-protein genes excluded); low-complexity reads removed; accuracy is
always reported.

---

## 1. Input reads (one ONT barcode at a time)

1. Read the FASTQ file(s) named by the run file's `raw_fastq_input` pattern(s), e.g.
   `barcodeNN/*.fastq.gz` (`barcodeNN` becomes `barcode06` etc.). Several patterns are
   separated by ` ; `; a pattern is relative to the run folder unless it is an absolute
   path. Concatenate all files (gzipped or plain). A read ID seen more than once is
   kept once, so pooled sources can never double-count a read. If the run file sets
   `max_reads`, only the first that many reads are used (files in name order, reads in
   file order) and the rest are ignored (260120: first 10,000 per barcode).
   With `--fastq-dir`, the dataset's single file from the repository FASTQ folder
   (`<chemistry>_<construct>/<YYMMDD>_BC<NN>.fastq.gz`: the same original reads, joined
   into one file) is read instead.
2. Trim the ONT adapter + native barcode from both ends of each read:
   - 5' end: find `CAGCACCT` in the first 80 bp; remove everything up to and including it.
   - 3' end: find `AGGTGCTG` in the last 80 bp (last occurrence); remove it and everything after.
   - If the sequence is not found, that end is left as is.
   - Reads that become empty are dropped.

## 2. Alignment

- Aligner: STARlong 2.7.11b, official binary (`tools/STARlong`), under Rosetta.
- How to download and build the reference files and the index: `README.md` §2.
- Index: dense combined human + mouse
  `references/indexes/star/combined_GRCh38.113_GRCm39.113/STAR_index`
  (chromosomes prefixed `HUMAN_` / `MOUSE_`). Never the sparse index.
- Each barcode's trimmed FASTQ is aligned on its own (no concatenating barcodes).
  One STARlong at a time, serially.
- Parameters (same for every dataset):
  ```
  --outSAMtype BAM Unsorted --outSAMattributes NH HI AS nM --outSAMunmapped Within
  --outSAMmultNmax 18446744073709551615 --outFilterMultimapNmax 10000
  --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3
  --outFilterMismatchNmax 1000 --outFilterMismatchNoverLmax 0.3
  --alignIntronMax 1000000 --alignSJDBoverhangMin 1 --alignEndsType Local
  --seedPerReadNmax 100000
  ```

## 3. Species of each read (from the alignment)

- Consider every alignment record of the read (primary, secondary, supplementary).
- **human**: all records on `HUMAN_` chromosomes. **mouse**: all on `MOUSE_`.
- **ambiguous** (records on both) and **unmapped** reads are not scored.

## 4. RNA category (only protein-coding mRNA is scored)

Decided in this order; the first match wins:

1. Any alignment record overlaps an rDNA locus (`references/annotations/rdna_loci.bed`) → rRNA.
2. Primary alignment on a `*_MT` chromosome → MT (mitochondrial).
3. Primary alignment overlaps an exon of a gene with biotype `rRNA`, `Mt_rRNA` or `rRNA_pseudogene` → rRNA.
4. Primary alignment overlaps an exon of a ribosomal-protein gene → other. These are
   `protein_coding` genes named RPL*/RPS*/MRPL*/MRPS* (case-insensitive), excluding
   the RPS6K* kinases. There are 342 such genes (human + mouse).
5. Primary alignment overlaps (≥ 1 bp) an exon of any other `protein_coding` gene → **mRNA**.
6. Anything else (lncRNA and other non-coding, intron-only, intergenic) → other.

"Overlaps" uses the alignment's aligned blocks (introns skipped by `N` are not
counted) and ignores strand. Only category **mRNA** is scored. Annotation:
`references/annotations/combined_genome.gtf`, which is identical to Ensembl
GRCh38.113 + GRCm39.113 (all 78,932 human and 78,298 mouse genes checked).

## 5. Low-complexity reads (removed)

Take the bases of the primary alignment that aligned to the genome (soft-clips
excluded). If their normalized dinucleotide Shannon entropy is **< 0.65**
(entropy / log2(16); sequences < 3 nt count as complex), the read is not scored.
This removes reads whose species call rests on a homopolymer or simple repeat.

## 6. Reading the barcode

1. The construct (run file) sets the **anchor**: the 12 nt immediately 5' of the
   code in the oligo. The code itself follows the anchor directly.

   | Construct | Chemistry | Oligo layout (5'→3') | Anchor | Code length |
   |---|---|---|---|---|
   | `bob-c-splint` | splint | BOB-C5 primer: bob-c `ATGTAGTCCGTACGCATGTC` + code + BOB_HR `ACATGGTAGC`, splint-ligated to the RT primer | `CGTACGCATGTC` | 4 |
   | `bob-a` | splint | BOB-A primer: TruSeq R1 `…GCTCTTCCGATCT` + code + `ACCTGGTAGC` | `CTCTTCCGATCT` | 8 |
   | `bob-polydt` | polydT | bob-c primer `ATGTAGTCCGTACGCATGTC` + code + polyT(-VN) | `CGTACGCATGTC` | 6 |
   | `bob-k-polydt` | polydT | TruSeq R1 `…GCTCTTCCGATCT` + code (3 × 4 bp) + 16N UMI + polyT | `CTCTTCCGATCT` | 12 |
   | `tso-bob` | 3'TSO | SMARTer backbone `AAGCAGTGGTATCAACGCAG` + code + rGrGrG | `GTATCAACGCAG` | 7 |
   | `bob-v1` | 3'TSO | TruSeq R1 handle `…CTCTTCCGATCT` + code + UMI/spacer + rGrGrG | `CTCTTCCGATCT` | 7 |
   | `nextera-tso` | 5'TSO | Nextera-R1 backbone `GTGACTGGAGTTCAGACGTG` + code + rGrGrG | `AGTTCAGACGTG` | 7 |
   | `nextera-tso (BOB-J)` | 5'TSO | 19-nt backbone `GTGACTGGAGTTCAGACGT` + code + rGrGrG | `GAGTTCAGACGT` | 8 / 10 |
   | `smarter-4bp` | 5'TSO | SMARTer primer `GCAGTGGTATCAACGCAGAG` + code (+ `TGG` in reads) | `ATCAACGCAGAG` | 4 |

   The code length is simply the length of the codes in the run file.
2. Search both strands for the anchor, **exact match** (no mismatches, no indels).
   The bases immediately 3' of each anchor hit form a candidate for each code.
3. A **barcode unit** = anchor + a candidate within 1 mismatch of either of the
   sample's two codes. Hits < 27 nt apart on the same strand count as one unit.
   (Tolerance 1 is used only to find units, so a near-miss second copy still
   counts as a chimera. It never makes a read scorable. With 4-bp codes that differ
   at 2 positions, a candidate can be within 1 of both codes; it is still one unit.)
4. Only the sample's **two codes** count. Every other bobcode (including the
   rest of any 24-code panel) is ignored.
5. The read is not scored if it has **0 units**, or **≥ 2 units** (TSO–TSO
   chimera, including a read carrying both codes).
6. For single-unit reads, the candidate must equal one of the two codes **exactly**
   (0 mismatches, no indels). Otherwise the read is not scored.

## 7. Scoring

- A read is **scored** if it passes 3–6: species human or mouse, mRNA, not low
  complexity, exactly one unit, exact code.
- **Correct** = the species assigned to its code (run file) equals its aligned species.
- **Accuracy** = correct / scored (%), with a 95% Wilson interval on the scored-read count.
- Report the human-RNA and mouse-RNA subsets as well (scored and correct for each).
- **Minimum size:** a dataset with fewer than **100** scored mRNA reads is excluded. The
  program still writes its numbers but marks the `status` column
  `EXCLUDE: fewer than 100 scored mRNA reads`; the dataset is then moved to
  `barcodes_excluded` in the run file and its FASTQ is not uploaded.

**Not required** for accuracy: GGG template-switch, UMI/spacer, barcode-to-insert
distance, insert length, polyT, C28, TruSeq ends, full molecule. There is no UMI
de-duplication.

## 8. Species direction note

Accuracy is **always** reported; nothing is withheld. For every ONT barcode the
code × aligned-species table over the scored reads is written to the output, and
if a code's declared species is not the majority of that code's reads (codes with
≥ 20 scored reads), the `status` column carries a note saying so. A human/mouse
mix-up in the run file would show as accuracy far below 50% and a note on both
codes; a note on one code only means the codes do not separate species well
(e.g. 260519). The direction in the run file is only ever changed by a person.

## 9. Outputs

`speciesmix_bobcode_metrics.py` writes one row per ONT barcode to `bobcode_metrics.tsv`:

| Column | Meaning |
|---|---|
| `chemistry`, `run`, `barcode`, `plot_no`, `construct`, `human_code`, `mouse_code` | from the run file |
| `sample_note` | the barcode's entry in the run file's `sample_notes` (`BC05: text \| BC06: text`), taken from the current run file each time the table is written; `condition variant` when the barcode has no description. The run file's `barcode_flags` column (same `BC05: text` format) holds per-barcode flags such as contamination notes; `make_supplementary_table.py` prints them in the Notes column |
| `reads_total` | reads in the raw input FASTQ(s) |
| `reads_mapped` | reads with ≥ 1 alignment (= mRNA + rRNA + MT + other) |
| `mRNA_reads`, `rRNA_reads`, `MT_reads`, `other_reads` | §4 categories over all mapped reads, any species |
| `human_mRNA_reads`, `mouse_mRNA_reads` | mRNA reads aligned only to human / only to mouse (before barcode filters) |
| `human_mRNA_accuracy_pct`, `human_mRNA_scored` | §7 accuracy on human mRNA reads, and the n it is based on |
| `mouse_mRNA_accuracy_pct`, `mouse_mRNA_scored` | the same for mouse |
| `total_mRNA_accuracy_pct`, `total_mRNA_scored`, `total_ci95_lo_pct`, `total_ci95_hi_pct` | all scored mRNA reads, with the 95% Wilson interval |
| `status` | `ok`, or a note from §8 |

Each barcode also gets a `metrics.json` file with the count removed at every filter
step, the code × species table, and provenance: SHA-256 checksums of the script and
both rules files, STARlong and samtools versions, index and GTF paths, and the date.

---

## 10. Known limitations

**splint: no collision in these libraries.** The splint code ATCG is the first 4 bases of the polydT
code ATCGCT behind the same bob-c anchor, so a polydT BOB-E-1 molecule would read as
ATCG. The BOB-E polydT primers were first used on 2026-05-01, after every splint run
(the last is 260424), so they cannot be in the splint libraries.

**polydT: no collision.** No other bob-c-anchored oligo in the lab's list is within
1 mismatch of a polydT code (closest: 2 mismatches). Where a library also carries the
other polydT pair as carry-over (260605 BC15–24, 260611), those molecules match neither
of the sample's codes and are simply not scored.

**TSO chemistries:**

Two codes are the same sequence as the barcode position of the standard (unbarcoded)
SMARTer TSO and exoTSO (`…CAACGCAGAGTACATrGrGrG`), so an unbarcoded TSO product cannot
be told apart from a coded one:

- **tso-bob human code AGTACAT** (3'TSO, 7 runs). In the data, accuracy is not lower than
  in bob-v1, which has no collision (median 99.37% in 28 tso-bob datasets, 99.02% in 45
  bob-v1 datasets), and mouse-RNA error exceeds human-RNA error in both (medians 1.11% vs
  0.30% in tso-bob, 1.26% vs 0.51% in bob-v1). The one tso-bob sample with exoTSO added on
  purpose (260605 BC02) has 99.85% accuracy (0.25% error on mouse RNA).
- **smarter-4bp mouse code TACA** (5'TSO, 260410): the 4 bases after the `ATCAACGCAGAG`
  anchor of the standard SMARTer TSO are `TACA`, so unbarcoded TSO scores as mouse. One
  plotted point (260410 BC05); not testable from the data.

The Nextera-backbone constructs (all other 5'TSO runs) have no such collision: the
exoTSO added in 260605 BC11–14 is a SMARTer oligo and lacks the Nextera anchor, so it
cannot be read as a code.

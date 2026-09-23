# Coverage tracks with junction arcs for the splicing loci (Fig. 7E, 7F)

Browser-style views of the splicing events in the reads: per locus the whole gene and a zoom on the event exon and its
neighbours, as coverage tracks for every replicate of control and treatment and as pseudobulk tracks per condition, the zoom
also with the splice junctions drawn as arcs that carry their fragment counts. Twenty-five loci: the 17 genes of the CHX top-10
lists of `../screen_splicing/` and the three CHX / NMD positive-control exons on control / CHX 1 / CHX 50 µg/mL tracks, the five
risdiplam positive-control exons on control / Ris 25 nM tracks. Fig. 7E shows AKAP8 and CLASRP (`zoom_arcs`, pseudobulk), Fig. 7F
FOXM1 and MADD (`zoom_arcs`, ris_pseudobulk).

The recipe is the benchmark's `gene_coverage_tracks.py` (one filled coverage track per BAM with the peak depth as its only y tick,
the canonical transcript below, 5' to 3' left to right, genomic and exon-only versions), with the gene model read from the GTF,
tracks coloured by condition, the zoom views, the four track sets and the junction arcs added.

## Data

The twelve wells `Hek_control_1/2/3`, `CHX_dose_1ug_mL-1/2/3`, `CHX_dose_50ug_mL-1/2/3`, `Ris_dose_25mM-1/2/3` of set
`bobseq_pe_native`, from the per-sample BAMs of the benchmark. Unit = fragments (MAPQ 255, a read pair counts once, not
deduplicated): the unit of the differential-splicing tables whose hits are drawn. Every track of a set is subsampled exactly
to the smallest member of its set (`results/track_depths.tsv`):

| set | tracks | matched to (fragments) |
|---|---|---|
| replicates | control 1, 2, 3; CHX 50 µg/mL 1, 2, 3 | 4,107,729 |
| pseudobulk | control, CHX 1 µg/mL, CHX 50 µg/mL (the three wells merged) | 13,315,023 |
| ris_replicates | control 1, 2, 3; Ris 25 nM 1, 2, 3 | 4,377,823 |
| ris_pseudobulk | control, Ris 25 nM | 16,865,099 |

Loci (`config/loci.tsv`, stage `loci`): the CHX top-10 loci of `../screen_splicing/results/top10_loci.tsv` (hits and the padded
non-hits RPL8, TMEM208 and NOP56, flagged), one per gene, with the top junction, its LSV and every member junction of that LSV
(from the row ids of S6); plus the eight positive controls of `src/positive_control_events.py` with their inclusion and skip
junctions and the event exon: CASP2 exon 9, SRSF3 exon 4 and TRA2B exon 2 (NMD-coupled, CHX track sets), FOXM1 exon 9, the
SLC25A17 poison exon, APLP2 exon 7, STRN3 exon 8 and the MADD exon-13 extension (risdiplam track sets).

## Views and rules (`src/gene_coverage_tracks.py`)

Coverage = depth per base from the aligned blocks of every alignment record (spliced gaps are not covered). Junctions = every N in
a CIGAR with >= 4 nt aligned on both sides on that read, counted once per read name. Views: `gene_with_introns` (the canonical
transcript span, introns to scale), `gene_exons_only` (exons concatenated, ticks = hg38 position of the exon starts), `zoom` and
`zoom_arcs` (every exon that carries a site of the LSV's junctions plus the event exon, widened to at least three exons with the
nearest canonical neighbours, terminal exons cut to 2 kb, ± 5 %; without and with arcs); `zoom_site_arcs` for CLASRP, a
446-nt window on its alternative 5' splice sites. Arcs (`zoom_arcs` only): one per junction of the junction catalogue that has
>= 2 fragments with both ends in view and >= 2 % of the strongest catalogue junction in that track's view; line width by
fragment count, the count printed at the apex; no exception for the LSV's own junctions. The catalogue is the junction set the
LSVs of the tables are defined on; `config/atlas_junctions_loci.tsv` holds its 3,963 junctions within 20 kb of the 25 loci, which
is everything the rule needs here. The event exon is drawn in red when it is not part of the canonical transcript.

## Results

`results/summary.tsv`: every LSV junction of every locus with its fragment count and share of the LSV on every track;
`results/values/<locus>_junctions_<set>.tsv`: every junction in the locus window, drawn or not; `results/gene_models.tsv` and
`results/windows.tsv`: the transcript drawn and the windows. On the pseudobulk tracks the top junction's share rises from control
through CHX 1 to CHX 50 for every hit (CLASRP 0.26, 0.43, 0.70; AKAP8 0.05, 0.23, 0.40); FOXM1 exon 9 and the SLC25A17 poison exon
are absent from control and appear under Ris 25 nM; the MADD distal 5' splice-site form goes down (5 of 9 to 1 of 15 exon-13
fragments), not up.

## Run

```
python3 run.py                                   # loci, bams (12 per-sample BAMs, filtering and exact subsampling, about 50 min), figures (about 30 min), summary
python3 run.py --only figures summary            # redraw from the cached track BAMs
python3 src/gene_coverage_tracks.py CLASRP AKAP8 --sets pseudobulk --views zoom_arcs
```
Needs the per-sample BAMs and the Ensembl 113 GTF (`../locations.py`), samtools, pysam, numpy, matplotlib. The per-base
coverage values are written to `results/values/<locus>_depth_<set>.tsv` on a run and are not part of the repository.

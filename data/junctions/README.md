# Per-sample splice-junction tables behind the junction panels (p5)

One file per sample and set, `<set>__<sample>_junctions.bed.gz`: BED6 with chromosome, intron start (0-based), intron end,
sample, molecules supporting the junction, strand ('?': the alignments carry no XS tag). Junctions were extracted from the
per-sample human BAMs of the benchmark after UMI-and-position deduplication (BOBseq 6-nt UMI, DRUG-seq and prime-seq
10 nt), MAPQ 255 molecules only, intron coordinates of every gapped alignment with at least 4-nt anchors on both sides; a
pair spanning a junction counts once. `sample_metadata.csv` gives, per sample and set, the method, cell line, condition,
read layout, records before and after deduplication and after the MAPQ filter, and the junction count.

Sets: drugseq and primeseq (50 nt), drugseq_native (52 nt), primeseq_native (50 nt), bobseq_50nt (read 2 truncated to
50 nt) and bobseq_pe_native (both mates, 2 x 150 nt); 103 sample-arms.

These are the tables the junction panels of the preprint were computed from. `benchmark/extract_junctions.sh` regenerates them from the
per-sample BAM store; `benchmark/scripts/junction_rarefaction*.py` read them from `<work>/data/junctions/`, so copy this
directory there (or point BOBSEQ_WORK at a directory containing it) to redraw the p5 panels without rerunning the alignments.

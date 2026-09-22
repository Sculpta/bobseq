#!/usr/bin/env python3
"""README.md of <fig set>/main_figure_panels/ for the tidied layout: the kept panels, values/, extra/, key numbers from the
panel sidecars. BM_COVSET=native|50nt."""

from settings import WORK
import os, sys, json, datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from palette import METH, COL  # method list + colours

P = WORK
SET = os.environ.get('BM_COVSET', 'native')
TAG = 'native' if SET == 'native' else 'matched'
D = f'{P}/{"preprint_figures_native" if SET == "native" else f"preprint_figures_{SET}"}/main_figure_panels'
J = lambda n: json.load(open(f'{D}/{n}'))
thr = J('p4_thresholds_numbers.json')
p6 = J(f'p6_picard_numbers_{TAG}.json')
p3 = J(f'p3_read_composition_{TAG}.json')
p3b = J(f'p3b_intronic_per_sample_{TAG}.json')
import statistics as st

setdesc = (
    'native read length: every method as sequenced, BOBseq 2 x 150 paired (both mates in the coverage panels, mate 1 for the pipeline metrics)'
    if SET == 'native'
    else 'read-length matched: every method one 50-nt single read per fragment (BOBseq read 2 truncated to 50 nt)'
)
L = [
    f'# Main-figure panels, {SET} set ({datetime.date.today()})',
    '',
    f'Set: {setdesc}. Samples: DRUG-seq 24 DMSO wells, prime-seq 8 HEK samples, BOBseq 18 HEK wells (Ris 500 nM excluded as input failures). One editable SVG per panel (svg.fonttype none, Arial), PNG preview beside it.',
    'Colours: '
    + ', '.join(f'{m} {COL[m]}' for m in METH)
    + ' (benchmark/scripts/palette.py). The raw values in values/ are the reference; replot.py restyles any panel from them.',
    '',
    '## Panels of the figure',
    '',
    f'- p0a_read_fate_{TAG}: reads into STAR surviving to mapped / unique / filtered (unique and exonic protein-coding), % of the well\'s reads into STAR; line = mean over the method\'s samples, ribbon = 95% t-interval of the mean, in-plot legend with n; per-sample values in values/.',
    f'- p1_molecules_vs_depth_{TAG}: molecules (UCI) per sample vs subsampled filtered reads, fig 3C style: thick line = mean over samples, ribbon = 95% t-interval of the mean (depths with >= 3 samples), linear axis, in-plot legend with n; per-sample values in values/.',
    f'- p1_molecules_vs_depth_log_{TAG}: the same panel with a log y axis (both versions kept); the ribbon is thinner than the line there because the interval is a few % of the mean.',
    f'- p2_genes_vs_depth_{TAG}: genes per sample vs subsampled filtered reads, ribosomal-protein genes (RPL/RPS/RPLP/MRPL/MRPS) not counted (column genes_B_noRP of rarefied.tsv, the count with them is kept beside it in values/); mean, ribbon and legend as p1.',
    f'- p3_read_composition_{TAG}: % of mapped reads (multimappers included) in the seven classes of the methods\' classification rules (mRNA, ribosomal-protein mRNA, exonic other biotype, intronic, intergenic, rRNA, mitochondrial), from the per-sample BAMs (composition_rules.py); the read-table version is in extra/ (*_readtable_*).',
    f'- p0c_reads_after_dedup_{TAG}: filtered reads per sample after de-duplication (molecules, UCI), one dot per sample, bar = median, log axis, no subsampling (only the de-duplicated count is a figure panel; the before/after pair p0cd is in extra/).',
    f'- p0d2_duplicate_rate_at_250k_{TAG}: duplicate rate per sample = 100 x (1 - molecules / filtered reads) with every sample subsampled to 250k filtered reads (matched depth; samples below 250k are not shown), one dot per sample, bar = median. The as-sequenced rate rises with depth and is not a figure panel; the paired version p0dd and p0d stay in extra/.',
    f'- p7a_read_length_{TAG}: aligned read length per method (BOBseq per mate), % of unique reads on canonical transcripts.',
    '- p7b_insert_length_bobseq: insert length of deduplicated BOBseq pairs on the mature transcript (both mates on the same canonical transcript; exon-aware), 18 benchmark wells pooled, per-well medians in the json (native only, the competitors are single-end).',
    f'- p7e_aligned_bases_total_{TAG}: total aligned bases per sample (sum of the aligned length over the unique reads on canonical transcripts, as sequenced; BOBseq both mates), filled = as sequenced, open = at 250k filtered reads (250k x the sample\'s aligned bases per fragment, which is depth-independent), bar = median (read number x mean aligned length).',
    f'- (extra/) p7f_aligned_bases_per_fragment_{TAG}: the total divided by the number of unique fragments; kept in extra/ only.',
    '- (extra/) p7c_quality_by_cycle: mean Phred per cycle for BOBseq read 1, BOBseq read 2 and a poly(dT) library read through the poly(A) (an internal poly(dT) BOBseq library); not a figure panel, it does not compare the benchmark methods.',
    '- (extra/) p7d_base_share_by_cycle: A and T share per cycle of the same reads; not a figure panel, the poly(A) block starts at a different cycle in every read so the per-cycle share smears.',
    f'- p5_junctions_vs_depth_min3_{TAG}: unique splice junctions with >= 3 supporting molecules vs deduplicated molecules per sample (per-sample junction BEDs of extract_junctions.sh, UMI-deduplicated MAPQ-255 molecules, pair spanning a junction = 1; rarefaction by binomial thinning), mean and 95% t-interval; other thresholds in junctions/all_junctions/.',
    f'- p0_samples_{TAG}: sample table (replicates, cell line, read layout per method), counts taken from the data.',
    f'- p3b_intronic_per_sample_{TAG}: intronic reads as % of mapped reads, one dot per sample, bar = median.',
    f'- p4i2_coverage_heatmaps_all_reads_1kb_noRP: {thr["genes"]} genes (>= 200 unique reads in every method, transcript >= 1 kb, {len(thr["rp_genes_removed"])} ribosomal-protein genes excluded), all unique reads, one row per gene scaled to its own peak, short to long.',
    f'- p4m2_coverage_vs_nt_from_polyA_mean_1kb_noRP: same genes, 200 reads per gene, coverage / gene mean, mean over genes, last 3,000 nt from the poly(A) site (25-nt running mean). Genes contributing at 1 / 2 / 3 kb: {thr["genes_contributing_at_nt"]["1000"]} / {thr["genes_contributing_at_nt"]["2000"]} / {thr["genes_contributing_at_nt"]["3000"]}.',
    f'- p6_picard_profile_{TAG}: Picard CollectRnaSeqMetrics normalised coverage (MINIMUM_LENGTH 1000, STRAND NONE, protein-coding transcripts without ribosomal-protein and MT genes), mean over samples with the 95% t-interval.',
    f'- p6b_picard_balance_{TAG}: 5\' to 3\' balance = 2 x coverage centroid per sample (1 = even); points + median (extra/ holds the box version).',
    '- legend_*: shared legends (method lines, method dots, composition classes).',
    '',
    '## Directories',
    '',
    '- values/: the raw numbers behind every kept panel, one TSV per panel (per sample, per gene or per bin).',
    '- extra/: every other panel and sidecar from the generators (p0b, p1/p2 strips and the other depth set, p4a-p4q coverage variants, p5 junction placeholders, p4/p5 number files). Regenerated by the pipeline, moved here by tidy_panels.py.',
    '',
    '## Key numbers',
    '',
    f'- Composition ({TAG}, the methods\' rules): rRNA % of mapped '
    + ', '.join(f'{m} {p3["per_method_pct_of_mapped"][m]["rRNA"]:.1f}' for m in METH)
    + '; mRNA (protein-coding, RP excluded) '
    + ', '.join(f'{m} {p3["per_method_pct_of_mapped"][m]["mRNA"]:.1f}' for m in METH)
    + '; intronic % per sample, median '
    + ', '.join(f'{m} {st.median(p3b[m].values()):.1f}' for m in METH)
    + '.',
    '- Coverage profile at 200 reads per gene (>= 1 kb, no RP): peak x gene mean '
    + ', '.join(f'{m} {thr["profile"][m]["peak"]:.2f}' for m in METH)
    + '; body 1-2 kb '
    + ', '.join(f'{m} {thr["profile"][m]["mean_1000_2000"]:.2f}' for m in METH)
    + '.',
    '- Picard: peak percentile / height '
    + ', '.join(f'{m} {p6["per_method"][m]["peak_percentile"]} / {p6["per_method"][m]["peak"]:.2f}' for m in METH)
    + '; balance mean (sd) '
    + ', '.join(
        f'{m} {p6["per_method"][m]["balance_mean"]:.3f} ({p6["per_method"][m]["balance_sd"]:.3f})' for m in METH
    )
    + '.',
    '',
    '## Generators (benchmark/)',
    '',
    '- sample_table_panel.py (p0), read_fate_panels.py (p0a), main_figure_panels.py (p1, p2 and the extra/ p3/p3b read-table versions, p4a-p4i, p5), composition_rules.py + composition_panels.py (p3, p3b), picard_custom.py + picard_unscaled.py (extra/ p4j-p4q), threshold_panels.py (p4i2/p4m2 threshold versions), picard_profile.py (p6, p6b; metrics cached in ../picard_metrics/), panel_legends.py, panel_values.py (values/), tidy_panels.py (extra/), validate_main_panels.py (VALIDATION_main_panels.md).',
    f'- The other set is in ../../preprint_figures_{"50nt" if SET == "native" else "native"}/main_figure_panels/ (BM_COVSET).',
]
open(f'{D}/README.md', 'w').write('\n'.join(L) + '\n')
print(SET, 'README written')

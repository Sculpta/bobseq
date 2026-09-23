"""Locations, track sets and frozen parameters of splicing_tracks (adapted from the benchmark's settings.py: the track sets are
the 24-plex conditions, and nothing else is hard-coded elsewhere).

Coverage and junction counts come from depth-matched BAMs built by prepare_bams.py from the per-sample BAMs of the
`bobseq_pe_native` wells (the per-sample BAM store of the benchmark, benchmark/build_per_sample_bams.sh): primary alignments
with MAPQ 255 (uniquely mapped), one BOBseq read pair = one fragment, NOT UMI-deduplicated -- the unit of the screen_splicing
differential-splicing tables whose hits are drawn here. Depth matching follows the benchmark: every track of a set is
subsampled to the same number of fragments (the smallest member of the set), so peak depths and junction counts are
comparable between the tracks of one figure.
"""
import os
import sys

CODE = os.path.dirname(os.path.abspath(__file__))                       # splicing_tracks/src
ANALYSIS = os.path.dirname(CODE)                                        # splicing_tracks
sys.path.insert(0, os.path.dirname(ANALYSIS))
from locations import GTF, SAMPLE_BAMS, table, open_text                # noqa: E402,F401
CONFIG = os.path.join(ANALYSIS, 'config')
RESULTS = os.path.join(ANALYSIS, 'results')
FIGURES = os.path.join(ANALYSIS, 'figures')
CACHE = os.path.join(ANALYSIS, '.cache')
BAMS = os.path.join(CACHE, 'bams')                                      # depth-matched track BAMs (rebuildable)

SCREEN_SPLICING = os.path.join(os.path.dirname(ANALYSIS), 'screen_splicing', 'results')   # top10_loci.tsv, ds_top10.tsv
JUNCTION_TABLE = 'S6'                                                   # S6_junction_counts_fragments.csv.gz: its row ids give the LSV members
RAW_BAMS = os.path.join(SAMPLE_BAMS, 'bobseq_pe_native')                # <well>.bam
SET = 'bobseq_pe_native'

# The wells behind the tracks (sample names of the sample sheet) and the track sets.
WELLS = {
    'control': ['Hek_control_1', 'Hek_control_2', 'Hek_control_3'],
    'chx_low': ['CHX_dose_1ug_mL-1', 'CHX_dose_1ug_mL-2', 'CHX_dose_1ug_mL-3'],
    'chx_high': ['CHX_dose_50ug_mL-1', 'CHX_dose_50ug_mL-2', 'CHX_dose_50ug_mL-3'],
    'ris_low': ['Ris_dose_25mM-1', 'Ris_dose_25mM-2', 'Ris_dose_25mM-3'],          # Ris 500 mM = the QC-failed wells, never drawn
}
COND_LABEL = {'control': 'control', 'chx_low': 'CHX 1 µg/mL', 'chx_high': 'CHX 50 µg/mL', 'ris_low': 'Ris 25 mM'}
# track id -> (label, condition, wells). '<x>replicates' = wells one by one, each subsampled to the smallest of the set;
# '<x>pseudobulk' = conditions, replicate BAMs merged, each subsampled to the smallest condition of the set. The CHX loci use
# replicates / pseudobulk (control, CHX 1, CHX 50); the risdiplam positive controls use ris_replicates / ris_pseudobulk
# (control, Ris 25 mM) -- config/loci.tsv column `sets` says which.
TRACK_SETS = {
    'replicates': [(f'{c}_{i}', f'{COND_LABEL[c]} {i}', c, [w]) for c in ('control', 'chx_high') for i, w in enumerate(WELLS[c], 1)],
    'pseudobulk': [(c, COND_LABEL[c], c, WELLS[c]) for c in ('control', 'chx_low', 'chx_high')],
    'ris_replicates': [(f'{c}_{i}', f'{COND_LABEL[c]} {i}', c, [w]) for c in ('control', 'ris_low') for i, w in enumerate(WELLS[c], 1)],
    'ris_pseudobulk': [(c, COND_LABEL[c], c, WELLS[c]) for c in ('control', 'ris_low')],
}

MAPQ = 255              # STAR: uniquely mapped
SEED = 0                # fragment subsampling
MIN_ANCHOR = 4          # nt aligned on both sides of a splice for a fragment to count for the junction (regtools -a 4, the rule of the junction tables)
MIN_ARC = 2             # junctions with fewer fragments than this in a track are not drawn
ARC_MIN_FRAC = 0.02     # ... nor junctions below this fraction of the strongest catalogue junction in the track's view, and only
                        # junctions of the catalogue (ATLAS_JUNCTIONS) are drawn at all; no exception for the LSV's own junctions
ATLAS_JUNCTIONS = os.path.join(CONFIG, 'atlas_junctions_loci.tsv')   # the junction catalogue the LSVs are defined on, for the loci drawn here
EXON_CAP = 2000         # terminal exons longer than this are cut to it in the zoom window (the benchmark rule)
ZOOM_PAD = 0.05         # zoom window padding, fraction of the window span
THREADS = int(os.environ.get('USV_THREADS', 4))

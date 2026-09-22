#!/usr/bin/env python3
"""The panels of the main figure: one list, used by tidy_panels.py (moves everything else to extra/) and by the
generators' save() (non-figure panels are written straight into extra/, so legacy panels no longer appear in the panel directory).
"""


def keep(tag):
    return {
        f'p0a_read_fate_{tag}',
        f'p1_molecules_vs_depth_{tag}',
        f'p1_molecules_vs_depth_log_{tag}',
        f'p2_genes_vs_depth_{tag}',
        f'p3_read_composition_{tag}',
        f'p3b_intronic_per_sample_{tag}',
        'p4i2_coverage_heatmaps_all_reads_1kb_noRP',
        'p4m2_coverage_vs_nt_from_polyA_mean_1kb_noRP',
        'p4_thresholds_numbers',
        f'p6_picard_profile_{tag}',
        f'p6b_picard_balance_{tag}',
        f'p6_picard_numbers_{tag}',
        f'p3_composition_rules_{tag}',
        f'p0_samples_{tag}',
        f'p0c_reads_after_dedup_{tag}',
        f'p0d2_duplicate_rate_at_250k_{tag}',
        f'p7a_read_length_{tag}',
        'p7b_insert_length_bobseq',
        f'p7e_aligned_bases_total_{tag}',
        f'p5_junctions_vs_depth_min3_{tag}',
        f'p5_annotated_junctions_vs_depth_min3_{tag}',
        f'p5_pooled_annotated_junctions_vs_depth_min3_{tag}',
        f'p5_pooled_all_junctions_vs_depth_min3_{tag}',
    }  # p7f (bases per fragment) stays in extra/: the per-fragment number is in p7a   # p7c_quality_by_cycle and p7d_base_share_by_cycle stay in extra/ (not a method comparison; the poly(A) block starts at a different cycle per read, so the per-cycle share smears)   # the single-quantity versions p0c / p0d / p0d2 live in extra/

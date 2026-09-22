#!/usr/bin/env python3
"""Curated figure tree of the preprint project: one main-figure dir, one supplement dir, one archive of dropped panels.
The generators keep writing into figures/_generator_output/<set>/... (their validated layout, reached through the preprint_figures_* symlinks); this script
copies each panel (svg + png + json, its values tables) to its place and rebuilds the three curated dirs from scratch on every run, so they never drift.
    main/benchmark, main/dataset_24plex, supplement/benchmark_native, supplement/benchmark_50nt, supplement/per_sample_24plex, archive_dropped/
Fails when a panel of a generator directory is not placed (new panel without a decision), when a placed panel is missing, or when a copy differs from its source.
Change a placement = edit PLACE below and rerun. Writes figures/README.md (index) and figures/MANIFEST_figures.tsv (destination, source, sha256).
"""

from settings import ROOT, PROJECT
import os, sys, glob, shutil, hashlib, datetime

FIG = f'{PROJECT}/figures'
B = f'{FIG}/_generator_output'
SRC = {
    'nat': f'{B}/native/main_figure_panels',
    'm50': f'{B}/50nt/main_figure_panels',
    'sup': f'{B}/native/supplement_per_sample_21',
    'mus': f'{B}/native/supplement_per_sample_mouse',
    'stk': f'{B}/native/supplement_stacked',
    'gex': f'{B}/native/gene_coverage_tracks',
}
PLACE = {
    'main/benchmark': dict(
        why='main figure, benchmark vs DRUG-seq and prime-seq (native read length)',
        panels=[
            ('nat', 'p0c_reads_after_dedup_native'),
            ('nat', 'p3_read_composition_native'),
            ('nat', 'p3b_intronic_per_sample_native'),
            ('nat', 'p4m2_coverage_vs_nt_from_polyA_mean_1kb_noRP'),
            ('nat', 'p5_junctions_vs_depth_min3_native'),
            ('nat', 'p5_annotated_junctions_vs_depth_min3_native'),
            ('nat', 'p5_pooled_annotated_junctions_vs_depth_min3_native'),
            ('nat', 'p5_pooled_all_junctions_vs_depth_min3_native'),
            ('nat', 'p6b_picard_balance_native'),
        ],
        values=[('nat', ('p0c_p0d_', 'p3_', 'p3b_', 'p4m2_', 'p5_', 'p6b_'))],
        files=[
            ('nat', 'p3_composition_rules_native.json'),
            ('nat', 'p4_thresholds_numbers.json'),
            ('nat', 'p6_picard_numbers_native.json'),
        ],
    ),
    'main/gene_examples': dict(
        why='main-figure gene examples on the pooled depth-matched tracks (8.0 M unique fragments per method), each with introns and exons only; full-coverage examples RELA, MLF2, BAG6, TUBG1, EEF2, MCM7, GIGYF1 and the prime-seq intronic examples GRB2, IRS4, FUS, PPP5C ; values/ = per-base depth per gene',
        panels=[
            ('gex', f'{g}_{v}')
            for g in ['RELA', 'MLF2', 'BAG6', 'TUBG1', 'EEF2', 'MCM7', 'GIGYF1', 'GRB2', 'IRS4', 'FUS', 'PPP5C']
            for v in ('with_introns', 'exons_only')
        ],
        values=[('gex', ('',))],
        files=[],
    ),
    'supplement/benchmark_native': dict(
        why='supplement, benchmark at native read length',
        panels=[
            ('nat', 'p0_samples_native'),
            ('nat', 'p0a_read_fate_native'),
            ('nat', 'p0d2_duplicate_rate_at_250k_native'),
            ('nat', 'p1_molecules_vs_depth_native'),
            ('nat', 'p1_molecules_vs_depth_log_native'),
            ('nat', 'p2_genes_vs_depth_native'),
            ('nat', 'p4i2_coverage_heatmaps_all_reads_1kb_noRP'),
            ('nat', 'p7e_aligned_bases_total_native'),
        ],
        values=[('nat', ('p0_samples_', 'p0a_', 'p0c_p0d_', 'p1_', 'p2_', 'p4i2_', 'p7e_p7f_'))],
        files=[('nat', 'p4_thresholds_numbers.json')],
    ),
    'supplement/benchmark_50nt': dict(
        why='supplement, the one 50-nt control kept: coverage vs distance from the poly(A) site with every method reduced to one 50-nt read, showing the coverage holds with shortened reads',
        panels=[('m50', 'p4m2_coverage_vs_nt_from_polyA_mean_1kb_noRP')],
        values=[('m50', ('p4m2_',))],
        files=[('m50', 'p4_thresholds_numbers.json')],
    ),
    'archive_dropped/benchmark_50nt': dict(
        why='the other 50-nt panels (only one 50-nt control is shown)',
        panels=[
            ('m50', x)
            for x in (
                'p0_samples_matched',
                'p0a_read_fate_matched',
                'p0c_reads_after_dedup_matched',
                'p0d2_duplicate_rate_at_250k_matched',
                'p1_molecules_vs_depth_matched',
                'p1_molecules_vs_depth_log_matched',
                'p2_genes_vs_depth_matched',
                'p3_read_composition_matched',
                'p3b_intronic_per_sample_matched',
                'p4i2_coverage_heatmaps_all_reads_1kb_noRP',
                'p5_junctions_vs_depth_min3_matched',
                'p6b_picard_balance_matched',
                'p7e_aligned_bases_total_matched',
            )
        ],
        values=[
            (
                'm50',
                ('p0_samples_', 'p0a_', 'p0c_p0d_', 'p1_', 'p2_', 'p3_', 'p3b_', 'p4i2_', 'p5_', 'p6b_', 'p7e_p7f_'),
            )
        ],
        files=[('m50', 'p3_composition_rules_matched.json'), ('m50', 'p6_picard_numbers_matched.json')],
    ),
    'supplement/per_sample_24plex': dict(
        why='supplement, per-sample panels of the 21 human BOBseq wells (mouse wells to be added)',
        panels=[('sup', 'S14_coverage_vs_nt_from_polyA_all24_native')],
        values=[('sup', ('polyA_profiles_all24_',))],
        files=[('sup', 'README_S14_all24.md')],
    ),
    'supplement/stacked_qc_all_samples': dict(
        why='supplement, the merged stacked QC figure: 14 metrics, narrow layout with staggered dots, columns = DRUG-seq, prime-seq, the BOBseq 24-plex by treatment (replicates stacked, line = mean), risdiplam 500 nM input failures, mouse wells in their own colour; matched-depth rows at 250k filtered reads',
        panels=[('stk', 'stacked_qc_all_samples')],
        values=[('stk', ('stacked_qc_all_',))],
        files=[],
    ),
    'covered_by_stacked_24plex/superseded_stacks': dict(
        why='the two separate stacked figures (24-plex only, benchmark only) merged into supplement/stacked_qc_all_samples; kept for reference',
        panels=[('stk', 'stacked_24plex_per_sample_qc'), ('stk', 'stacked_benchmark_native_per_sample_qc')],
        values=[('stk', ('stacked_24plex_', 'stacked_benchmark_'))],
        files=[],
    ),
    'covered_by_stacked_24plex/human': dict(
        why='single-metric per-sample panels whose content is now a row of supplement/stacked_24plex; kept for reference, not part of the supplement',
        panels=[
            ('sup', x)
            for x in (
                'S01_read_fate_native',
                'S02_reads_after_dedup_native',
                'S03_duplicate_rate_native',
                'S06_composition_native',
                'S07_intronic_native',
                'S08_read_length_native',
                'S10_aligned_bases_native',
                'S12_picard_balance_native',
            )
        ],
        values=[],
        files=[],
    ),
    'covered_by_stacked_24plex/mouse': dict(
        why='the same for the three mouse wells',
        panels=[
            ('mus', x + '_native_mouse')
            for x in (
                'S01_read_fate',
                'S02_reads_after_dedup',
                'S03_duplicate_rate',
                'S06_composition',
                'S07_intronic',
                'S08_read_length',
                'S10_aligned_bases',
                'S12_picard_balance',
            )
        ],
        values=[],
        files=[],
    ),
    'archive_dropped': dict(
        why='dropped from the preprint; last version kept for reference, still produced by the generators',
        panels=[
            ('nat', 'p7a_read_length_native'),
            ('nat', 'p7b_insert_length_bobseq'),
            ('nat', 'p6_picard_profile_native'),
            ('m50', 'p7a_read_length_matched'),
            ('m50', 'p6_picard_profile_matched'),
            ('sup', 'S03b_duplicate_rate_at_250k_native'),
            ('sup', 'S04_molecules_vs_depth_native'),
            ('sup', 'S04b_molecules_vs_depth_log_native'),
            ('sup', 'S05_genes_vs_depth_native'),
            ('sup', 'S11_picard_profile_native'),
            ('sup', 'S09_insert_length_native'),
            ('mus', 'S09_insert_length_native_mouse'),
            ('sup', 'S15_junctions_vs_depth_min3_native'),
        ],
        values=[
            ('nat', ('p7a_', 'p7b_', 'p6_picard_profile_')),
            ('m50', ('p7a_', 'p6_picard_profile_')),
            ('sup', ('junctions_',)),
        ],
        files=[],
    ),
}


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def find(src, name):
    for d in (SRC[src], f'{SRC[src]}/extra'):
        if os.path.exists(f'{d}/{name}'):
            return f'{d}/{name}'


errors = []
man = []
placed = {}
for dest, spec in PLACE.items():
    for src, stem in spec['panels']:
        placed.setdefault((src, stem), []).append(dest)
for k, v in placed.items():
    if len(v) > 1:
        errors.append(f'{k[1]} placed twice: {v}')
for src, d in SRC.items():
    for f in sorted(glob.glob(f'{d}/*.svg')):
        stem = os.path.basename(f)[:-4]
        if stem.startswith('legend_'):
            continue
        if (src, stem) not in placed:
            errors.append(f'not placed: {stem} ({d})')
for src, stem in placed:
    for ext in ('svg', 'png'):
        if not os.path.exists(f'{SRC[src]}/{stem}.{ext}'):
            errors.append(f'missing in the generator output: {stem}.{ext}')
if errors:
    sys.exit('organize_figures: nothing written\n  ' + '\n  '.join(errors))
for top in (
    'main',
    'supplement',
    'archive_dropped',
    'covered_by_stacked_24plex',
    'covered_by_all24_polyA_panel',
    'legends',
):  # names no longer produced are listed too, so a stale copy of such a folder is removed
    if os.path.isdir(f'{FIG}/{top}'):
        shutil.rmtree(f'{FIG}/{top}')


def put(srcpath, destdir):
    os.makedirs(destdir, exist_ok=True)
    t = f'{destdir}/{os.path.basename(srcpath)}'
    shutil.copy2(srcpath, t)
    if sha(t) != sha(srcpath):
        sys.exit(f'copy differs: {t}')
    man.append((os.path.relpath(t, FIG), os.path.relpath(srcpath, FIG), sha(t)))


for dest, spec in PLACE.items():
    for src, stem in spec['panels']:
        for ext in ('svg', 'png', 'json'):
            p = f'{SRC[src]}/{stem}.{ext}'
            if os.path.exists(p):
                put(p, f'{FIG}/{dest}')
    for src, name in spec['files']:
        p = find(src, name)
        if p:
            put(p, f'{FIG}/{dest}')
        else:
            sys.exit(f'number file missing: {name}')
    for src, prefixes in spec['values']:
        for p in sorted(glob.glob(f'{SRC[src]}/values/*.tsv')):
            if any(os.path.basename(p).startswith(x) for x in prefixes):
                put(p, f'{FIG}/{dest}/values')
for p in sorted(glob.glob(f"{SRC['nat']}/legend_*")):
    put(p, f'{FIG}/legends')
# the methods text of the panels is maintained in the repository (figures/METHODS_panels.md), not generated
_methods = os.path.join(ROOT, 'figures', 'METHODS_panels.md')
if os.path.exists(_methods):
    shutil.copy2(_methods, f'{FIG}/METHODS_panels.md')
elif os.path.exists(f"{SRC['nat']}/METHODS_main_figure_panels.md"):
    shutil.copy2(f"{SRC['nat']}/METHODS_main_figure_panels.md", f'{FIG}/METHODS_panels.md')
with open(f'{FIG}/MANIFEST_figures.tsv', 'w') as f:
    f.write('file\tsource\tsha256\n')
    [f.write('\t'.join(r) + '\n') for r in man]
L = [
    f'# Figures of the BOBseq preprint ({datetime.date.today()})',
    '',
    'Curated by `benchmark/scripts/organize_figures.py` from the generator output; do not edit files here, rerun the generator and then the organizer. ',
    '',
    '| directory | what | panels |',
    '|---|---|---|',
]
for dest, spec in PLACE.items():
    L.append(f"| `{dest}/` | {spec['why']} | " + ', '.join(s for _, s in spec['panels']) + ' |')
L += [
    '| `legends/` | shared method and composition legends (SVG + PNG) | |',
    '| `_generator_output/` | what the generator scripts write (both sets, extra/ with legacy panels, per-directory README, PANEL_SOURCES, validation reports); the curated dirs above are copies of it | |',
    '',
    'Each panel is one editable SVG (text as text) plus a PNG; `values/` beside the panels holds the raw values per panel for replotting; `METHODS_panels.md` gives the definition of every panel; `MANIFEST_figures.tsv` lists every curated file with its source and sha256.',
    'Terms: filtered reads = uniquely mapped reads on an exon of a protein-coding gene; molecules = unique (UMI, gene) pairs after the over-use filter and collision correction (UCI in the tables).',
    '',
]
open(f'{FIG}/README.md', 'w').write('\n'.join(L))
print(
    f'organize_figures: {len(man)} files in the curated tree; '
    + '; '.join(f"{d} {len(s['panels'])} panels" for d, s in PLACE.items())
)

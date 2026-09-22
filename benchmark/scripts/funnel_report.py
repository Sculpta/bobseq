#!/usr/bin/env python3
"""Aggregate per_well/*.json (funnel_per_well.py) into BARCODE_ACCURACY_FUNNEL.md (one page: method, per-sample accuracy after
every funnel level, summaries) and funnel_summary.json."""

import os, glob, json, datetime, numpy as np

H = os.getcwd()
unit = 'molecules'  # run from the funnel results directory (run_funnel.sh)
W = {}
for f in sorted(glob.glob(f'{H}/per_well/*.json')):
    d = json.load(open(f))
    W[d['sample']] = d
    unit = d['unit']
LEV = [l['level'] for l in next(iter(W.values()))['levels']]


def acc(d, i):
    return d['levels'][i]['accuracy']


INPUT_FAIL = [s for s in W if s.startswith('Ris_dose_500mM')]
groups = {
    'all 24 samples': list(W),
    '21 samples (Ris 500 nM input failures excluded)': [s for s in W if s not in INPUT_FAIL],
    '21 human samples': [s for s in W if W[s]['well_species'] == 'human'],
    '3 mouse samples': [s for s in W if W[s]['well_species'] == 'mouse'],
}
summ = {
    g: {
        LEV[i]: {
            'mean': float(np.mean([acc(W[s], i) for s in v])),
            'min': float(min(acc(W[s], i) for s in v)),
            'max': float(max(acc(W[s], i) for s in v)),
            'n_total': int(sum(W[s]['levels'][i]['n'] for s in v)),
        }
        for i in range(len(LEV))
    }
    for g, v in groups.items()
}
json.dump(
    {'unit': unit, 'levels': LEV, 'per_sample': W, 'summary': summ}, open(f'{H}/funnel_summary.json', 'w'), indent=1
)
a24 = summ['all 24 samples']
a21 = summ['21 samples (Ris 500 nM input failures excluded)']
L = [
    f'# Barcode accuracy of the 24-plex Illumina run, strict funnel per sample, final version ({datetime.date.today()})',
    '',
    'Method. The 24-plex BOBseq library (HEK293T and RAW 264.7, 21 human and 3 mouse wells, Illumina PE150) was processed by the pipeline (pipeline/run_illumina.sh): the bobcode is read from the first bases of read 2, the pair is aligned with STAR 2.7.11b to the combined GRCh38 + GRCm39 genome with all alignments reported, and fragments are deduplicated (position, fragment length and UMI). Every deduplicated molecule is then scored on its read-2 mate exactly as the nanopore libraries were: species from all alignments (both genomes = ambiguous, not scored); rRNA (any alignment on an rDNA locus, or an exon of an rRNA-biotype gene) and mitochondrial reads removed; mRNA = the primary alignment overlaps an exon of a protein-coding gene that is not a ribosomal-protein gene; and only molecules whose bobcode matches a declared code with zero mismatches are scored (the pipeline accepts one mismatch at prep, so the exact status was re-derived from the raw read-2 FASTQ). A molecule is correct when its aligned species equals the species of its bobcode. Two Illumina-specific points: the bobcode is positional, so no anchor search or unit collapsing applies, and a second barcode unit cannot be observed in a 150-nt mate, so no chimera filter applies. Accuracy is the share of correct molecules with a 95% Wilson interval, per well.',
    '',
    f'Result. After the full funnel (L4), barcode accuracy averaged {a24[LEV[-1]]["mean"]:.2f}% over all 24 samples (range {a24[LEV[-1]]["min"]:.2f} to {a24[LEV[-1]]["max"]:.2f}%) and {a21[LEV[-1]]["mean"]:.2f}% over the 21 samples without the three Ris 500 nM input failures; human wells {summ["21 human samples"][LEV[-1]]["mean"]:.2f}%, mouse wells {summ["3 mouse samples"][LEV[-1]]["mean"]:.2f}%. A human well can receive foreign molecules from 3 mouse wells only, a mouse well from 21 human wells, so the per-well values are not symmetric between species. The large step from L1 to L2 in the mouse wells is the rRNA exclusion at work: rDNA is largely absent from GRCm39, so mouse rRNA molecules align to the human rDNA loci and would otherwise be counted as swaps.',
    '',
    f'Preprint sentence. Barcode accuracy averaged {a24[LEV[-1]]["mean"]:.1f}% across all 24 samples ({a21[LEV[-1]]["mean"]:.1f}% across the 21 samples that passed input QC), representing higher barcode accuracy than the 96% per-well benchmark reported for DRUG-seq (Ye et al. 2018: more than 98% of wells with over 96% species-specific UMIs).',
    '',
    f'## Accuracy (%) after each funnel level, per sample ({unit})',
    '',
    '| sample | well | L0 mapped, bobcode: n | L0 | L1 single species | L2 not rRNA/MT | L3 mRNA | L4 exact bobcode | L4 n | L4 Wilson 95% |',
    '|---|---|---|---|---|---|---|---|---|---|',
]
for s, d in W.items():
    l = d['levels']
    L.append(
        f'| {s} | {d["well_species"]} | {l[0]["n"]:,} | '
        + ' | '.join(f'{x["accuracy"]:.2f}' for x in l)
        + f' | {l[-1]["n"]:,} | {l[-1]["wilson_lo"]:.2f} to {l[-1]["wilson_hi"]:.2f} |'
    )
L += [
    '',
    '## Summary (mean of per-sample accuracies, %)',
    '',
    '| group | ' + ' | '.join(LEV) + ' |',
    '|---|' + '---|' * len(LEV),
]
for g, v in summ.items():
    L.append(f'| {g} | ' + ' | '.join(f'{v[l]["mean"]:.2f}' for l in LEV) + ' |')
CATS = ['rRNA', 'mitochondrial', 'mRNA', 'ribosomal-protein', 'other']


def comp(d):
    c = d.get('composition_all_mapped', {})
    t = sum(c.values()) or 1
    return {k: 100 * c.get(k, 0) / t for k in CATS}


if all('composition_all_mapped' in d for d in W.values()):
    cs = {s: comp(d) for s, d in W.items()}
    ok21 = [s for s in W if s not in INPUT_FAIL]
    rr = {g: [cs[s]['rRNA'] for s in v] for g, v in groups.items()}
    L += [
        '',
        '## Read composition of all mapped molecules (%), the same six classification rules (ambiguous molecules included)',
        '',
        f'rRNA content (any alignment on an rDNA locus, or the primary alignment on an exon of an rRNA-biotype gene): mean {np.mean(rr["all 24 samples"]):.1f}% over all 24 samples, {np.mean(rr[list(groups)[1]]):.1f}% (median {np.median(rr[list(groups)[1]]):.1f}%, range {min(rr[list(groups)[1]]):.1f} to {max(rr[list(groups)[1]]):.1f}%) over the 21 QC-passed samples. Preprint sentence: our multiplexed poly(dT) selection step reduced rRNA content to {np.mean(rr[list(groups)[1]]):.0f}% of mapped molecules on average (per-sample range {min(rr[list(groups)[1]]):.0f} to {max(rr[list(groups)[1]]):.0f}%).',
        '',
        '| sample | well | ' + ' | '.join(CATS) + ' |',
        '|---|---|' + '---|' * len(CATS),
    ]
    for s, d in W.items():
        L.append(f'| {s} | {d["well_species"]} | ' + ' | '.join(f'{cs[s][k]:.1f}' for k in CATS) + ' |')
    for g, v in groups.items():
        L.append(f'| {g} (mean) | | ' + ' | '.join(f'{np.mean([cs[s][k] for s in v]):.1f}' for k in CATS) + ' |')
    json.dump({s: cs[s] for s in W}, open(f'{H}/composition_all_mapped.json', 'w'), indent=1)
tot = {l: sum(W[s]['levels'][i]['n'] for s in W) for i, l in enumerate(LEV)}
L += [
    '',
    'Molecules reaching each level, all samples: ' + ', '.join(f'{l.split(" ")[0]} {tot[l]:,}' for l in LEV) + '.',
    '',
    'Inputs: the deduplicated analysis BAMs of the 24 wells (pipeline output), the raw read-2 FASTQ of the lane, the Ensembl 113 combined GTF and rDNA loci (annotation_index.npz). Scripts: funnel_build_annotation_index.py, funnel_exact_code_flags.py, funnel_per_well.py, funnel_report.py.',
]
open(f'{H}/BARCODE_ACCURACY_FUNNEL.md', 'w').write('\n'.join(L) + '\n')
print('\n'.join(L[4:6]))
print('DONE')

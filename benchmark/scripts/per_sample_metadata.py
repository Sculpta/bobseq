#!/usr/bin/env python3
"""metadata.tsv of the per-sample BAM store (benchmark_uniform/per_sample_bams): one row per sample and set with the BAM path
(relative to the store), read layout, cell line, condition, replicate, whether the sample is in the benchmark, and the source of
the raw data. Read by the panel scripts (competitor sample lists) and by extract_junctions.sh.
    per_sample_metadata.py
"""
import csv
import json
import os
import re

from settings import CONFIG, RUN_JSON, UNIFORM

B = f'{UNIFORM}/per_sample_bams'
rows = []


def add(method, sample, layout50, layout_native, **kw):
    for d, layout, set_name in ((method, layout50, 'matched_50nt'), (method + '_native', layout_native, 'native')):
        if os.path.exists(f'{B}/{d}/{sample}.bam'):
            rows.append(dict(method=method, set=set_name, sample=sample, bam=f'{d}/{sample}.bam', read_layout=layout, **kw))


# DRUG-seq: the 24 DMSO wells of the plate, named by plate well and barcode
dmso = set(open(f'{CONFIG}/drugseq_dmso_wells.txt').read().split())
plate = {r['well_index']: r for r in csv.DictReader(open(f'{CONFIG}/drugseq_plate_VH02001704.csv'))}
for bc in sorted(dmso):
    r = plate[bc]
    add('drugseq', f"DMSO_{r['plate_well']}_{bc}", '50 nt, 1 read (truncated)', '52 nt, 1 read', cell_line='U-2 OS',
        condition='DMSO vehicle, 24 h', replicate=r['plate_well'], in_benchmark='yes',
        source='GSE176150 / SRR14730306, plate VH02001704_S4 (Li et al. 2022)', bobcode='')
# prime-seq: the eight samples, named as in the read names (underscore removed)
ena = {}
for r in csv.DictReader(open(f'{CONFIG}/primeseq_samples.tsv'), delimiter='\t'):
    ena.setdefault(r['sample'], []).append(r['run_accession'])
for s in sorted(ena):
    add('primeseq', s.replace('_', ''), '50 nt, 1 read', '50 nt, 1 read (= the matched set)', cell_line='HEK293T',
        condition='RNA-extraction set (Janjic et al. 2022, Fig. 4), Incubation + ProtK, 10,000 cells', replicate=s, in_benchmark='yes',
        source=f"E-MTAB-10142 / {', '.join(ena[s])}", bobcode='')
# BOBseq: every human well of the 24-plex; in_benchmark = the 18 wells outside the exclusion table
cfg = json.load(open(RUN_JSON))
excluded = {l.split('\t')[0]: l.split('\t')[2] for l in open(f'{CONFIG}/units_excluded_24plex.tsv').read().splitlines()[1:] if l.strip()}
run = cfg['run']


def safe(n):
    return re.sub(r'_+', '_', ''.join(c if (c.isalnum() or c in '-_.') else '_' for c in n)).strip('_')


for code, lab in cfg['bobcode_labels'].items():
    if cfg['tso_species_map'][code] != 'human':
        continue
    cond, rep = lab.rsplit('-', 1) if '-' in lab else lab.rsplit(' ', 1)
    s = safe(lab)
    inb = 'no (' + excluded[code.upper()].split(':')[0] + ')' if code.upper() in excluded else 'yes'
    if os.path.exists(f'{B}/bobseq_50nt/{s}.bam'):
        rows.append(dict(method='bobseq', set='matched_50nt', sample=s, bam=f'bobseq_50nt/{s}.bam', read_layout='50 nt, read 2 only (truncated)',
                         cell_line='HEK', condition=cond, replicate=rep, in_benchmark=inb, source=f'{run}, bobcode {code.upper()}', bobcode=code.upper()))
    if os.path.exists(f'{B}/bobseq_pe_native/{s}.bam'):
        rows.append(dict(method='bobseq', set='native', sample=s, bam=f'bobseq_pe_native/{s}.bam', read_layout='2 x 150 nt, paired',
                         cell_line='HEK', condition=cond, replicate=rep, in_benchmark=inb, source=f'{run}, bobcode {code.upper()}', bobcode=code.upper()))
cols = ['method', 'set', 'sample', 'bam', 'read_layout', 'cell_line', 'condition', 'replicate', 'in_benchmark', 'source', 'bobcode']
with open(f'{B}/metadata.tsv', 'w') as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter='\t')
    w.writeheader()
    w.writerows(rows)
import collections

print(len(rows), 'rows', dict(collections.Counter((r['method'], r['set']) for r in rows)))

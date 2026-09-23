#!/usr/bin/env python3
"""Stage summary: the numbers behind the figures, one row per locus and LSV junction -> results/summary.tsv.
For every locus: the junctions of its LSV (or the positive control's inclusion / skip junctions) with the fragment count
on every track of both sets (from results/values/<locus>_junctions_<set>.tsv), plus the junction's share of the LSV's
fragments per track (count / sum over the LSV's junctions on that track -- a PSI in the sense of the junction tables, on the
depth-matched fragments). The top junction of the screen_splicing table is flagged.
"""
import csv, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from settings import CONFIG, RESULTS

VALUES = os.path.join(RESULTS, 'values')
loci = list(csv.DictReader(open(os.path.join(CONFIG, 'loci.tsv')), delimiter='\t'))
sets = sorted({r['set'] for r in csv.DictReader(open(os.path.join(RESULTS, 'track_depths.tsv')), delimiter='\t')})
rows, tracks = [], {}
for loc in loci:
    members = [m.split(':')[1] for m in loc['members'].split(';')]            # "start-end" (0-based intron start)
    top = loc['junction'].split(':')[1]
    counts = {}
    for st in sets:
        f = os.path.join(VALUES, f"{loc['locus']}_junctions_{st}.tsv")
        if not os.path.exists(f):
            continue
        for r in csv.DictReader(open(f), delimiter='\t'):
            key = f"{int(r['intron_start_1based']) - 1}-{r['intron_end']}"
            for t in [c for c in r if c not in ('chrom', 'intron_start_1based', 'intron_end', 'length')]:
                counts[(key, st, t)] = int(r[t]); tracks.setdefault(st, []).append(t) if t not in tracks.get(st, []) else None
    for m in members:
        row = dict(locus=loc['locus'], dir=loc['dir'], source=loc['source'], contrast=loc['contrast'], lsv=loc['lsv'] or 'positive control', junction=f"{loc['junction'].split(':')[0]}:{m}:{loc['strand']}", top=int(m == top))
        for st in sets:
            for t in tracks.get(st, []):
                n = counts.get((m, st, t), 0); tot = sum(counts.get((k, st, t), 0) for k in members)
                row[f'{st}:{t}'] = n; row[f'{st}:{t}:share'] = round(n / tot, 3) if tot else ''
        rows.append(row)
cols = []
for r in rows:
    cols += [c for c in r if c not in cols]
for r in rows:
    for c in cols:
        r.setdefault(c, '')
with open(os.path.join(RESULTS, 'summary.tsv'), 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in rows:
        f.write('\t'.join(str(r[c]) for c in cols) + '\n')
print(f"{len(rows)} LSV junctions over {len(loci)} loci -> results/summary.tsv")

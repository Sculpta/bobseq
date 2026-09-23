#!/usr/bin/env python3
"""Stage loci: the loci to draw -> config/loci.tsv.

  * the CHX differential-splicing top-10 lists of ../screen_splicing (contrasts chx_dose and chx_high; hits and the padded
    non-hits alike, flagged), one locus per gene, with the top junction, its LSV and every member junction of that LSV
    (from the row ids of the junction table -- the sites that define the zoom window);
  * the eight positive controls of positive_control_events.py with their inclusion and skip junctions and the event exon:
    the three CHX / NMD exons (CASP2 exon 9, SRSF3 exon 4, TRA2B exon 2; drawn on the CHX track sets) and the five
    risdiplam exons (FOXM1 exon 9, SLC25A17 poison exon, APLP2 exon 7, STRN3 exon 8, MADD exon 13 extension; drawn on the
    control / Ris 25 mM track sets). Column `sets` names the track sets of a locus, `dir` its figure directory.
Junction coordinates are Ensembl contig, 0-based half-open intron [start, end); exons likewise.
"""
import os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from settings import CONFIG, SCREEN_SPLICING, JUNCTION_TABLE, table, open_text
import positive_control_events as pc

rows = []
# --- screen_splicing CHX top-10 loci ---------------------------------------------------------------------------------
import csv
members = {}
with open_text(table(JUNCTION_TABLE)) as fh:
    next(fh)
    for line in fh:
        lsv, jct = line.split(',', 1)[0].split('|')
        members.setdefault(lsv, []).append(jct)
seen = {}
for r in csv.DictReader(open(os.path.join(SCREEN_SPLICING, 'top10_loci.tsv')), delimiter='\t'):
    if r['contrast'] not in ('chx_dose', 'chx_high'):
        continue
    g = r['gene']
    jct = r['junction_id'].split('|')[1]                      # "19:45068058-45068419:+"
    entry = dict(locus=g, gene=g, source='screen_splicing', contrast=r['contrast'], rank=int(r['rank']), hit=int(r['hit']),
                 dpsi=float(r['dpsi']), padj=float(r['padj']), lsv=r['lsv'], strand=r['strand'],
                 junction=jct, members=';'.join(members[r['lsv']]), event_exon='', event_label='', event_transcript='')
    if g in seen:                                             # the same gene in both CHX lists: keep the better-ranked entry, note both
        seen[g]['contrast'] += f";{r['contrast']}"; seen[g]['rank'] = min(seen[g]['rank'], int(r['rank'])); seen[g]['hit'] = max(seen[g]['hit'], int(r['hit']))
        continue
    seen[g] = entry; rows.append(entry)
# --- positive controls -----------------------------------------------------------------------------------------------
for g, contrast in (('CASP2', 'chx_high'), ('SRSF3', 'chx_high'), ('TRA2B', 'chx_high'),
                    ('FOXM1', 'ris_low'), ('SLC25A17', 'ris_low'), ('APLP2', 'ris_low'), ('STRN3', 'ris_low'), ('MADD', 'ris_low')):
    ev, ge = pc.EVENTS[g], pc.GENES[g]
    jcts = [f"{ge['chrom']}:{s}-{e}:{ge['strand']}" for _, (s, e) in ev['incl'] + ev['skip']]
    ex0, ex1 = ev['exon'][0] - 1, ev['exon'][1]               # frozen 1-based inclusive -> 0-based half-open
    m = re.search(r'[A-Z0-9]+-\d{3}', ev['note'])
    rows.append(dict(locus=g, gene=g, source='positive_controls', contrast=contrast, rank=0, hit=1, dpsi=float('nan'), padj=float('nan'),
                     lsv='', strand=ge['strand'], junction=jcts[0], members=';'.join(jcts),
                     event_exon=f"{ge['chrom']}:{ex0}-{ex1}", event_label=ev['label'], event_transcript=m.group(0) if m else ''))
# a tighter window on the alternative splice site itself (view zoom_site_arcs), 1-based inclusive, per locus: CLASRP, the
# alternative 5'ss pair 45,068,054 / 45,068,058 and the 45,068,120 site, from the middle exon (45,068,015-45,068,058) to the
# shared acceptor at 45,068,419 and 20 nt into that exon
SITE_WINDOW = {'CLASRP': '19:45067995-45068440'}
# figure directory per locus = gene + source tag: _chx_d (chx_dose list), _chx_h (chx_high list),
# _chx_dh (both lists), _chx_pc (CHX / NMD positive control)
TAG = {'chx_dose': '_chx_d', 'chx_high': '_chx_h', 'chx_dose;chx_high': '_chx_dh', 'chx_high;chx_dose': '_chx_dh'}
for r in rows:
    ris = r['contrast'] == 'ris_low'
    r['dir'] = r['locus'] + (('_ris_pc' if ris else '_chx_pc') if r['source'] == 'positive_controls' else TAG[r['contrast']])
    r['sets'] = 'ris_replicates;ris_pseudobulk' if ris else 'replicates;pseudobulk'
    r['site_window'] = SITE_WINDOW.get(r['locus'], '')
os.makedirs(CONFIG, exist_ok=True)
cols = ['locus', 'dir', 'sets', 'gene', 'source', 'contrast', 'rank', 'hit', 'dpsi', 'padj', 'lsv', 'strand', 'junction', 'members', 'event_exon', 'event_label', 'event_transcript', 'site_window']
with open(os.path.join(CONFIG, 'loci.tsv'), 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in rows:
        f.write('\t'.join(str(r[c]) for c in cols) + '\n')
print(f"{len(rows)} loci -> config/loci.tsv: " + ', '.join(r['locus'] + ('' if r['hit'] else ' (ns)') for r in rows))

#!/usr/bin/env python3
"""Mouse wells: per-well tables from their arm-path read tables, same rules as supp21_tables.py for the Ris 500 nM wells, species = mouse.
perwell_mapped: reads.tsv rows per well without an N in the UMI (mapped primary reads, any contig), unique = NH 'U'
slice:          region E, NH U, MOUSE contigs -> <arm>/stats_E_U_mouse3/wupg.tsv;  basic/reads_per_well.tsv, basic/perwell_B.tsv (umi_b_lib.b_stats, genes_noRP added)
"""

from settings import WORK
import os, sys, collections

H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, H)
from umi_b_lib import b_stats

P = WORK
ARM = f'{P}/benchmark_uniform/bob57_24plex_pe_native'
SP = f'{ARM}/mouse3'
ST = f'{ARM}/stats_E_U_mouse3'
os.makedirs(f'{ST}/basic', exist_ok=True)
WELLS = {'TAAGACG': 'RAW_control_1', 'GTCCATC': 'RAW_control_2', 'GATGACAGTAT': 'RAW_control_3'}
mp = collections.Counter()
uq = collections.Counter()
dr = collections.Counter()
nsl = collections.Counter()
spc = collections.Counter()
wupg = open(f'{ST}/wupg.tsv', 'w')
for n in WELLS.values():
    for l in open(f'{SP}/{n}/reads.tsv'):
        f = l.rstrip('\n').split('\t')
        if f[0] not in WELLS:
            continue
        if 'N' in f[1]:
            dr[f[0]] += 1
            continue
        mp[f[0]] += 1
        uq[f[0]] += f[4] == 'U'
        spc[(f[0], f[2][:5], f[5])] += 1
        if f[5] == 'E' and f[4] == 'U' and f[2].startswith('MOUSE_'):
            wupg.write('\t'.join((f[0], f[1], f[2], f[3], f[6], f[7])) + '\n')
            nsl[f[0]] += 1
wupg.close()
with open(f'{SP}/perwell_mapped.tsv', 'w') as f:
    [f.write(f'{w}\t{mp[w]}\t{uq[w]}\n') for w in WELLS]
with open(f'{SP}/n_umi_dropped_per_well.tsv', 'w') as f:
    [f.write(f'{w}\t{dr[w]}\n') for w in WELLS]
with open(f'{SP}/region_by_species.tsv', 'w') as f:
    f.write('well\tspecies\tregion\trows\n')
    [f.write(f'{w}\t{s}\t{r}\t{v}\n') for (w, s, r), v in sorted(spc.items())]
with open(f'{ST}/basic/reads_per_well.tsv', 'w') as f:
    [f.write(f'{w}\t{nsl[w]}\n') for w in sorted(WELLS)]
Um = collections.defaultdict(list)
Gm = collections.defaultdict(list)
for l in open(f'{ST}/wupg.tsv'):
    f = l.rstrip('\n').split('\t')
    Um[f[0]].append(f[1])
    Gm[f[0]].append(f[4])
with open(f'{ST}/basic/perwell_B.tsv', 'w') as out:
    out.write(
        'well\tB_raw\tB_kept\tumis_excluded\tB_corr\tgenes_B\tgenes_saturated\ttags_seen\ttags_homopolymer\ttags_overused\tgenes_B_noRP\n'
    )
    for w in sorted(Um):
        b = b_stats(Um[w], Gm[w])
        out.write(
            f"{w}\t{b['B_raw']}\t{b['B_kept']}\t{b['excluded']}\t{b['B_corr']}\t{b['genes']}\t{b['saturated']}\t{b['tags']}\t{b['hp']}\t{b['ov']}\t{b['genes_noRP']}\n"
        )
        print(
            WELLS[w],
            'mapped',
            mp[w],
            'unique',
            uq[w],
            'filtered',
            len(Um[w]),
            'UCI',
            b['B_corr'],
            'genes noRP',
            b['genes_noRP'],
            flush=True,
        )
print('MOUSE3 TABLES DONE')

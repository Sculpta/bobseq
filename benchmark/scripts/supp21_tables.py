#!/usr/bin/env python3
"""Supplement (21 human wells): per-well tables for the three risdiplam 500 nM wells, from their arm-path read tables (build_read_table.sh), with the
rules the benchmark uses for the 18 wells. Usage: supp21_tables.py native|50nt
  perwell_mapped: reads.tsv rows per well without an N in the UMI (mapped primary reads), unique = NH 'U'
  slice:          region E, NH U, human contigs -> <arm>/stats_E_U_supp21/wupg.tsv (well umi contig pos gene strand)  (as slice.sh)
  rarefied.tsv:   rarefy.py <arm> E_U_supp21;  basic/perwell_B.tsv + reads_per_well.tsv: umi_b_lib.b_stats per well
  50 nt only:     per-well BAMs like per_sample_bams/bobseq_50nt (primary mapped records on human contigs, bare names, sorted, indexed)
                  -> positions_canonical_50nt/bobseq-ris-500mM-<k>.npz (build_positions.py), composition rules (matched), Picard (matched).
"""

from settings import CODE, WORK, PROJECT, COVERAGE, REFFLAT
import os, sys, json, subprocess, collections, glob, numpy as np, pysam

H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, H)
from umi_b_lib import b_stats

SET = sys.argv[1]
P = WORK
U = f'{P}/benchmark_uniform'
PR = PROJECT
ARM = f'{U}/bob57_24plex_pe_native' if SET == 'native' else f'{U}/bob57_24plex_pe_hs'
SP = f'{ARM}/supp21'
ST = f'{ARM}/stats_E_U_supp21'
os.makedirs(f'{ST}/basic', exist_ok=True)
WELLS = {'CTCTTAC': 'Ris_dose_500mM-1', 'ACACCAT': 'Ris_dose_500mM-2', 'GCGAATG': 'Ris_dose_500mM-3'}
tables = [f'{SP}/{n}/reads.tsv' for n in WELLS.values()] if SET == 'native' else [f'{SP}/reads.tsv']
mp = collections.Counter()
uq = collections.Counter()
dr = collections.Counter()
wupg = open(f'{ST}/wupg.tsv', 'w')
nsl = collections.Counter()
for t in tables:
    for l in open(t):
        f = l.rstrip('\n').split('\t')
        if f[0] not in WELLS:
            continue
        if 'N' in f[1]:
            dr[f[0]] += 1
            continue
        mp[f[0]] += 1
        uq[f[0]] += f[4] == 'U'
        if f[5] == 'E' and f[4] == 'U' and f[2].startswith('HUMAN_'):
            wupg.write('\t'.join((f[0], f[1], f[2], f[3], f[6], f[7])) + '\n')
            nsl[f[0]] += 1
wupg.close()
with open(f'{SP}/perwell_mapped.tsv', 'w') as f:
    [f.write(f'{w}\t{mp[w]}\t{uq[w]}\n') for w in WELLS]
with open(f'{SP}/n_umi_dropped_per_well.tsv', 'w') as f:
    [f.write(f'{w}\t{dr[w]}\n') for w in WELLS]
with open(f'{ST}/basic/reads_per_well.tsv', 'w') as f:
    [f.write(f'{w}\t{nsl[w]}\n') for w in sorted(WELLS)]
print(SET, 'mapped/unique/usable per well:', {WELLS[w]: (mp[w], uq[w], nsl[w]) for w in WELLS}, flush=True)
r = subprocess.run(['python3', f'{CODE}/rarefy.py', ARM, 'E_U_supp21'], capture_output=True, text=True)
print(r.stdout.strip()[-200:], r.stderr.strip()[-300:], flush=True)
Um = collections.defaultdict(list)
Gm = collections.defaultdict(list)
for l in open(f'{ST}/wupg.tsv'):
    f = l.rstrip('\n').split('\t')
    Um[f[0]].append(f[1])
    Gm[f[0]].append(f[4])
with open(f'{ST}/basic/perwell_B.tsv', 'w') as out:
    out.write(
        'well\tB_raw\tB_kept\tumis_excluded\tB_corr\tgenes_B\tgenes_saturated\ttags_seen\ttags_homopolymer\ttags_overused\n'
    )
    for w in sorted(Um):
        b = b_stats(Um[w], Gm[w])
        out.write(
            f"{w}\t{b['B_raw']}\t{b['B_kept']}\t{b['excluded']}\t{b['B_corr']}\t{b['genes']}\t{b['saturated']}\t{b['tags']}\t{b['hp']}\t{b['ov']}\n"
        )
        print(WELLS[w], 'usable', len(Um[w]), 'UCI', b['B_corr'], 'genes', b['genes'], flush=True)
if SET == '50nt':
    PB = f'{U}/per_sample_bams'
    os.makedirs(f'{PB}/bobseq_50nt_supp21', exist_ok=True)
    src = pysam.AlignmentFile(f'{SP}/Aligned.out.bam', 'rb')
    hd = src.header.to_dict()
    keep = [(i, sq) for i, sq in enumerate(hd['SQ']) if sq['SN'].startswith('HUMAN_')]
    newh = {'HD': {'VN': '1.6', 'SO': 'unsorted'}, 'SQ': [{'SN': sq['SN'][6:], 'LN': sq['LN']} for _, sq in keep]}
    tid_map = {i: k for k, (i, _) in enumerate(keep)}
    outs = {
        w: pysam.AlignmentFile(f'{PB}/bobseq_50nt_supp21/{n}.unsorted.bam', 'wb', header=newh) for w, n in WELLS.items()
    }
    cnt = collections.Counter()
    for a in src.fetch(until_eof=True):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.reference_id not in tid_map:
            continue
        w = a.query_name.rsplit('_', 2)[1]
        if w not in outs:
            continue
        b = pysam.AlignedSegment(outs[w].header)
        b.query_name = a.query_name
        b.flag = a.flag
        b.reference_id = tid_map[a.reference_id]
        b.reference_start = a.reference_start
        b.mapping_quality = a.mapping_quality
        b.cigar = a.cigar
        b.next_reference_id = -1
        b.next_reference_start = -1
        b.template_length = 0
        b.query_sequence = a.query_sequence
        b.query_qualities = a.query_qualities
        b.set_tags(a.get_tags())
        outs[w].write(b)
        cnt[w] += 1
    for w, o in outs.items():
        o.close()
        n = WELLS[w]
        pysam.sort('-o', f'{PB}/bobseq_50nt_supp21/{n}.bam', f'{PB}/bobseq_50nt_supp21/{n}.unsorted.bam')
        pysam.index(f'{PB}/bobseq_50nt_supp21/{n}.bam')
        os.remove(f'{PB}/bobseq_50nt_supp21/{n}.unsorted.bam')
        print(n, '50-nt per-well BAM', cnt[w], 'primary human records', flush=True)
    for k, n in enumerate(WELLS.values(), 1):
        r = subprocess.run(
            [
                'python3',
                f'{CODE}/build_positions.py',
                f'bobseq-ris-500mM-{k}',
                f'{PB}/bobseq_50nt_supp21/{n}.bam',
                f'{COVERAGE}/positions_canonical_50nt/bobseq-ris-500mM-{k}.npz',
            ],
            capture_output=True,
            text=True,
        )
        print(r.stdout.strip()[-160:], r.stderr.strip()[-200:], flush=True)
    os.environ['BM_COVSET'] = '50nt'
    import composition_rules as cr

    res = {}
    for n in WELLS.values():
        m, name, c = cr.classify_bam(('BOBseq', n, f'bobseq_50nt_supp21/{n}.bam', False))
        res[name] = c
        print(name, c, flush=True)
    json.dump(res, open(f'{PR}/results/supplement_21/composition_rules_ris500_50nt.json', 'w'), indent=1)
    MET = f'{P}/preprint_figures_50nt/coverage_architecture/picard_metrics/matched'
    REFFLAT = f'{COVERAGE}/picard_refflat_pc_noRP_noMT_bare.txt'
    for n in WELLS.values():
        out = f'{MET}/{n}.RNA_Metrics.txt'
        if os.path.exists(out) and '## HISTOGRAM' in open(out).read():
            print(n, 'picard cached')
            continue
        cmd = [
            'picard',
            'CollectRnaSeqMetrics',
            '-I',
            f'{PB}/bobseq_50nt_supp21/{n}.bam',
            '-O',
            f'{MET}/{n}.RNA_Metrics.txt',
            '--REF_FLAT',
            REFFLAT,
            '--STRAND_SPECIFICITY',
            'NONE',
            '--MINIMUM_LENGTH',
            '1000',
            '--VALIDATION_STRINGENCY',
            'LENIENT',
        ]
        r = subprocess.run(cmd, capture_output=True, text=True)
        open(f'{MET}/{n}.log', 'w').write(r.stdout + r.stderr)
        print(n, 'picard', 'done' if r.returncode == 0 else 'FAILED', flush=True)
print('SUPP21 TABLES DONE', SET)

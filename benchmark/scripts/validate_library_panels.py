#!/usr/bin/env python3
"""Independent checks of the library panels: p0cd / p0dd (reads before and after dedup, duplicate
rate), p7a (read length), p7b (insert length). Every number is recomputed from the raw inputs without the generator
code, and the provenance (which run, which wells, which files) is asserted. Output: <native panel dir>/VALIDATION_library_panels.md
"""

from settings import WORK, PROJECT, COVERAGE
import os, sys, csv, json, glob, subprocess, numpy as np, pysam

H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, H)
from palette import METH
from umi_b_lib import b_stats
import coverage_lib as cl

P = WORK
U = f'{P}/benchmark_uniform'
PR = PROJECT
D = f'{P}/preprint_figures_native/main_figure_panels'
V = f'{D}/values'
ARMS = {'DRUG-seq': 'drugseq_native', 'prime-seq': 'primeseq_native', 'BOBseq': 'bob57_24plex_pe_native'}
BAMD = {'DRUG-seq': 'drugseq_native', 'prime-seq': 'primeseq_native', 'BOBseq': 'bobseq_pe_native'}
EXP_N = {'DRUG-seq': 24, 'prime-seq': 8, 'BOBseq': 18}
L = []
fails = 0


def check(name, ok, detail=''):
    global fails
    fails += not ok
    L.append(f"{'PASS' if ok else 'FAIL'} {name}{(': ' + detail) if detail else ''}")
    print(L[-1], flush=True)


# ---------- provenance ----------
readme = open(f'{U}/{ARMS["BOBseq"]}/README.native.txt').read()
check('BOBseq arm is the paired-end run (README.native.txt)', 'paired-end 2x150' in readme, readme.strip()[:90])
hdr = pysam.AlignmentFile(f'{U}/{ARMS["BOBseq"]}/Aligned.out.bam', 'rb').header.to_dict()
cl_ = next((x.get('CL', '') for x in hdr.get('PG', []) if x.get('ID') == 'STAR'), '')
check(
    'BOBseq arm BAM: STAR ran with two read files (paired-end)',
    cl_.count('.fastq.gz') >= 2
    or '--readFilesIn' in cl_
    and len(cl_.split('--readFilesIn')[1].split('--')[0].split()) == 2,
    cl_.split('--readFilesIn')[1].split('--')[0].strip()[:120] if '--readFilesIn' in cl_ else 'no CL',
)
exc = {
    l.split('\t')[0]
    for l in open(f'{U}/{ARMS["BOBseq"]}/units_excluded.tsv')
    if l.strip() and not l.startswith('#') and not l.startswith('unit')
}
keep = {l.strip() for l in open(f'{U}/{ARMS["BOBseq"]}/wells_keep.txt') if l.strip()}
wells = {
    m: {l.split('\t')[0] for l in open(f'{U}/{ARMS[m]}/stats_E_U/basic/reads_per_well.tsv') if l.strip()} for m in METH
}
check(
    'BOBseq wells in the read table = wells_keep minus units_excluded (18)',
    wells['BOBseq'] == (keep - exc) and len(wells['BOBseq']) == 18,
    f'{len(wells["BOBseq"])} wells, excluded {sorted(exc)}',
)
for m in METH:
    check(f'{m}: number of wells in the per-well tables', len(wells[m]) == EXP_N[m], f'{len(wells[m])}')
# ---------- p0cd / p0dd ----------
tsv = list(csv.DictReader(open(f'{V}/p0c_p0d_dedup_per_sample_native.tsv'), delimiter='\t'))
pj = json.load(open(f'{D}/p0c_reads_after_dedup_native.json'))
dj = json.load(
    open(
        next(
            x
            for x in (f'{D}/p0dd_duplicate_rate_native.json', f'{D}/extra/p0dd_duplicate_rate_native.json')
            if os.path.exists(x)
        )
    )
)
d2j = json.load(open(f'{D}/p0d2_duplicate_rate_at_250k_native.json'))
for m in METH:
    rows = [r for r in tsv if r['method'] == m]
    check(
        f'{m}: values TSV has one row per benchmark well', {r['sample'] for r in rows} == wells[m], f'{len(rows)} rows'
    )
    rar = {
        (r['well'], int(r['depth'])): r
        for r in csv.DictReader(open(f'{U}/{ARMS[m]}/stats_E_U/rarefied.tsv'), delimiter='\t')
    }
    nat = {w: float(next(r for (ww, d), r in rar.items() if ww == w)['reads_native']) for w in wells[m]}
    check(
        f'{m}: filtered reads per well == rarefied.tsv reads_native (the p1 / p0b source)',
        all(abs(float(r['filtered_reads']) - nat[r['sample']]) < 0.5 for r in rows),
    )
    B = {r['well']: r for r in csv.DictReader(open(f'{U}/{ARMS[m]}/stats_E_U/basic/perwell_B.tsv'), delimiter='\t')}
    check(
        f'{m}: molecules per well == perwell_B.tsv B_corr (the p1 currency)',
        all(abs(float(r['molecules_UCI']) - float(B[r['sample']]['B_corr'])) < 0.5 for r in rows),
    )
    check(
        f'{m}: duplicate rate == 100 x (1 - molecules / filtered reads)',
        all(
            abs(float(r['duplicate_rate_pct']) - 100 * (1 - float(r['molecules_UCI']) / float(r['filtered_reads'])))
            < 0.01
            for r in rows
        ),
    )
    ok = True
    n250 = 0
    for r in rows:
        k = (r['sample'], 250000)
        if k in rar and rar[k]['at_native'] == '0':
            n250 += 1
            ok &= abs(float(r['duplicate_rate_pct_at_250k']) - 100 * (1 - float(rar[k]['B_corr']) / 250000)) < 0.01
        else:
            ok &= r['duplicate_rate_pct_at_250k'] == ''
    check(
        f'{m}: duplicate rate at 250k == 1 - rarefied B_corr@250k / 250k, wells below 250k blank',
        ok,
        f'{n250} wells at 250k',
    )
    check(
        f'{m}: p0c json median molecules_UCI == median of the values TSV',
        abs(pj[m]['median'] - np.median([float(r['molecules_UCI']) for r in rows])) < 0.5,
        f"{pj[m]['median']:,.0f}",
    )
    check(
        f'{m}: p0d2 (figure panel) json median == median of the values TSV at 250k',
        abs(
            d2j[m]['median']
            - np.median([float(r['duplicate_rate_pct_at_250k']) for r in rows if r['duplicate_rate_pct_at_250k'] != ''])
        )
        < 0.01
        and d2j[m]['n'] == sum(1 for r in rows if r['duplicate_rate_pct_at_250k'] != ''),
        f"{d2j[m]['median']:.1f}% (n={d2j[m]['n']})",
    )
    for key in ('duplicate_rate_pct', 'duplicate_rate_pct_at_250k'):
        check(
            f'{m}: p0dd json median {key} == median of the values TSV',
            abs(dj[m][key]['median'] - np.median([float(r[key]) for r in rows if r[key] != ''])) < 0.01,
            f"{dj[m][key]['median']:.1f}%",
        )
    # independent recomputation of filtered reads and B_corr from the read-table slice for two wells
    pick = sorted(wells[m])[:2]
    cnt = {w: 0 for w in pick}
    um = {w: [] for w in pick}
    ge = {w: [] for w in pick}
    for line in open(f'{U}/{ARMS[m]}/stats_E_U/wupg.tsv'):
        w = line[: line.index('\t')]
        if w in cnt:
            f = line.rstrip('\n').split('\t')
            cnt[w] += 1
            um[w].append(f[1])
            ge[w].append(f[4])
    for w in pick:
        bs = b_stats(um[w], ge[w])
        r = next(r for r in rows if r['sample'] == w)
        check(
            f'{m} well {w}: filtered reads recounted from wupg.tsv == panel value',
            cnt[w] == int(float(r['filtered_reads'])),
            f'{cnt[w]:,}',
        )
        check(
            f'{m} well {w}: UCI recomputed with umi_b_lib.b_stats == panel value',
            bs['B_corr'] == int(float(r['molecules_UCI'])),
            f"{bs['B_corr']:,} vs {float(r['molecules_UCI']):,.0f}",
        )
# ---------- p7a read length ----------
R = json.load(open(f'{COVERAGE}/read_length_native.json'))
S = [l.rstrip('\n').split('\t') for l in open(f'{COVERAGE}/samples.tsv') if l.strip()]
check(
    'p7a: positions sample table = the benchmark sample set',
    {m: sum(1 for s in S if s[2] == m) for m in METH} == EXP_N,
    str({m: sum(1 for s in S if s[2] == m) for m in METH}),
)
for m in METH:
    lab, bam, _ = next(s for s in S if s[2] == m)
    z = np.load(f'{COVERAGE}/positions_canonical/{lab}.npz')
    al = z['alen'].astype(int)
    mate = z['mate']
    f = f'{U}/per_sample_bams/{bam}'
    AL = {1: [], 2: [], 0: []}
    n = 0
    for a in pysam.AlignmentFile(f, 'rb').fetch(until_eof=True):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality != 255:
            continue
        AL[1 if a.is_read1 else (2 if a.is_read2 else 0)].append(sum(ln for op, ln in a.cigartuples if op in (0, 7, 8)))
        n += 1
        if n >= 400000:
            break
    if m == 'BOBseq':
        check(
            'p7a BOBseq: per-sample BAM holds both mates',
            len(AL[1]) > 0 and len(AL[2]) > 0,
            f'mate1 {len(AL[1]):,}, mate2 {len(AL[2]):,} in the first 400k unique records',
        )
        for k in (1, 2):
            check(
                f'p7a BOBseq mate {k} ({lab}): aligned-length median from the BAM == npz median',
                abs(np.median(AL[k]) - np.median(al[mate == k])) <= 2,
                f'BAM {np.median(AL[k]):.0f} vs npz {np.median(al[mate == k]):.0f}',
            )
    else:
        v = AL[0] + AL[1] + AL[2]
        check(
            f'p7a {m} ({lab}): aligned-length median from the BAM == npz median',
            abs(np.median(v) - np.median(al)) <= 2,
            f'BAM {np.median(v):.0f} vs npz {np.median(al):.0f}',
        )
    tot = int(R[m]['n_reads']) if m in R else None
    check(
        f'p7a {m}: json read count == sum over the set\'s npz',
        tot == sum(len(np.load(f"{COVERAGE}/positions_canonical/{s[0]}.npz")['alen']) for s in S if s[2] == m),
        f'{tot:,}',
    )
# ---------- p7e / p7f aligned bases ----------
AB = list(csv.DictReader(open(f'{V}/p7e_p7f_aligned_bases_per_sample_native.tsv'), delimiter='\t'))
ej = json.load(open(f'{D}/p7e_aligned_bases_total_native.json'))
fj = json.load(
    open(
        next(
            x
            for x in (
                f'{D}/p7f_aligned_bases_per_fragment_native.json',
                f'{D}/extra/p7f_aligned_bases_per_fragment_native.json',
            )
            if os.path.exists(x)
        )
    )
)
check('p7e: one row per benchmark sample', {m: sum(1 for r in AB if r['method'] == m) for m in METH} == EXP_N)
ok_t = ok_f = True
for r in AB:
    z = np.load(f"{COVERAGE}/positions_canonical/{r['sample']}.npz")
    al = z['alen'].astype(np.int64)
    mate = z['mate']
    nf = int((mate == 1).sum()) if (mate > 0).any() else len(al)
    ok_t &= int(r['total_aligned_bases']) == int(al.sum()) and int(r['reads']) == len(al)
    ok_f &= abs(float(r['bases_per_fragment']) - al.sum() / nf) < 0.01 and int(r['fragments']) == nf
check('p7e: total aligned bases per sample == sum of aligned lengths in the position tables (all 50 samples)', ok_t)
check('p7f: bases per fragment == total / fragments (BOBseq fragments = mate-1 records) for all 50 samples', ok_f)
for m in METH:
    check(
        f'p7e {m}: json median (Gb) == median of the values TSV',
        abs(
            ej[m]['total_aligned_bases']['median']
            - np.median([int(r['total_aligned_bases']) for r in AB if r['method'] == m]) / 1e9
        )
        < 1e-6,
        f"{ej[m]['total_aligned_bases']['median']:.3f} Gb",
    )
    check(
        f'p7e {m}: matched-depth lane == 250k x bases per fragment, json median == values TSV',
        all(
            abs(float(r['aligned_bases_at_250k']) - 250000 * float(r['bases_per_fragment'])) < 1
            for r in AB
            if r['method'] == m
        )
        and abs(
            ej[m]['aligned_bases_at_250k']['median']
            - np.median([float(r['aligned_bases_at_250k']) for r in AB if r['method'] == m]) / 1e9
        )
        < 1e-6,
        f"{ej[m]['aligned_bases_at_250k']['median']:.4f} Gb",
    )
    check(
        f'p7f {m}: json median == median of the values TSV',
        abs(fj[m]['median'] - np.median([float(r['bases_per_fragment']) for r in AB if r['method'] == m])) < 0.01,
        f"{fj[m]['median']:.1f} nt",
    )
lab, bam, _ = next(s_ for s_ in S if s_[2] == 'BOBseq')
tot = 0
n1 = 0
n = 0
for a in pysam.AlignmentFile(f'{U}/per_sample_bams/{bam}', 'rb').fetch(until_eof=True):
    if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality != 255:
        continue
    n += 1
    tot += sum(ln for op, ln in a.cigartuples if op in (0, 7, 8))
    n1 += a.is_read1
    if n >= 400000:
        break
r = next(r for r in AB if r['sample'] == lab)
z = np.load(f"{COVERAGE}/positions_canonical/{lab}.npz")
al = z['alen'].astype(np.int64)
check(
    f'p7f BOBseq ({lab}): bases per fragment from the BAM (first 400k unique reads, both mates) within 3% of the panel value',
    abs((tot / max(n1, 1)) - float(r['bases_per_fragment'])) / float(r['bases_per_fragment']) < 0.03,
    f'BAM {tot / max(n1, 1):.0f} vs panel {float(r["bases_per_fragment"]):.0f} nt',
)
# ---------- p7b insert length ----------
ib = f'{PR}/data/bams/pipeline_dedup_pe/dedup_pairs_6NHH'
inv = {
    l.split('\t')[0].split('/')[-1]: int(l.split('\t')[1])
    for l in open(f'{PR}/data/gcs_inventory.tsv')
    if 'dedup_pairs_6NHH/' in l
}
check(
    'p7b: 24 staged pair BAMs have the sizes of the backup (pipeline_dedup_pe/dedup_pairs_6NHH)',
    all(os.path.getsize(f'{ib}/{k}') == v for k, v in inv.items() if k.endswith('.bam')),
    f'{sum(1 for k in inv if k.endswith(".bam"))} BAMs',
)
jw = {os.path.basename(f)[:-5] for f in glob.glob(f'{PR}/results/library_metrics/insert_per_well/*.json')}
pb = (
    json.load(open(f'{D}/p7b_insert_length_bobseq.json'))
    if os.path.exists(f'{D}/p7b_insert_length_bobseq.json')
    else None
)
if pb:
    used = set(pb['wells'])
    check(
        'p7b: wells used = 18 benchmark wells (no Ris 500 nM, no mouse)',
        len(used) == 18 and not any(w.startswith('Ris_dose_500mM') or w.startswith('RAW') for w in used),
        f'{len(used)} of {len(jw)} per-well results',
    )
    tx, idx = cl.build()
    Hh = {}
    same = 0
    tot = 0
    diff = []
    for a in pysam.AlignmentFile(f'{ib}/Hek_control_1_star_pairs.bam', 'rb').fetch(until_eof=True):
        if (
            a.is_unmapped
            or a.is_secondary
            or a.is_supplementary
            or a.mapping_quality != 255
            or not a.reference_name.startswith('HUMAN_')
        ):
            continue
        spliced = any(op == 3 for op, _ in a.cigartuples)
        p5 = a.reference_end - 1 if a.is_reverse else a.reference_start
        r = cl.position(idx, a.reference_name[6:], p5)
        if r is None:
            continue
        q = a.query_name
        m = Hh.pop(q, None)
        if m is None:
            Hh[q] = (r[0], r[1], abs(a.template_length), spliced)
            continue
        if m[0] != r[0]:
            continue
        tot += 1
        if not spliced and not m[3]:
            ins = abs(m[1] - r[1]) + 1
            same += ins == m[2]
            diff.append(ins - m[2])
        if tot >= 200000:
            break
    check(
        'p7b: for unspliced pairs the transcript-space insert equals the genomic TLEN (projection check, Hek_control_1)',
        len(diff) > 1000 and same / len(diff) >= 0.95,
        f'{same:,} of {len(diff):,} equal ({100*same/max(len(diff),1):.1f}%), median |diff| {np.median(np.abs(diff)) if diff else "na"}',
    )
    h = np.array(
        json.load(open(f'{PR}/results/library_metrics/insert_per_well/Hek_control_1.json'))['hist_1nt_to_3000'], float
    )
    cum = np.cumsum(h) / h.sum()
    check(
        'p7b: per-well json median (Hek_control_1) == json insert_median',
        abs(int(np.searchsorted(cum, 0.5)) - pb['wells']['Hek_control_1']['insert_median']) <= 1,
        f"{pb['wells']['Hek_control_1']['insert_median']:.0f}",
    )
else:
    check('p7b: panel json present', False, 'not built yet')
open(f'{D}/VALIDATION_library_panels.md', 'w').write(
    '# Validation of the library panels (p0cd, p0dd, p7a, p7b), native set, '
    + __import__('datetime').date.today().isoformat()
    + '\n\n'
    + '\n'.join('- ' + x for x in L)
    + f'\n\n{fails} FAILED of {len(L)} checks\n'
)
print(f'{fails} FAILED of {len(L)} checks')

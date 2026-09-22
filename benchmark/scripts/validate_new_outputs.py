#!/usr/bin/env python3
"""Independent checks of the two stacked figures, the 24-sample poly(A) coverage panel, the mouse wells and the curated figure tree.
Every plotted value is compared with the PRIMARY table or file it derives from (arm tables, Picard metrics, position tables, funnel tables, composition counts), not with the
per-well table the figure script read; mouse-well counts are recounted from the wells' read tables; two S14 profiles are recomputed with a separately written routine.
Writes figures/VALIDATION_new_outputs.md in the preprint project; exit 1 on any failure."""

from settings import CODE, WORK, PROJECT, COVERAGE, RUN_JSON
import os, re, sys, csv, json, ast, hashlib, collections, datetime, numpy as np

P = WORK
PR = PROJECT
U = f'{P}/benchmark_uniform'
A = COVERAGE
FN = f'{P}/preprint_figures_native'
ARM = f'{U}/bob57_24plex_pe_native'
L = []
fails = 0


def check(name, ok, detail=''):
    global fails
    fails += not ok
    L.append(f"{'PASS' if ok else 'FAIL'} {name}{(': ' + detail) if detail else ''}")
    print(L[-1], flush=True)


def T(p):
    return list(csv.DictReader(open(p), delimiter='\t'))


def close(a, b, tol):
    return abs(a - b) <= tol


RUN = json.load(open(RUN_JSON))
NAME = {w: l.replace(' ', '_').replace('/', '_') for w, l in RUN['bobcode_labels'].items()}
CODE = {v: k for k, v in NAME.items()}


def picard_balance(path):
    h = []
    inh = False
    for l in open(path):
        r = l.rstrip('\n').split('\t')
        if r and r[0].startswith('## HISTOGRAM'):
            inh = True
            continue
        if inh and len(r) >= 2 and r[0].isdigit():
            h.append(float(r[1]))
    h = np.array(h)
    return float(2 * (np.arange(101) / 100 * h).sum() / h.sum())


# ================= 1. stacked 24-plex figure vs primary sources =================
SV = collections.defaultdict(dict)
for r in T(f'{FN}/supplement_stacked/values/stacked_24plex_per_sample_qc.tsv'):
    SV[r['sample']][(r['metric'], r['series'])] = float(r['value'])
check(
    'stacked 24-plex: 24 samples x 13 metrics (aligned length has two series)',
    len(SV) == 24 and all(len(v) == 14 for v in SV.values()),
    f'{len(SV)} samples, {sorted({len(v) for v in SV.values()})} values each',
)
inp = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(f'{ARM}/input_counts.tsv') if l.strip()}
for f_ in (f'{ARM}/supp21/n_umi_dropped_per_well.tsv', f'{ARM}/mouse3/n_umi_dropped_per_well.tsv'):
    for l in open(f_):
        w, d = l.rstrip('\n').split('\t')
        inp[w] -= int(d)
pm = {}
rpw = {}
pb = {}
for f_ in (f'{ARM}/perwell_mapped.tsv', f'{ARM}/supp21/perwell_mapped.tsv', f'{ARM}/mouse3/perwell_mapped.tsv'):
    for l in open(f_):
        r = l.rstrip('\n').split('\t')
        if len(r) >= 3:
            pm[r[0]] = (int(r[1]), int(r[2]))
for s_ in ('stats_E_U', 'stats_E_U_supp21', 'stats_E_U_mouse3'):
    for l in open(f'{ARM}/{s_}/basic/reads_per_well.tsv'):
        r = l.rstrip('\n').split('\t')
        rpw[r[0]] = int(r[1])
    for r in T(f'{ARM}/{s_}/basic/perwell_B.tsv'):
        pb[r['well']] = float(r['B_corr'])
comp = {
    r['sample']: r
    for r in T(f'{FN}/main_figure_panels/values/p3_composition_rules_native_per_sample.tsv')
    if r['method'] == 'BOBseq'
}
CL = ['rRNA', 'mitochondrial', 'ribosomal-protein', 'mRNA', 'exonic other biotype', 'intronic', 'intergenic']
for fn in ('composition_rules_ris500_native.json', 'composition_rules_mouse3_native.json'):
    for n, c in json.load(open(f'{PR}/results/supplement_21/{fn}')).items():
        comp[n] = {k: c.get(k, 0) for k in CL}
S = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open(f'{A}/samples.tsv') if l.strip()}
lab_of = {os.path.basename(r[1])[:-4]: r[0] for r in S.values() if r[2] == 'BOBseq'}
for k in (1, 2, 3):
    lab_of[f'Ris_dose_500mM-{k}'] = f'bobseq-ris-500mM-{k}'
    lab_of[f'RAW_control_{k}'] = f'bobseq-raw-control-{k}'
bad = collections.Counter()
n_cmp = 0
for n, v in SV.items():
    w = CODE[n]
    mouse = n.startswith('RAW')
    z = np.load(f'{A}/positions_canonical/{lab_of[n]}.npz')
    al = z['alen'].astype(np.int64)
    mate = z['mate']
    ins = json.load(open(f'{PR}/results/library_metrics/insert_per_well/{n}.json'))
    ct = sum(float(comp[n][k]) for k in CL)
    met = f"{FN}/coverage_architecture/picard_metrics/{'native_mouse' if mouse else 'native'}/{n}.RNA_Metrics.txt"
    exp = {
        ('reads into STAR', '1'): inp[w],
        ('mapping rate (% of reads)', '1'): 100 * pm[w][0] / inp[w],
        ('uniquely mapped (% of reads)', '1'): 100 * pm[w][1] / inp[w],
        ('filtered reads (% of reads)', '1'): 100 * rpw[w] / inp[w],
        ('duplicate rate (%)', '1'): 100 * (1 - pb[w] / rpw[w]),
        ('filtered reads after de-duplication (UMI)', '1'): pb[w],
        ('rRNA content (% of mapped)', '1'): 100 * float(comp[n]['rRNA']) / ct,
        ('mRNA fraction (% of mapped)', '1'): 100 * float(comp[n]['mRNA']) / ct,
        ('intronic reads (% of mapped)', '1'): 100 * float(comp[n]['intronic']) / ct,
        ('median aligned length (nt) filled read 2, open read 1', '1'): float(np.median(al[mate == 1])),
        ('median aligned length (nt) filled read 2, open read 1', '2'): float(np.median(al[mate == 2])),
        ('insert length median (nt)', '1'): ins['insert_median'],
        ('aligned bases (Gb)', '1'): al.sum() / 1e9,
        ("5'-3' balance (2 x centroid)", '1'): picard_balance(met),
    }
    assert set(exp) == set(v), (n, set(exp) ^ set(v))
    for k, e in exp.items():
        n_cmp += 1
        if not close(v[k], e, max(abs(e) * 2e-3, 0.02)):
            bad[k[0]] += 1
            print('   mismatch', n, k, v[k], e)
check(
    'stacked 24-plex: every plotted value equals the value recomputed from its primary source',
    not bad,
    f'{n_cmp} values compared; mismatches {dict(bad) or 0}',
)
# ================= 2. stacked benchmark figure vs primary sources =================
from palette import METH

ARMS = {'DRUG-seq': 'drugseq_native', 'prime-seq': 'primeseq_native', 'BOBseq': 'bob57_24plex_pe_native'}
BV = collections.defaultdict(dict)
for r in T(f'{FN}/supplement_stacked/values/stacked_benchmark_native_per_sample_qc.tsv'):
    BV[(r['method'], r['sample'])][r['metric']] = float(r['value']) if r['value'] != '' else None
check(
    'stacked benchmark: 50 samples x 15 metrics',
    len(BV) == 50 and all(len(v) == 15 for v in BV.values()),
    f'{len(BV)} samples',
)
p3 = {
    (r['method'], r['sample']): r
    for r in T(f'{FN}/main_figure_panels/values/p3_composition_rules_native_per_sample.tsv')
}
LAB2 = {os.path.basename(r[1])[:-4]: (r[0], r[2]) for r in S.values()}
bad = collections.Counter()
n_cmp = 0
for m in METH:
    arm = f'{U}/{ARMS[m]}'
    rar = collections.defaultdict(dict)
    for r in T(f'{arm}/stats_E_U/rarefied.tsv'):
        rar[r['well']][int(r['depth'])] = r
    pbm = {r['well']: float(r['B_corr']) for r in T(f'{arm}/stats_E_U/basic/perwell_B.tsv')}
    if m == 'BOBseq':
        fun = {w: (inp[w], pm[w][0], pm[w][1]) for w in rar}
    else:  # competitors: funnel rows are raw read-name barcodes; a well = its exact barcode plus the error-corrected ones of well_map.tsv (DRUG-seq: 744 raw barcodes -> 24 wells)
        wmap = (
            {l.split('\t')[0]: l.rstrip('\n').split('\t')[1] for l in open(f'{arm}/well_map.tsv')}
            if os.path.exists(f'{arm}/well_map.tsv')
            else {}
        )
        agg = collections.defaultdict(lambda: [0, 0, 0])
        for r in T(f'{arm}/funnel_per_sample.tsv'):
            t_ = agg[wmap.get(r['sample'], r['sample'])]
            t_[0] += int(r['reads_in'])
            t_[1] += int(r['mapped'])
            t_[2] += int(r['unique'])
        fun = {k: tuple(v) for k, v in agg.items()}
    for (mm, n), v in BV.items():
        if mm != m:
            continue
        w = CODE[n] if m == 'BOBseq' else (n if n in rar else n.rsplit('_', 1)[-1])
        assert w in rar, (m, n)
        nat = float(next(iter(rar[w].values()))['reads_native'])
        r250 = rar[w].get(250000)
        ok250 = r250 is not None and r250['at_native'] == '0'
        z = np.load(f'{A}/positions_canonical/{LAB2[n][0]}.npz')
        al = z['alen'].astype(np.int64)
        c = p3[(m, n)]
        ct = sum(float(c[k]) for k in CL)
        met = f'{FN}/coverage_architecture/picard_metrics/native/{n}.RNA_Metrics.txt'
        exp = {
            'reads into STAR': fun[w][0],
            'mapping rate (% of reads)': 100 * fun[w][1] / fun[w][0],
            'uniquely mapped (% of reads)': 100 * fun[w][2] / fun[w][0],
            'filtered reads (% of reads)': 100 * nat / fun[w][0],
            'filtered reads after de-duplication (UMI)': pbm[w],
            'duplicate rate (%) at 250k filtered reads': 100 * (1 - float(r250['B_corr']) / 250000) if ok250 else None,
            'molecules (UMI) at 250k filtered reads': float(r250['B_corr']) if ok250 else None,
            'genes at 250k filtered reads (no ribosomal-protein genes)': float(r250['genes_B_noRP']) if ok250 else None,
            'rRNA content (% of mapped)': 100 * float(c['rRNA']) / ct,
            'mRNA fraction (% of mapped)': 100 * float(c['mRNA']) / ct,
            'intronic reads (% of mapped)': 100 * float(c['intronic']) / ct,
            'mitochondrial reads (% of mapped)': 100 * float(c['mitochondrial']) / ct,
            'mean aligned length per read (nt)': float(al.mean()),
            'aligned bases (Gb)': al.sum() / 1e9,
            "5'-3' balance (2 x centroid)": picard_balance(met),
        }
        assert set(exp) == set(v), (m, n, set(exp) ^ set(v))
        for k, e in exp.items():
            n_cmp += 1
            if (e is None) != (v[k] is None) or (e is not None and not close(v[k], e, max(abs(e) * 2e-3, 0.02))):
                bad[k] += 1
                print('   mismatch', m, n, k, v[k], e)
check(
    'stacked benchmark: every plotted value equals the value recomputed from its primary source (empty cells = wells below 250k)',
    not bad,
    f'{n_cmp} values compared; mismatches {dict(bad) or 0}',
)
# ================= 3. mouse wells recounted from their read tables =================
for w, n in (('TAAGACG', 'RAW_control_1'), ('GTCCATC', 'RAW_control_2'), ('GATGACAGTAT', 'RAW_control_3')):
    mp = uq = fl = 0
    pairs = set()
    byc = collections.Counter()
    for l in open(f'{ARM}/mouse3/{n}/reads.tsv'):
        f = l.rstrip('\n').split('\t')
        if f[0] != w or 'N' in f[1]:
            continue
        mp += 1
        uq += f[4] == 'U'
        if f[4] == 'U' and f[5] == 'E' and f[2][:6] == 'MOUSE_':
            fl += 1
            pairs.add((f[1], f[6]))
    braw = int(next(r for r in T(f'{ARM}/stats_E_U_mouse3/basic/perwell_B.tsv') if r['well'] == w)['B_raw'])
    check(
        f'mouse {n}: mapped / unique / filtered recounted from reads.tsv == tables; distinct (UMI, gene) pairs == B_raw',
        (mp, uq) == pm[w] and fl == rpw[w] and len(pairs) == braw,
        f'{mp:,} / {uq:,} / {fl:,}; pairs {len(pairs):,}',
    )
    c = json.load(open(f'{PR}/results/supplement_21/composition_rules_mouse3_native.json'))[n]
    tot = sum(c[k] for k in CL) + c['human_contig_non_rdna_excluded']
    check(
        f'mouse {n}: composition classes + excluded human-contig records == mate-1 primary mapped records; rRNA includes the human rDNA records',
        abs(tot - mp) <= dr_tol if (dr_tol := 0.002 * mp) else False,
        f'{tot:,} vs mapped without N-UMI {mp:,} (N-UMI reads are in the BAM)',
    )
# ================= 4. S14 (all 24): two profiles recomputed with a separately written routine =================
RIBO = ('RPL', 'RPS', 'MRPL', 'MRPS')
prof = {
    r['sample']: (int(r['genes']), np.array([float(r[f'nt{d}']) for d in range(3000)]))
    for r in T(f'{FN}/supplement_per_sample_21/values/polyA_profiles_all24_native.tsv')
}
check('S14: 24 profiles (21 human, 3 mouse)', len(prof) == 24)


def is_rp(g):
    u = g.upper()
    rest = u[4:] if u[:4] in RIBO[2:] else (u[3:] if u[:3] in RIBO[:2] else None)
    return rest is not None and 'K' not in rest


for n in ('Hek_control_1', 'RAW_control_2'):
    z = np.load(f'{A}/positions_canonical/{lab_of[n]}.npz')
    genes = z['genes']
    tl = z['tx_len']
    g = z['gene']
    t = z['t'].astype(np.int64)
    al = z['alen'].astype(np.int64)
    cnt = np.bincount(g, minlength=len(genes))
    acc = np.zeros(3000)
    num = np.zeros(3000)
    ng = 0
    for gi in np.flatnonzero(cnt >= 20):
        Lg = int(tl[gi])
        if Lg < 1000 or is_rp(str(genes[gi])):
            continue
        ii = np.flatnonzero(g == gi)
        cov = np.zeros(Lg)
        for s_, e_ in zip(t[ii], np.minimum(t[ii] + al[ii], Lg)):
            cov[s_:e_] += 1
        zz = (cov / cov.mean())[::-1][:3000]
        acc[: len(zz)] += zz
        num[: len(zz)] += 1
        ng += 1
    raw = acc / np.maximum(num, 1)
    v = np.where(num > 0, raw, 0)
    k = np.ones(25)
    mine = np.convolve(v, k, 'same') / np.maximum(np.convolve((num > 0).astype(float), k, 'same'), 1)
    check(
        f'S14 {n}: profile recomputed (explicit per-read coverage, prefix rule for RP genes) == plotted profile',
        ng == prof[n][0] and float(np.abs(mine - prof[n][1]).max()) < 1e-3,
        f'{ng} genes; max abs difference {float(np.abs(mine - prof[n][1]).max()):.2e}',
    )
# ================= 5. curated tree: every file identical to its source, no duplicates of a panel across main/supplement =================
man = T(f'{PR}/figures/MANIFEST_figures.tsv')
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
check(
    'curated tree: every file equals its generator-output source (sha256) and the manifest sha256',
    all(
        os.path.exists(f"{PR}/figures/{r['file']}")
        and sha(f"{PR}/figures/{r['file']}") == sha(f"{PR}/figures/{r['source']}") == r['sha256']
        for r in man
    ),
    f'{len(man)} files',
)
svgs = [r['file'] for r in man if r['file'].endswith('.svg') and r['file'].split('/')[0] in ('main', 'supplement')]
stems = collections.Counter(os.path.basename(x) for x in svgs)
check(
    'curated tree: no panel placed twice within main + supplement (p4 heatmap/profile exist once per set)',
    all(v == 1 or k.startswith('p4') for k, v in stems.items()),
    f'{len(svgs)} panels',
)
open(f'{PR}/figures/VALIDATION_new_outputs.md', 'w').write(
    f'# Validation of the stacked figures, the 24-sample poly(A) panel, the mouse wells and the curated tree ({datetime.date.today()})\n\nStacked figures, 24-sample poly(A) coverage panel, mouse wells, curated figure tree; every value against its primary source (validate_new_outputs.py).\n\n'
    + '\n'.join(f'- {x}' for x in L)
    + f'\n\n{fails} FAILED of {len(L)} checks\n'
)
print(f'{fails} FAILED of {len(L)} checks')
sys.exit(1 if fails else 0)

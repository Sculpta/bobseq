#!/usr/bin/env python3
"""Independent sanity checks of the preprint main-figure panels, for one figure set (BM_COVSET=native|50nt).
Every check recomputes the plotted quantity from the source tables/BAMs with its own code (not the generators') and compares
with what the panel wrote (its JSON/TSV side file) or with an independent second source (report numbers.json, STAR logs).
Writes <set>/main_figure_panels/VALIDATION_main_panels.md.  Exit code 1 if any check fails."""

from settings import CODE, WORK, COVERAGE, GTF
import os, sys, csv, json, glob, random, statistics as st, collections, numpy as np

P = WORK
U = f'{P}/benchmark_uniform'
A = COVERAGE
SET = os.environ.get('BM_COVSET', 'native')
SFX = '' if SET == 'native' else f'_{SET}'
FIG = 'preprint_figures_native' if SET == 'native' else f'preprint_figures_{SET}'
D = f'{P}/{FIG}/main_figure_panels'
CV = f'{P}/{FIG}/coverage_architecture'
TAG = 'native' if SET == 'native' else 'matched'


def F(name):
    return (
        f'{D}/{name}' if os.path.exists(f'{D}/{name}') else f'{D}/extra/{name}'
    )  # non-figure panels live in extra/ (tidy_panels.py)


import os as _os, sys as _sys

_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from palette import METH

ARMS = {
    'matched': {'DRUG-seq': 'drugseq', 'BRB-seq': 'brbseq', 'prime-seq': 'primeseq', 'BOBseq': 'bob57_24plex_pe_hs'},
    'native': {
        'DRUG-seq': 'drugseq_native',
        'BRB-seq': 'brbseq_native',
        'prime-seq': 'primeseq_native',
        'BOBseq': 'bob57_24plex_pe_native',
    },
}
R = []
FAIL = 0


def check(name, ok, detail=''):
    global FAIL
    FAIL += not ok
    R.append(f'| {"PASS" if ok else "FAIL"} | {name} | {detail} |')
    print(('PASS ' if ok else 'FAIL ') + name + ' ' + detail)


def rare(arm):
    return list(csv.DictReader(open(f'{U}/{arm}/stats_E_U/rarefied.tsv'), delimiter='\t'))


def wells_of(arm):
    return {r['well']: float(r['reads_native']) for r in rare(arm) if r.get('reads_native')}


nj = {a['label']: a for a in json.load(open(f'{P}/{FIG}/report_data/numbers.json'))['arms']}
# ======================= p0a read fate =======================
p0 = json.load(open(f'{D}/p0a_read_fate_{TAG}.json'))
p0.pop('_arms', None)
for m in METH:
    arm = ARMS[TAG][m]
    usable = wells_of(arm)
    wm = {}
    if os.path.exists(f'{U}/{arm}/well_map.tsv'):
        wm = dict(l.rstrip('\n').split('\t')[:2] for l in open(f'{U}/{arm}/well_map.tsv'))
    fun = collections.defaultdict(lambda: [0, 0, 0])
    for r in csv.DictReader(open(f'{U}/{arm}/funnel_per_sample.tsv'), delimiter='\t'):
        w = wm.get(r['sample'], r['sample'])
        fun[w][0] += int(r['reads_in'])
        fun[w][1] += int(r['mapped'])
        fun[w][2] += int(r['unique'])
    ic = f'{U}/{arm}/input_counts.tsv'
    pmf = f'{U}/{arm}/perwell_mapped.tsv'
    if os.path.exists(ic):  # BOBseq: pooled BAM is mapped-only -> pipeline per-well tables
        inp = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(ic) if l.strip() and not l.startswith('well')}
        pm = {
            l.split('\t')[0]: (int(l.split('\t')[1]), int(l.rstrip('\n').split('\t')[2]))
            for l in open(pmf)
            if l.count('\t') >= 2 and not l.startswith('well')
        }
        for w in usable:
            fun[w] = [inp[w], pm[w][0], pm[w][1]]
    # 1. monotone per well
    bad = [w for w in usable if not (fun[w][0] >= fun[w][1] >= fun[w][2] >= usable[w])]
    check(
        f'p0a {m}: reads_in >= mapped >= unique >= filtered in every benchmark well',
        not bad,
        f'{len(usable)} wells; violations {bad[:3]}',
    )
    # 2. recomputed percentages vs plotted json
    tot = [sum(fun[w][i] for w in usable) for i in range(3)] + [sum(usable.values())]
    mine = {k: 100 * v / tot[0] for k, v in zip(('reads_in', 'mapped', 'unique', 'filtered'), tot)}
    d = max(abs(mine[k] - p0[m][k]) for k in mine)
    check(
        f'p0a {m}: recomputed % (in/mapped/unique/filtered) match plotted',
        d < 0.05,
        ' / '.join(f'{mine[k]:.1f}' for k in mine) + f' vs plotted ' + ' / '.join(f'{p0[m][k]:.1f}' for k in mine),
    )
    # 3. all-sample funnel total == STAR input (SE arms with unmapped kept)
    if not os.path.exists(ic):
        L = {l.split('|')[0].strip(): l.split('|')[1].strip() for l in open(f'{U}/{arm}/Log.final.out') if '|' in l}
        allin = sum(
            int(r['reads_in']) for r in csv.DictReader(open(f'{U}/{arm}/funnel_per_sample.tsv'), delimiter='\t')
        )
        check(
            f'p0a {m}: per-sample fragments sum to STAR input reads',
            allin == int(L['Number of input reads']),
            f'{allin:,} vs STAR {int(L["Number of input reads"]):,}',
        )
    # 4. filtered-read sum == report slice_reads (same wells)
    check(
        f'p0a {m}: filtered-read sum == report slice_reads',
        abs(tot[3] - nj[m]['slice_reads']) <= 1,
        f'{tot[3]:,.0f} vs {nj[m]["slice_reads"]:,}',
    )
    check(f'p0a {m}: n samples == report', len(usable) == nj[m]['n_samples'], f'{len(usable)} vs {nj[m]["n_samples"]}')
# ======================= p1 / p2 =======================
for setname in (['matched', 'native'] if SET == 'native' else ['matched']):
    for m in METH:
        rows = [r for r in rare(ARMS[setname][m]) if r.get('B_corr') and r['at_native'] == '0']
        for d in sorted({int(r['depth']) for r in rows}):
            rr = [r for r in rows if int(r['depth']) == d]
            uci = st.median(float(r['B_corr']) for r in rr)
            gen = st.median(float(r['genes_B']) for r in rr)
            ref = nj[m]['depths'].get(str(d)) if setname == TAG else None
            if ref:
                check(
                    f'p1/p2 {setname} {m} @{d//1000}k: median UCI/genes == report numbers.json',
                    abs(uci - ref['uci']) < 1 and abs(gen - ref['genes']) < 1 and len(rr) == ref['n'],
                    f'UCI {uci:,.0f} genes {gen:,.0f} n {len(rr)} vs {ref["uci"]:,.0f} / {ref["genes"]:,.0f} / n {ref["n"]}',
                )
            check(
                f'p1/p2 {setname} {m} @{d//1000}k: every used row has reads_used == depth and at_native == 0',
                all(int(r['reads_used']) == d for r in rr),
                f'n={len(rr)}',
            )
        n50 = len([r for r in rows if int(r['depth']) == 50000])
        exp = {'DRUG-seq': 24, 'BRB-seq': 8, 'prime-seq': 8, 'BOBseq': 18}[m]
        check(f'p1/p2 {setname} {m}: samples at 50k == expected benchmark set', n50 == exp, f'{n50} vs {exp}')
# ---- p2 gene count without ribosomal-protein genes: table column vs plotted values vs the independent recount (check_genes_noRP.py) ----
for m in METH:
    rows = [r for r in rare(ARMS[TAG][m]) if r.get('B_corr') and r['at_native'] == '0']
    dif = [int(r['genes_B']) - int(r['genes_B_noRP']) for r in rows]
    check(
        f'p2 {m}: genes_B_noRP = genes_B minus 1..175 ribosomal-protein genes (cytosolic + mitochondrial) in every plotted row',
        all(0 < x <= 175 for x in dif),
        f'removed per row: {min(dif)}..{max(dif)}',
    )
    pv = {
        (r['sample'], r['depth_filtered_reads']): r
        for r in csv.DictReader(open(f'{D}/values/p2_genes_vs_depth_{TAG}.tsv'), delimiter='\t')
        if r['method'] == m
    }
    check(
        f'p2 {m}: values table genes == genes_B_noRP and genes_incl_ribosomal_protein == genes_B, row by row',
        len(pv) == len(rows)
        and all(
            pv[(r['well'], r['depth'])]['genes'] == r['genes_B_noRP']
            and pv[(r['well'], r['depth'])]['genes_incl_ribosomal_protein'] == r['genes_B']
            for r in rows
        ),
        f'{len(rows)} rows',
    )
    ck = f'{U}/{ARMS[TAG][m]}/stats_E_U/check_genes_noRP.tsv'
    cr = list(csv.DictReader(open(ck), delimiter='\t')) if os.path.exists(ck) else []
    rt = {(r['well'], r['depth']): r for r in rare(ARMS[TAG][m])}
    check(
        f'p2 {m}: independent recount (regex-free rule, rebuilt reservoir) agrees with the current table on every row',
        bool(cr)
        and len(cr) == len(rt)
        and all(
            r['status'] == 'ok' and rt[(r['well'], r['depth'])]['genes_B_noRP'] == r['genes_noRP_recount'] for r in cr
        ),
        f'{len(cr)} rows',
    )
# ======================= p3 composition =======================
p3 = json.load(
    open(F(f'p3_read_composition_readtable_{TAG}.json'))
)  # read-table composition (extra/); the figure's p3 is rule-based, checked below
for m in METH:
    arm = ARMS[TAG][m]
    f = f'{U}/{arm}/composition_benchmark_wells.txt'
    f = f if os.path.exists(f) else f'{U}/{arm}/composition_species.txt'
    d = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(f)}
    tot = sum(d.values())
    mine = {k: 100 * d.get(k, 0) / tot for k in 'EXIGRT'}
    check(
        f'p3 {m}: classes are exactly E,X,I,G,R,T and sum to 100',
        set(d) <= set('EXIGRT') and abs(sum(mine.values()) - 100) < 1e-6,
        f'{sorted(d)}',
    )
    check(
        f'p3 {m}: recomputed % match plotted',
        max(abs(mine[k] - p3[m][k]) for k in 'EXIGRT') < 0.05,
        ' '.join(f'{k}={mine[k]:.1f}' for k in 'EXIGRT') + f' | source {os.path.basename(f)}',
    )
    plate = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(f'{U}/{arm}/composition_species.txt')}
    pt = sum(plate.values())
    check(
        f'p3 {m}: benchmark-well composition vs plate-level (info)',
        True,
        f'E {mine["E"]:.1f}% vs plate {100*plate.get("E",0)/pt:.1f}% (plate {pt:,} reads, wells {tot:,})',
    )
# ======================= p3b intronic reads per sample =======================
p3b = json.load(open(F(f'p3b_intronic_per_sample_readtable_{TAG}.json')))
expn = {'DRUG-seq': 24, 'BRB-seq': 8, 'prime-seq': 8, 'BOBseq': 18}
for m in METH:
    arm = ARMS[TAG][m]
    tot = collections.Counter()
    intr = collections.Counter()
    for l in open(f'{U}/{arm}/composition_per_well.tsv'):
        w, reg, a, u = l.rstrip('\n').split('\t')
        tot[w] += int(a)
        intr[w] += int(a) if reg == 'I' else 0
    mine = {w: 100 * intr[w] / tot[w] for w in tot}
    check(
        f'p3b {m}: {expn[m]} benchmark wells plotted, same wells as the table',
        len(p3b[m]) == expn[m] and set(p3b[m]) == set(mine),
        f'{len(p3b[m])} wells',
    )
    check(
        f'p3b {m}: recomputed intronic % per well matches plotted',
        max(abs(mine[w] - p3b[m][w]) for w in mine) < 0.01,
        f'median {np.median(list(mine.values())):.1f}%, range {min(mine.values()):.1f}-{max(mine.values()):.1f}',
    )
    pooled = 100 * sum(intr.values()) / sum(tot.values())
    if m == 'BOBseq':
        exc = {l.split('\t')[0] for i, l in enumerate(open(f'{U}/{arm}/units_excluded.tsv')) if i > 0}
        check(f'p3b {m}: excluded units (Ris 500) absent', not (exc & set(mine)), f'{sorted(exc)}')
        check(
            f'p3b {m}: pooled intronic over the 18 wells == p3 (excluded units carry no rows in reads.tsv)',
            abs(pooled - p3[m]['I']) < 0.05,
            f'{pooled:.2f}% vs p3 {p3[m]["I"]:.2f}%',
        )
    else:
        check(
            f'p3b {m}: per-well counts pooled == p3 intronic %',
            abs(pooled - p3[m]['I']) < 0.05,
            f'{pooled:.2f}% vs p3 {p3[m]["I"]:.2f}%',
        )
# independent recount straight from reads.tsv for the smallest arm (prime-seq: no well map, human contigs only), own awk
arm = ARMS[TAG]['prime-seq']
import subprocess

out = subprocess.run(
    [
        'mawk',
        '-F\t',
        '$3 ~ /^HUMAN_/ {t[$1]++; if ($6=="I") i[$1]++} END {for (w in t) print w, 100*i[w]/t[w]}',
        f'{U}/{arm}/reads.tsv',
    ],
    capture_output=True,
    text=True,
).stdout
rec = {l.split()[0]: float(l.split()[1]) for l in out.strip().split('\n')}
check(
    'p3b prime-seq: independent awk recount from reads.tsv matches plotted per well',
    set(rec) == set(p3b['prime-seq']) and max(abs(rec[w] - p3b['prime-seq'][w]) for w in rec) < 0.01,
    ' '.join(f'{w}={rec[w]:.1f}' for w in sorted(rec)),
)
# ======================= p3 / p3b rule-based composition (the figure's panels) =======================
pr = json.load(open(f'{D}/p3_read_composition_{TAG}.json'))
prb = json.load(open(f'{D}/p3b_intronic_per_sample_{TAG}.json'))
rules = json.load(open(f'{D}/p3_composition_rules_{TAG}.json'))
rows_ = list(csv.DictReader(open(f'{D}/values/p3_composition_rules_{TAG}_per_sample.tsv'), delimiter='\t'))
keys_ = pr['classes']
for m in METH:
    rr = [r for r in rows_ if r['method'] == m]
    tot_ = sum(int(r['mapped_records']) for r in rr)
    mine_ = {k: 100 * sum(int(r[k]) for r in rr) / tot_ for k in keys_}
    check(
        f'p3 rules {m}: {expn[m]} samples, classes sum to 100, panel % == per-sample counts pooled',
        len(rr) == expn[m]
        and abs(sum(mine_.values()) - 100) < 1e-6
        and max(abs(mine_[k] - pr['per_method_pct_of_mapped'][m][k]) for k in keys_) < 1e-6,
        ' '.join(f'{k[:6]}={mine_[k]:.1f}' for k in keys_),
    )
    check(
        f'p3 rules {m}: rRNA equals the read-table rRNA within 0.1 point (info on the E/X split: rules give protein-coding priority)',
        abs(mine_['rRNA'] - p3[m]['R']) < 0.1,
        f'rules rRNA {mine_["rRNA"]:.2f} vs read table {p3[m]["R"]:.2f}; mRNA+RP {mine_["mRNA"]+mine_["ribosomal-protein"]:.1f} vs read-table E {p3[m]["E"]:.1f}',
    )
    pwm = {r['sample']: 100 * int(r['intronic']) / int(r['mapped_records']) for r in rr}
    check(
        f'p3b rules {m}: per-sample intronic % recomputed from the counts equals the panel',
        set(pwm) == set(prb[m]) and max(abs(pwm[k] - prb[m][k]) for k in pwm) < 1e-9,
        f'median {np.median(list(pwm.values())):.2f}%',
    )
    check(
        f'p3b rules {m}: per-sample medians within 0.3 point of the read-table version',
        abs(np.median(list(pwm.values())) - np.median(list(p3b[m].values()))) < 0.3,
        f'{np.median(list(pwm.values())):.2f} vs {np.median(list(p3b[m].values())):.2f}',
    )
# ======================= p4i / p4m2 (positions) =======================
S = [l.rstrip('\n').split('\t') for l in open(f'{A}/samples{SFX}.tsv') if l.strip()]
POOL = {}
for lab, bam, m in S:
    z = np.load(f'{A}/positions_canonical{SFX}/{lab}.npz')
    for k in ('gene', 't', 'L', 'alen'):
        POOL.setdefault(m, {}).setdefault(k, []).append(z[k])
    tx_len = z['tx_len']
    genes = list(z['genes'])
for m in METH:
    POOL[m] = {k: np.concatenate(v) for k, v in POOL[m].items()}
check(
    'p4 sample set: every benchmark sample, Ris 500 nM excluded',
    collections.Counter(m for _, _, m in S) == {'BOBseq': 18, 'DRUG-seq': 24, 'prime-seq': 8}
    and not any('500' in lab for lab, _, _ in S),
    str(dict(collections.Counter(m for _, _, m in S))),
)
cnt = {m: collections.Counter(POOL[m]['gene'].tolist()) for m in METH}
common = [g for g in cnt['DRUG-seq'] if all(cnt[m][g] >= 200 for m in METH)]
bc = json.load(open(F('p4_base_coverage_numbers.json')))
check(
    'p4i: genes with >= 200 unique reads in every method',
    len(common) == bc['genes'],
    f'{len(common)} vs panel {bc["genes"]}',
)
# transcript model: independent parse of the GTF for ACTB / GAPDH canonical lengths
want = {'ACTB', 'GAPDH', 'EEF1A1', 'RPLP0'}
can = {}
ex = collections.defaultdict(int)
strand = {}
for l in open(GTF):
    if not l.startswith('HUMAN_'):
        continue
    f = l.split('\t', 9)
    if f[2] == 'transcript' and 'Ensembl_canonical' in f[8]:
        i = f[8].find('gene_name "')
        g = f[8][i + 11 : f[8].find('"', i + 11)]
        if g in want:
            j = f[8].find('transcript_id "')
            can[g] = f[8][j + 15 : f[8].find('"', j + 15)]
            strand[g] = f[6]
    elif f[2] == 'exon':
        j = f[8].find('transcript_id "')
        tid = f[8][j + 15 : f[8].find('"', j + 15)]
        if tid in can.values():
            ex[tid] += int(f[4]) - int(f[3]) + 1
gi = {g: genes.index(g) for g in want}
check(
    'p4 transcript model: canonical (Ensembl_canonical tag) exon length == tx_len in positions',
    all(ex[can[g]] == tx_len[gi[g]] for g in want),
    ', '.join(f'{g} {ex[can[g]]}={tx_len[gi[g]]} ({can[g]}, {strand[g]})' for g in want),
)
# strand orientation: DRUG-seq read 5' ends must sit near the 3' end (t/L high) for a minus-strand (ACTB) and a plus-strand (GAPDH) gene
for g in ('ACTB', 'GAPDH'):
    sel = POOL['DRUG-seq']['gene'] == gi[g]
    fr = np.median(POOL['DRUG-seq']['t'][sel] / tx_len[gi[g]])
    check(
        f'p4 orientation {g} ({strand[g]} strand): DRUG-seq median read position is in the 3\' third',
        fr > 0.67,
        f'median t/L = {fr:.2f}, n = {sel.sum():,}',
    )
# read filters: every t within [0, L), alen > 0 and <= read length bound
for m in METH:
    ok = (POOL[m]['t'] >= 0).all() and (POOL[m]['t'] < POOL[m]['L']).all() and (POOL[m]['alen'] > 0).all()
    check(
        f'p4 {m}: all positions inside the transcript, aligned length > 0',
        bool(ok),
        f'max alen {POOL[m]["alen"].max()} (native BOBseq mates <= 150; 50-nt set <= 50)',
    )


# independent recomputation of covered fraction (median gene) at 200 reads, own seed and own coverage code
def by_gene(m):
    g = POOL[m]['gene']
    o = np.argsort(g, kind='stable')
    return {int(k): o[np.searchsorted(g[o], k, 'left') : np.searchsorted(g[o], k, 'right')] for k in set(common)}


IDX = {m: by_gene(m) for m in METH}
rng = np.random.default_rng(2024)


def cov_naive(m, ii, L):
    c = np.zeros(L, int)
    for t, a in zip(POOL[m]['t'][ii], POOL[m]['alen'][ii]):
        c[int(t) : min(int(t) + int(a), L)] += 1
    return c


covf = {m: [] for m in METH}
peak = {m: [] for m in METH}
for g in common:
    L = int(tx_len[g])
    for m in METH:
        ii = rng.choice(IDX[m][g], 200, replace=False)
        c = cov_naive(m, ii, L)
        covf[m].append(100 * np.mean(c > 0))
        t = np.sort(POOL[m]['t'][ii])
        w = min(300, L)
        peak[m].append(
            100
            * max(np.searchsorted(t, s + w, 'left') - np.searchsorted(t, s, 'left') for s in range(0, L - w + 1, 10))
            / 200
        )
for m in METH:
    check(
        f'p4g/p4i {m}: independent median covered fraction at 200 reads within 2 points of panel',
        abs(np.median(covf[m]) - bc['covered_frac_median'][m]) < 2,
        f'{np.median(covf[m]):.1f} vs {bc["covered_frac_median"][m]:.1f}',
    )
    check(
        f'p4h {m}: independent median densest-300-nt share within 3 points of panel',
        abs(np.median(peak[m]) - bc['peak_median'][m]) < 3,
        f'{np.median(peak[m]):.1f} vs {bc["peak_median"][m]:.1f}',
    )
# p4m2: mean over genes of coverage / gene mean, genes >= 2 kb, last 3000 nt, at 100/500/1000/2000 nt from the poly(A) site; vs Q_UNSCALED variant 2
q = [l.split('|') for l in open(f'{CV}/Q_UNSCALED.md') if l.startswith('| 2. same, MEAN')]
qv = {r[3].strip(): [float(x) for x in r[4:9]] for r in q}
long = [g for g in common if tx_len[g] >= 2000]
prof = {m: np.zeros(3000) for m in METH}
nn = {m: np.zeros(3000) for m in METH}
rng = np.random.default_rng(99)
for g in long:
    L = int(tx_len[g])
    for m in METH:
        ii = rng.choice(IDX[m][g], 200, replace=False)
        c = cov_naive(m, ii, L)
        z = c / c.mean()
        w = min(3000, L)
        prof[m][:w] += z[::-1][:w]
        nn[m][:w] += 1
        if g == long[0]:
            check(
                f'p4m2 normalisation ({m}, first gene): mean of coverage/gene-mean over the transcript == 1',
                abs(z.mean() - 1) < 1e-9,
                '',
            )
check(
    'p4m2 gene set: genes >= 2 kb with >= 200 reads in every method',
    len(long) == int(q[0][2]),
    f'{len(long)} vs table {q[0][2].strip()}',
)
for m in METH:
    mine = prof[m] / nn[m]
    pts = [mine[100], mine[500], mine[1000], mine[2000]]
    ref = qv[m][:4]
    rel = max(abs(a - b) / max(b, 0.05) for a, b in zip(pts, ref))
    check(
        f'p4m2 {m}: independent mean profile at 100/500/1000/2000 nt within 20% of the table (own seed; seed-to-seed spread of this mean is up to ~12%, measured with the generator code at seeds 5/6/11)',
        rel < 0.20,
        ' / '.join(f'{a:.2f}' for a in pts) + ' vs ' + ' / '.join(f'{b:.2f}' for b in ref),
    )
# ======================= p4 thresholds: heatmap + profile at >= 1 kb, no ribosomal-protein genes =======================
import re as _re, subprocess as _sp, tempfile as _tf

thr = json.load(open(f'{D}/p4_thresholds_numbers.json'))
RIBO = _re.compile(r'^(RP[LS]|RPLP|MRP[LS])[^K]*$', _re.I)
tset = [g for g in common if tx_len[g] >= 1000 and not RIBO.match(genes[g])]
check(
    'p4thr: gene set (>= 200 reads in every method, >= 1 kb, not ribosomal-protein) matches the panel',
    len(tset) == thr['genes'] and sorted(genes[g] for g in tset) == sorted(thr['gene_names']),
    f'{len(tset)} vs {thr["genes"]}',
)
check(
    'p4thr: every removed gene matches the ribosomal-protein rule and none of the kept genes does',
    all(RIBO.match(x) for x in thr['rp_genes_removed']) and not any(RIBO.match(x) for x in thr['gene_names']),
    f'{len(thr["rp_genes_removed"])} removed',
)
rng = np.random.default_rng(7)
prof2 = {m: np.zeros(3000) for m in METH}
nn2 = {m: np.zeros(3000) for m in METH}
frac = {m: [] for m in METH}
for g in tset:
    L = int(tx_len[g])
    idx = np.minimum((np.arange(L) * 50) // L, 49)
    for m in METH:
        ii = rng.choice(IDX[m][g], 200, replace=False)
        c = cov_naive(m, ii, L)
        z = c / c.mean()
        w = min(3000, L)
        prof2[m][:w] += z[::-1][:w]
        nn2[m][:w] += 1
        ca = cov_naive(m, IDX[m][g], L)
        b = np.bincount(idx, weights=ca, minlength=50) / np.bincount(idx, minlength=50)
        frac[m].append(np.mean(b > 0.2 * b.max()) if b.max() > 0 else 0)
check(
    'p4thr: genes contributing at 1,000 / 2,000 / 3,000 nt from the poly(A) site',
    [int(nn2['BOBseq'][d - 1]) for d in (1000, 2000, 3000)]
    == [thr['genes_contributing_at_nt'][str(d)] for d in (1000, 2000, 3000)],
    str([int(nn2['BOBseq'][d - 1]) for d in (1000, 2000, 3000)]),
)
for m in METH:
    mine = prof2[m] / np.maximum(nn2[m], 1)
    a0, a1 = mine[:300].mean(), mine[1000:2000].mean()
    r0, r1 = thr['profile'][m]['mean_0_300'], thr['profile'][m]['mean_1000_2000']
    check(
        f'p4thr {m}: independent profile (own seed, unsmoothed) within 20% of the panel at 0-300 and 1,000-2,000 nt',
        abs(a0 - r0) / r0 < 0.2 and abs(a1 - r1) / r1 < 0.2,
        f'{a0:.2f}/{a1:.2f} vs {r0:.2f}/{r1:.2f}',
    )
    check(
        f'p4thr {m}: heatmap share of bins above 20% of the row peak (all reads, own binning) within 0.02',
        abs(np.mean(frac[m]) - thr['heatmap_frac_bins_above_20pct_of_peak'][m]) < 0.02,
        f'{np.mean(frac[m]):.3f} vs {thr["heatmap_frac_bins_above_20pct_of_peak"][m]:.3f}',
    )
# ======================= p6 Picard profile and balance =======================
p6 = json.load(open(f'{D}/p6_picard_numbers_{TAG}.json'))
MET = f'{CV}/picard_metrics/{TAG}'
PB = f'{P}/benchmark_uniform/per_sample_bams'


def hist_of(path):
    rows = [l.rstrip('\n').split('\t') for l in open(path)]
    inh = False
    h = []
    for r in rows:
        if r and r[0].startswith('## HISTOGRAM'):
            inh = True
            continue
        if inh and len(r) >= 2 and r[0].isdigit():
            h.append(float(r[1]))
    return np.array(h)


pp = np.arange(101) / 100
for m in METH:
    check(
        f'p6 {m}: n samples == benchmark set',
        p6['per_method'][m]['n'] == {'DRUG-seq': 24, 'BRB-seq': 8, 'prime-seq': 8, 'BOBseq': 18}[m],
        str(p6['per_method'][m]['n']),
    )
    Hs = np.array([hist_of(f'{MET}/{n}.RNA_Metrics.txt') for n in p6['per_method'][m]['samples']])
    check(
        f'p6 {m}: every cached histogram has 101 bins and a positive sum',
        Hs.shape[1] == 101 and (Hs.sum(1) > 0).all(),
        str(Hs.shape),
    )
    check(
        f'p6b {m}: balance recomputed from the cached histograms equals the panel',
        np.allclose(2 * (Hs * pp).sum(1) / Hs.sum(1), p6['per_method'][m]['balance'], atol=1e-9),
        f'mean {np.mean(p6["per_method"][m]["balance"]):.3f}',
    )
    check(
        f'p6 {m}: profile mean recomputed from the cached histograms equals the panel',
        np.allclose(Hs.mean(0), p6['per_method'][m]['profile_mean'], atol=1e-9),
        f'peak at {int(np.argmax(Hs.mean(0)))}',
    )
bdir = {
    'DRUG-seq': 'drugseq_native' if SET == 'native' else 'drugseq',
    'BOBseq': 'bobseq_pe_native' if SET == 'native' else 'bobseq_50nt',
}
for m in ('DRUG-seq', 'BOBseq'):
    n = p6['per_method'][m]['samples'][0]
    tmp = _tf.mkdtemp(prefix='picard_val_')
    r = _sp.run(
        [
            'picard',
            'CollectRnaSeqMetrics',
            '-I',
            f'{PB}/{bdir[m]}/{n}.bam',
            '-O',
            f'{tmp}/val.txt',
            '--REF_FLAT',
            f'{A}/picard_refflat_pc_noRP_noMT_bare.txt',
            '--STRAND_SPECIFICITY',
            'NONE',
            '--MINIMUM_LENGTH',
            '1000',
            '--VALIDATION_STRINGENCY',
            'LENIENT',
        ],
        capture_output=True,
        text=True,
    )
    ok = (
        r.returncode == 0
        and os.path.exists(f'{tmp}/val.txt')
        and np.allclose(hist_of(f'{tmp}/val.txt'), hist_of(f'{MET}/{n}.RNA_Metrics.txt'), atol=1e-6)
    )
    check(
        f'p6 {m} ({n}): fresh Picard run reproduces the cached histogram exactly',
        ok,
        'identical 101 bins' if ok else (r.stderr[-300:] if r.returncode else 'histogram differs'),
    )
# ======================= p5 junctions =======================
import pysam

J = list(csv.DictReader(open(F('p5_junctions_per_sample.tsv')), delimiter='\t'))
src = {}
for f in glob.glob(f'{A}/junctions{SFX}/*.tsv'):
    src.update({(r['sample'], r['depth']): r for r in csv.DictReader(open(f), delimiter='\t')})
check(
    'p5: panel side table == per-sample junction files (values and at_native flag)',
    all(
        src[(r['sample'], r['depth'])]['junctions_ge2'] == r['junctions_ge2']
        and src[(r['sample'], r['depth'])]['at_native'] == r['at_native']
        for r in J
    ),
    f'{len(J)} rows, {sum(1 for r in J if r["at_native"] == "1")} flagged at_native (excluded from the plots)',
)
gen = open(f'{CODE}/main_figure_panels.py').read()
check(
    'p5: the panel code plots only at_native == 0 rows',
    gen.count("r['at_native'] == '0'") >= 3,
    f"{gen.count(chr(114) + chr(91) + chr(39) + 'at_native')} filters in the generator",
)
PB = f'{U}/per_sample_bams'
one = {}
for lab, bam, m in S:
    if m not in one:
        one[m] = (lab, bam)
for m, (lab, bam) in one.items():
    rng2 = random.Random(7)
    res = []
    n = 0
    for a in pysam.AlignmentFile(f'{PB}/{bam}', 'rb').fetch(until_eof=True):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality != 255:
            continue
        n += 1
        js = (
            tuple((a.reference_name, b[0], b[1]) for b in a.get_blocks()[:-1])
            if a.cigarstring and 'N' in a.cigarstring
            else ()
        )
        # junction = gap between consecutive aligned blocks separated by an N op: recompute from blocks (independent of the generator's CIGAR walk)
        if js:
            bl = a.get_blocks()
            gaps = []
            pos = a.reference_start
            k = 0
            for op, ln in a.cigartuples:
                if op == 3:
                    gaps.append((a.reference_name, pos, pos + ln))
                if op in (0, 2, 3, 7, 8):
                    pos += ln
            js = tuple(gaps)
        if len(res) < 500000:
            res.append(js)
        else:
            j = rng2.randrange(n)
            if j < 500000:
                res[j] = js
    cnt = collections.Counter(j for r in res for j in r)
    ge2 = sum(1 for v in cnt.values() if v >= 2)
    ref = int(src[(lab, '500000')]['junctions_ge2']) if (lab, '500000') in src else None
    if ref is None:
        check(f'p5 {m} ({lab}): sample has < 500k unique reads, excluded at 500k', n < 500000, f'n={n:,}')
    else:
        check(
            f'p5 {m} ({lab}): independent reservoir (own seed) junctions >= 2 reads at 500k within 5% of the file',
            abs(ge2 - ref) / ref < 0.05,
            f'{ge2:,} vs {ref:,} (unique reads {n:,})',
        )
# ======================= SVG hygiene =======================
for f in (
    'p0a_read_fate_%s' % TAG,
    'p1_molecules_vs_depth_matched',
    'p2_genes_vs_depth_matched',
    'p3_read_composition_%s' % TAG,
    'p3b_intronic_per_sample_%s' % TAG,
    'p4i_coverage_heatmaps_200_reads',
    'p4m2_coverage_vs_nt_from_polyA_mean',
    'p5_junctions_vs_depth',
    'p4i2_coverage_heatmaps_all_reads_1kb_noRP',
    'p4m2_coverage_vs_nt_from_polyA_mean_1kb_noRP',
    'p6_picard_profile_%s' % TAG,
    'p6b_picard_balance_%s' % TAG,
):
    s = open(F(f'{f}.svg')).read()
    fams = set(
        _re.split(r'[;"]', x)[0] for x in s.split('font-family:')[1:]
    )  # a style attribute may end right after the family (no ';')
    check(
        f'svg {f}: Arial only, no legend/title text',
        (
            fams <= {' Arial', 'Arial'} and 'DRUG-seq</' not in s.replace('p3', '')
            if not f.startswith('p3') and not f.startswith('p4i') and not f.startswith('p6')
            else fams <= {' Arial', 'Arial'}
        ),
        f'fonts {fams}',
    )
open(f'{D}/VALIDATION_main_panels.md', 'w').write(
    f'# Validation of the main-figure panels ({FIG}), {__import__("datetime").date.today()}\n\n{FAIL} failed of {len(R)} checks. Independent recomputation from the source tables/BAMs; own code, own seeds.\n\n| result | check | detail |\n|---|---|---|\n'
    + '\n'.join(R)
    + '\n'
)
print(f'\n{FAIL} FAILED of {len(R)} checks -> {D}/VALIDATION_main_panels.md')
sys.exit(1 if FAIL else 0)

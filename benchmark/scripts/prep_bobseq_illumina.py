#!/usr/bin/env python3
"""BOBseq Illumina reads -> the benchmark's common input format, `<readid>_<BOBCODE>_<UMI>`.

This is the ONLY BOBseq-specific step in the benchmark. Everything after it (truncate to
50 nt, STAR, read table, slices, molecule definitions) is the same code the competitor arms
go through, which is what makes the comparison a fact rather than a claim.

Deliberately different from pipeline/scripts/prep_reads_illumina.py:
  * no poly(A) trimming and no minimum-insert drop beyond a sanity floor: the competitor
    arms get neither, and the common path truncates to 50 nt anyway
  * the non-templated G-run is stripped so the read starts at cDNA like every other arm
  * UMI is chosen from the chemistry, never hardcoded:
        R1_26N_polydT  -> 32N (26 nt R1 + 6 nt TSO)   or 6N with --umi 6N
        R1long-6N / 9N -> 6N  (TSO only: the R1 hexamer is the priming site, not a UMI)
        --umi full     -> the whole degenerate stretch between the bobcode and the G-run
                          (6N + spacer: 8 / 11 / 14 nt by code set, 7 for the 11-mer; every
                          position measured degenerate on the 24-plex; this is what the
                          pipeline uses)
  * bobcode: positional at R2[0:L], at most one substitution, ties dropped as ambiguous,
    the same rule the competitor whitelists use
"""

import argparse, gzip, json, sys, collections


def opener(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)


ap = argparse.ArgumentParser()
ap.add_argument('r1')
ap.add_argument('r2')
ap.add_argument('out')
ap.add_argument('--run-json', required=True, help='run configuration JSON (bobcodes, offsets, RT primers), e.g. config/run_24plex.json')
ap.add_argument('--umi', default='auto', choices=['auto', '6N', '32N', 'full'])
ap.add_argument('--max-mm', type=int, default=1)
ap.add_argument('--min-insert', type=int, default=30)
ap.add_argument('--stats', default=None)
a = ap.parse_args()

cfg = json.load(open(a.run_json))
codes = {k.upper(): v for k, v in (cfg.get('tso_species_map') or {}).items()}
offs = {k.upper(): int(v) for k, v in (cfg.get('grun_offsets') or {}).items()}
if not codes:
    sys.exit('prep_bobseq: no tso_species_map in the run json')
rt = ' '.join(cfg.get('rt_primers_used') or []).lower()
if a.umi == 'auto':
    umi_mode = '32N' if '26n' in rt else ('6N' if any(t in rt for t in ('6n', '9n', 'random', 'hexamer')) else None)
    if umi_mode is None:
        sys.exit(f'prep_bobseq: cannot place the UMI for rt_primers_used={rt!r}; pass --umi')
else:
    umi_mode = a.umi
by_len = sorted({len(c) for c in codes}, reverse=True)


def call(seq):
    """(code, mismatches) or (None, None); longest codes first, ties ambiguous."""
    for L in by_len:
        w = seq[:L]
        if len(w) < L:
            continue
        best, bd, tie = None, a.max_mm + 1, False
        for c in (x for x in codes if len(x) == L):
            d = sum(1 for x, y in zip(w, c) if x != y)
            if d < bd:
                best, bd, tie = c, d, False
            elif d == bd:
                tie = True
        if best is not None and bd <= a.max_mm:
            return (None, None) if tie else (best, bd)
    return None, None


st = collections.Counter()
with opener(a.r1) as f1, opener(a.r2) as f2, gzip.open(a.out, 'wt', compresslevel=1) as out:
    while True:
        h1 = f1.readline()
        if not h1:
            break
        s1 = f1.readline().rstrip('\n')
        f1.readline()
        f1.readline()
        h2 = f2.readline()
        s2 = f2.readline().rstrip('\n')
        f2.readline()
        q2 = f2.readline().rstrip('\n')
        st['reads'] += 1
        code, mm = call(s2)
        if code is None:
            st['no_bobcode'] += 1
            continue
        st['mm0' if mm == 0 else 'mm1'] += 1
        L = len(code)
        umi2 = s2[L : L + (offs.get(code, 8) if umi_mode == 'full' else 6)]
        umi = umi2 if umi_mode in ('6N', 'full') else s1[:26] + umi2
        if len(umi) < (offs.get(code, 8) if umi_mode == 'full' else 6 if umi_mode == '6N' else 32) or 'N' in umi:
            st['umi_bad'] += 1
            continue
        gstart = L + offs.get(code, 8)
        g = 0
        while gstart + g < len(s2) and s2[gstart + g] == 'G':
            g += 1
        ins, q = s2[gstart + g :], q2[gstart + g :]
        if len(ins) < a.min_insert:
            st['too_short'] += 1
            continue
        rid = h2[1:].split()[0]
        out.write(f'@{rid}_{code}_{umi}\n{ins}\n+\n{q}\n')
        st['written'] += 1
        st[f'well:{code}'] += 1
st['umi_mode'] = umi_mode
sys.stderr.write(
    '  '
    + ' | '.join(f'{k}={v:,}' if isinstance(v, int) else f'{k}={v}' for k, v in st.items() if not k.startswith('well:'))
    + '\n'
)
if a.stats:
    json.dump(dict(st), open(a.stats, 'w'), indent=1)

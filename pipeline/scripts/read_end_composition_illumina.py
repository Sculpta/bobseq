#!/usr/bin/env python3
"""Standardized read START / distal-END composition for ONT cDNA libraries.

For a barcode's FASTQ it reports, in four views:
  1a. true FIRST element   -- literal 5' terminus (ONT adapter included)
  1b. element AFTER ONT    -- first library element past the ONT block (orientation)
  2a. true LAST element    -- literal 3' terminus (ONT adapter included)
  2b. element BEFORE ONT   -- last library element before the ONT block

Standard category groupings (fixed across chemistries):
    ONT adapter + barcode   ONT native adapter + AAGGTTAA/CAGCACCT flanks + the
                            (auto-detected, well-specific) 24 nt native barcode
    P5 / i5                 P5 (+ auto-detected i5 index)          -- RT/P5 end
    P7 / i7                 P7 (+ auto-detected i7 index)          -- P7 end
    TSO / R2 (...)          TruSeq Read2 handle + species 7-mer + GGG  -- P7 end
    R1 (TruSeq1)            TruSeq Read1 handle                    -- RT/P5 end
    polydT                  poly-T/A tract
    likely cDNA             no adapter near this terminus -> inferred insert
                            (NOT verified against an alignment, hence "likely")
    unresolved              adapter block detectable nearby but terminus unmatched

MATCHING.  ONT reads carry ~5% substitution AND indel error, and quality collapses
over the last ~20 bp, so a single long exact match is far too insensitive and a
Hamming-fuzzy match cannot cross indels (and, on a 7 bp motif, invents matches).
Instead every element is reduced to DENSE 9-mer exact seeds (stride 1, both
strands) held in a hash index: one substitution or indel kills only the seeds
spanning it while the rest still hit, and lookup is O(read length), not
O(#seeds).  A read that stops PART-WAY into an element (no full 9-mer at the
terminus) is caught by a separate terminal-stub check.  Seeds shared by >1 group
(e.g. the CTCTTCCGATCT anchor common to R1 and R2) are dropped -- they cannot
discriminate and would assign by dict order.

Sample-specific sequences (ONT native barcode, i5, i7) are LEARNED from the data,
never hardcoded, so this runs unchanged on any plate well.  Diagnostics (detected
sequences, failures, module constants) travel with the result for provenance.

Usage:  python3 read_end_composition_illumina.py <fastq> [chem]   (default chem = 24plex)
"""
import os
import re
import sys
import json
import collections

# ---- module constants ------------------------------------------------------
ONT_ADAPTER = 'AATGTACTTCGTTCAGTTACGTATTGCT'   # ONT native adapter
NB_FLANK_5 = 'AAGGTTAA'                          # native-barcode 5' flank
NB_FLANK_3 = 'CAGCACCT'                          # native-barcode 3' flank (library-proximal)

ONT_GROUP = 'ONT adapter + barcode'
TSO_GROUP = 'TSO / R2 (TruSeq2 + 7mer + GGG)'
P5_GROUP = 'P5 / i5'
P7_GROUP = 'P7 / i7'
R1_GROUP = 'R1 (TruSeq1)'
POLYT = 'polydT'
LIKELY_CDNA = 'likely cDNA'
UNRESOLVED = 'unresolved (degraded adapter)'

# ORDER is the canonical stacking / legend order and the complete category set.
ORDER = [ONT_GROUP, P5_GROUP, P7_GROUP, TSO_GROUP, R1_GROUP,
         POLYT, LIKELY_CDNA, UNRESOLVED]

# Category colours: hue encodes molecule END (magenta = P7 side, green = P5/RT
# side); grey = ONT, blue = likely cDNA, black = unresolved.
HUES = {ONT_GROUP: '#9e9e9e', P7_GROUP: '#c51b8a', TSO_GROUP: '#f78fc2',
        P5_GROUP: '#00701f', R1_GROUP: '#4bb062', POLYT: '#a8ddb5',
        LIKELY_CDNA: '#1f6fd0', UNRESOLVED: '#000000'}

# Per-chemistry element sequences. `groups` = universal motifs (never a
# sample-specific index); `indexes` = where to LEARN i5/i7 from the data.
CHEMISTRIES = {
    # 24plex: FULL Illumina; the TSO backbone IS the TruSeq Read2 handle.
    '24plex': {
        'groups': {
            P5_GROUP:  ['AATGATACGGCGACCACCGAGATCTACAC'],                 # P5
            P7_GROUP:  ['CAAGCAGAAGACGGCATACGAGAT'],                     # P7
            # Handle motif only. Species 7-mers are below the 9-nt seed length and hit by chance in
            # the terminal windows (measured: true-start TSO 0.4 -> 5.5% with one run's codes), and a
            # fixed list belongs to one run, not the chemistry.
            TSO_GROUP: ['GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'],           # full R2/TSO
            R1_GROUP:  ['ACACTCTTTCCCTACACGACGCTCTTCCGATCT'],            # full R1
        },
        'indexes': {
            'i5': {'anchor': 'AATGATACGGCGACCACCGAGATCTACAC', 'side': 'after',
                   'group': P5_GROUP, 'len': 12},
            'i7': {'anchor': 'CAAGCAGAAGACGGCATACGAGAT', 'side': 'after',
                   'group': P7_GROUP, 'len': 12},
        },
    },
    # tso_half_illumina: HALF Illumina (P5 side only); plain TSO primer on the other end.
    'tso_half_illumina': {
        'groups': {
            P5_GROUP:  ['AATGATACGGCGACCACCGAGATCTACAC'],
            P7_GROUP:  ['CAAGCAGAAGACGGCATACGAGAT'],                     # confirm absent
            TSO_GROUP: ['AAGCAGTGGTATCAACGCAG', 'ACCTTGA'],
            R1_GROUP:  ['CTACACGACGCTCTTCCGATCT'],
        },
        'indexes': {
            'i5': {'anchor': 'AATGATACGGCGACCACCGAGATCTACAC', 'side': 'after',
                   'group': P5_GROUP, 'len': 10},
        },
    },
}

K = 9                    # dense seed length
W = 60                   # terminal window scanned for the true-end views
HEAD, SPAN = 90, 130     # start-view: ONT search window, downstream library window
TAIL = 110               # distal-boundary search window (ONT block is <~68 bp)
STUB_MIN = {ONT_GROUP: 5}   # min terminal-stub length per group
STUB_MIN_DEFAULT = 8
CHEM_MIN_ANCHOR = 0.30      # <this fraction of reads with a P5/P7 anchor -> warn

POLY = re.compile(r'T{10,}|A{10,}')
GRUN = re.compile(r'G{3,}|C{3,}')


def rc(s):
    return s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def dense_seeds(motifs):
    """All overlapping K-mers (stride 1, both strands); motifs < K kept whole."""
    out = set()
    for m in motifs:
        for ref in (m, rc(m)):
            if len(ref) < K:
                out.add(ref)
            else:
                for i in range(len(ref) - K + 1):
                    out.add(ref[i:i + K])
    return out


# ---- learning sample-specific sequences from the data ----------------------
def detect_between(fastq, f5, f3, lo, hi, sample=8000, minfrac=0.15):
    """Consensus of the constant sequence between two fixed flanks (both strands).
    Used for the ONT native barcode. Returns the sequence or None."""
    mid = collections.Counter()
    with _reads(fastq, sample) as it:
        for s in it:
            for hay in (s[:160], rc(s[-160:])):
                i = hay.find(f5); j = hay.find(f3)
                if 0 <= i and i + len(f5) <= j:
                    g = j - (i + len(f5))
                    if lo <= g <= hi:
                        mid[hay[i + len(f5): j]] += 1
                        break
    if not mid:
        return None
    seq, c = mid.most_common(1)[0]
    return seq if c >= minfrac * sum(mid.values()) else None


def detect_adjacent(fastq, anchor, side, length=12, sample=8000, minfrac=0.15):
    """Consensus of `length` bases immediately after/before an anchor (both
    strands). Used for i5 (after P5) and i7 (before P7). Returns seq or None."""
    cnt = collections.Counter()
    with _reads(fastq, sample) as it:
        for s in it:
            for strand in (s, rc(s)):
                i = strand.find(anchor)
                if i < 0:
                    continue
                seg = (strand[i + len(anchor): i + len(anchor) + length]
                       if side == 'after' else strand[max(0, i - length): i])
                if len(seg) == length:
                    cnt[seg] += 1
                break
    if not cnt:
        return None
    seq, c = cnt.most_common(1)[0]
    return seq if c >= minfrac * sum(cnt.values()) else None


class _reads:
    """Iterate up to `limit` read sequences from a FASTQ (context manager)."""
    def __init__(self, fastq, limit):
        self.fastq, self.limit = fastq, limit
    def __enter__(self):
        self.fh = open(self.fastq)
        return self._gen()
    def _gen(self):
        n = 0
        while n < self.limit:
            if not self.fh.readline():
                break
            s = self.fh.readline().strip(); self.fh.readline(); self.fh.readline()
            n += 1
            yield s
    def __exit__(self, *a):
        self.fh.close()


# ---- classifier ------------------------------------------------------------
def build(chem, detected=None):
    """Compile a classifier for `chem`, folding in any detected sample-specific
    sequences (native barcode -> ONT group; i5/i7 -> their P5/P7 group)."""
    detected = detected or {}
    spec = CHEMISTRIES[chem]
    groups = {g: list(ms) for g, ms in spec['groups'].items()}
    for name, ix in spec.get('indexes', {}).items():
        if detected.get(name):
            groups[ix['group']].append(detected[name])
    ont = [ONT_ADAPTER, NB_FLANK_5, NB_FLANK_3]
    if detected.get('native_bc'):
        ont.append(detected['native_bc'])
    groups[ONT_GROUP] = ont

    seed_sets = {g: dense_seeds(ms) for g, ms in groups.items()}
    seen = collections.Counter(s for ss in seed_sets.values() for s in ss)
    kidx, shorts = {}, []
    for g, ss in seed_sets.items():
        for s in ss:
            if seen[s] != 1:                      # non-discriminative -> drop
                continue
            (kidx.__setitem__(s, g) if len(s) == K else shorts.append((s, g)))
    return {
        'chem': chem,
        'motifs': groups,                         # full sequences per group
        'kidx': kidx,                             # discriminative K-mer -> group
        'shorts': shorts,                         # (motif<K, group)
        'stub': {g: [(m, rc(m)) for m in ms] for g, ms in groups.items()},
        'ont_seeds': sorted(seed_sets[ONT_GROUP]),
        'group_seeds': {g: ss for g, ss in seed_sets.items()},
    }


def _scan(win, clf, from_end, skip=None):
    """Outermost recognised element in `win` as (distance_from_terminus, group).
    K-mer hits via one hash-lookup pass; short motifs by find/rfind."""
    kidx = clf['kidx']; L = len(win); bd, bg = 10 ** 9, None
    if from_end:
        for i in range(L - K, -1, -1):            # rightmost hit = outermost
            g = kidx.get(win[i:i + K])
            if g is not None and g != skip:
                bd, bg = L - (i + K), g
                break
    else:
        for i in range(0, L - K + 1):             # leftmost hit = outermost
            g = kidx.get(win[i:i + K])
            if g is not None and g != skip:
                bd, bg = i, g
                break
    for q, g in clf['shorts']:
        if g == skip:
            continue
        j = win.rfind(q) if from_end else win.find(q)
        if j >= 0:
            d = (L - (j + len(q))) if from_end else j
            if d < bd:
                bd, bg = d, g
    return bd, bg


def terminal_stub(seq, stub, three_prime):
    """Read stops part-way INTO an element: its terminal k bases equal that
    element's leading (3') / trailing (5') k bases. Longest, per-group-min."""
    bk, bn = 0, None
    for g, refs in stub.items():
        mn = STUB_MIN.get(g, STUB_MIN_DEFAULT)
        for m, m_rc in refs:
            for ref in (m, m_rc):
                for k in range(min(len(ref), 24), mn - 1, -1):
                    if k <= bk:
                        break
                    obs = seq[-k:] if three_prime else seq[:k]
                    if len(obs) < k:                       # length guard
                        continue
                    exp = ref[:k] if three_prime else ref[-k:]
                    if sum(a != b for a, b in zip(obs, exp)) <= (1 if k >= 12 else 0):
                        bk, bn = k, g
                        break
    return bn


def _residual(win, clf):
    """Nothing matched: adapter still nearby -> unresolved; else -> likely cDNA."""
    ont = clf['group_seeds'][ONT_GROUP]
    return UNRESOLVED if any(q in win for q in ont) else LIKELY_CDNA


def _terminal(seq, clf, three_prime):
    """Identity of one literal terminus (scan-first, stub fallback)."""
    win = seq[-W:] if three_prime else seq[:W]
    d, g = _scan(win, clf, three_prime)
    m = list(POLY.finditer(win))
    if m:
        mm = m[-1] if three_prime else m[0]
        dd = (len(win) - mm.end()) if three_prime else mm.start()
        if dd < d:
            d, g = dd, POLYT
    if g is not None and d == 0:                  # element sits AT the terminus
        return g
    stub = terminal_stub(seq, clf['stub'], three_prime)  # partial element at end
    if stub:
        return stub
    if g is not None:
        return g
    gr = list(GRUN.finditer(win))
    if gr:
        mg = gr[-1] if three_prime else gr[0]
        if (len(win) - mg.end() if three_prime else mg.start()) <= 6:
            return TSO_GROUP                       # GGG folded into the TSO block
    return _residual(seq[-80:] if three_prime else seq[:80], clf)


def classify_5p(seq, clf):
    return _terminal(seq, clf, False)


def classify_end(seq, clf):
    return _terminal(seq, clf, True)


def end_boundary(seq, clf):
    """(index where the LIBRARY ends, ONT block present at 3'?). Trims the ONT
    block so 2b reports the adjacent library element, not the adapter."""
    tail = seq[-TAIL:]; off = len(seq) - len(tail)
    for q in (rc(NB_FLANK_3), NB_FLANK_3):
        i = tail.find(q)
        if i >= 0:
            return off + i, True
    b = None
    for q in clf['ont_seeds']:
        i = tail.find(q)
        if i >= 0 and (b is None or i < b):
            b = i
    return (off + b, True) if b is not None else (len(seq), False)


def classify_start(seq, clf):
    """(has_ont, has_flank, first library group after the ONT block)."""
    head = seq[:HEAD]; hit = None
    for q in clf['ont_seeds']:
        i = head.find(q)
        if i >= 0 and (hit is None or i < hit):
            hit = i + len(q)                      # innermost adapter edge
    if hit is None:
        return False, False, _residual(seq[:80], clf)
    aft = seq[hit: hit + SPAN]; fe = None
    for q in (NB_FLANK_3, rc(NB_FLANK_3)):
        i = aft.find(q)
        if i >= 0 and (fe is None or i < fe):
            fe = i + len(q)
    win = seq[hit + (fe or 0):][:SPAN]
    d, g = _scan(win, clf, False, skip=ONT_GROUP)
    m = POLY.search(win)
    if m and m.start() < d:
        g = POLYT
    return True, fe is not None, (g or _residual(win, clf))


# ---- diagnostics + driver --------------------------------------------------
def _constants(chem):
    return {
        'ont_adapter': ONT_ADAPTER, 'nb_flank_5': NB_FLANK_5,
        'nb_flank_3': NB_FLANK_3, 'seed_k': K, 'window_true_end': W,
        'window_start': HEAD, 'chemistry': chem,
        'group_motifs': CHEMISTRIES[chem]['groups'],
    }


def validate_chemistry(fastq, clf, sample=3000):
    """Fraction of sampled reads carrying a P5 or P7 seed -- a sanity check that
    the chemistry actually matches the data."""
    anchors = clf['group_seeds'][P5_GROUP] | clf['group_seeds'][P7_GROUP]
    n = hit = 0
    with _reads(fastq, sample) as it:
        for s in it:
            n += 1
            if any(q in s for q in anchors):
                hit += 1
    return hit / n if n else 0.0


def analyse(fastq, chem='24plex'):
    if chem not in CHEMISTRIES:
        raise SystemExit(f"unknown chemistry '{chem}'; choose {list(CHEMISTRIES)}")
    warnings = []
    native_bc = detect_between(fastq, NB_FLANK_5, NB_FLANK_3, 20, 28)
    if not native_bc:
        warnings.append("ONT native barcode NOT detected -- the ~32 nt barcode block is "
                        "unmodeled, so 'unresolved' will be inflated. Is this library "
                        "natively barcoded (or is the flank sequence different)?")
    detected = {'native_bc': native_bc}
    for name, ix in CHEMISTRIES[chem].get('indexes', {}).items():
        seq = detect_adjacent(fastq, ix['anchor'], ix['side'], ix.get('len', 12))
        detected[name] = seq
        if not seq:
            warnings.append(f"{name} index not detected; classification still anchors on "
                            f"{ix['group']}.")
    clf = build(chem, detected)

    frac = validate_chemistry(fastq, clf)
    if frac < CHEM_MIN_ANCHOR:
        warnings.append(f"POSSIBLE CHEMISTRY MISMATCH: only {100*frac:.0f}% of sampled reads "
                        f"carry a P5/P7 anchor for chem='{chem}'. Wrong chemistry argument?")

    ts = collections.Counter(); ao = collections.Counter()
    te = collections.Counter(); bo = collections.Counter()
    n = n_ont = n_flank = n_ont_end = 0
    with open(fastq) as fh:
        while True:
            if not fh.readline():
                break
            s = fh.readline().strip(); fh.readline(); fh.readline()
            n += 1
            ts[classify_5p(s, clf)] += 1
            te[classify_end(s, clf)] += 1
            b, had = end_boundary(s, clf); n_ont_end += had
            lib = s[:b]
            bo[classify_end(lib, clf) if len(lib) >= 20 else UNRESOLVED] += 1
            ont, fl, g = classify_start(s, clf)
            n_ont += ont; n_flank += fl; ao[g] += 1

    return {
        'n': n, 'native_bc': native_bc,
        'n_ont': n_ont, 'n_flank': n_flank, 'n_ont_end': n_ont_end,
        'true_start': ts, 'after_ont': ao, 'true_end': te, 'before_ont': bo,
        'diagnostics': {
            'chem': chem, 'native_bc': native_bc,
            'indexes': {k: detected.get(k) for k in CHEMISTRIES[chem].get('indexes', {})},
            'p5p7_anchor_frac': round(frac, 3),
            'pct_ont_start': round(100 * n_ont / n, 1) if n else 0,
            'pct_ont_end': round(100 * n_ont_end / n, 1) if n else 0,
            'warnings': warnings, 'constants': _constants(chem),
        },
    }


def to_json(result):
    return {k: (dict(v) if isinstance(v, collections.Counter) else v)
            for k, v in result.items()}


if __name__ == '__main__':
    fq = sys.argv[1]
    chem = sys.argv[2] if len(sys.argv) > 2 else '24plex'
    r = analyse(fq, chem); n = r['n']; dg = r['diagnostics']
    print(f"{os.path.basename(fq)}   reads = {n:,}   chemistry = {chem}")
    print(f"  native barcode : {r['native_bc'] or 'NONE'}")
    print(f"  i5 / i7        : {dg['indexes'].get('i5') or '-'} / {dg['indexes'].get('i7') or '-'}")
    print(f"  P5/P7 anchor   : {100*dg['p5p7_anchor_frac']:.0f}% of reads")
    print(f"  ONT start/end  : {dg['pct_ont_start']}% / {dg['pct_ont_end']}%")
    for w in dg['warnings']:
        print(f"  WARNING: {w}")
    for title, c in (('1a true FIRST', r['true_start']), ('1b AFTER ont', r['after_ont']),
                     ('2a true LAST', r['true_end']), ('2b BEFORE ont', r['before_ont'])):
        print(f'\n{title}')
        for g in ORDER:
            if c[g]:
                print(f'   {g:<34}{c[g]:>8,} ({100*c[g]/n:5.1f}%)')
    out = os.path.splitext(fq)[0] + '_composition.json'
    json.dump(to_json(r), open(out, 'w'), indent=2)
    print(f'\nJSON -> {out}')

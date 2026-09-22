#!/usr/bin/env python3
"""Barcode-accuracy filter funnel for the cross-barcode master summary.

Applies the 4-filter swap-inclusion recipe to a barcode's STAR BAM and reports
barcode accuracy + reads remaining at each of 5 cumulative steps:

  0. primary-mapped          (all primary, species-callable, barcode-bearing reads)
  1. + protein-coding mRNA   (exonic protein_coding; excludes MT/rRNA/intergenic/intron)
  2. + exact barcode + single TSO  (7-mer Hamming-dist 0; drop reads carrying both
                                    a human and a mouse barcode, i.e. TSO-TSO chimeras)
  3. + excess <=150 bp       (soft-clip excess = readlen - aligned insert)
  4. + insert 100-500 bp     (aligned insert length)

Accuracy at a step = of the barcode-bearing reads surviving it, the fraction
whose barcode-called species matches the STAR-aligned species (Wilson 95% CI).

DETECTION IS RUN-SPECIFIC: the species-barcode backbone + 7-mers vary by run
(e.g. AGTACAT/ACCTTGA on the GCAGTGGTATCAACGCAG backbone for the standard
TSO-bob runs; ATCGAAA/CAGTTGA on the R2 backbone GTGACTGGAGTTCAGACGTG for
other runs). So we read the backbone, 7-mers, and detection functions
from the run's OWN override_bob_tso_7mer module (with getattr fallbacks for the
RECOVER_/RETAIN_BACKBONE_MIN constants that older overrides lack). Built-in
AGTACAT/ACCTTGA defaults are used only if that import fails entirely.
"""
import os
import re
import math
import statistics
import bisect
import subprocess

import star_taxonomy_illumina as st
import override_24plex_illumina as _ox24   # 24plex variant detection

# --- run-specific barcode detection (sourced from the run's override) --------
_DEF_BACKBONE = 'GCAGTGGTATCAACGCAG'
_DEF_BC_HUMAN, _DEF_BC_MOUSE = 'AGTACAT', 'ACCTTGA'
_COMP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}


def _vrc(s):
    return ''.join(_COMP.get(b, 'N') for b in reversed(s))


def _ham(a, b):
    return sum(1 for x, y in zip(a, b) if x != y) if len(a) == len(b) else 99


try:
    import override_bob_tso_7mer_illumina as _ot
    BC_HUMAN = getattr(_ot, 'BC_HUMAN', _DEF_BC_HUMAN)
    BC_MOUSE = getattr(_ot, 'BC_MOUSE', _DEF_BC_MOUSE)
    rc = getattr(_ot, 'rc', _vrc)
    RECOVER_BACKBONE_MIN = getattr(_ot, 'RECOVER_BACKBONE_MIN', 2)
    RETAIN_BACKBONE_MIN = getattr(_ot, 'RETAIN_BACKBONE_MIN', 5)
    FUZZY = getattr(_ot, 'FUZZY_BARCODE_MAX_MM', 2)
    _bbmatch = getattr(_ot, '_backbone_tail_match', None)
    _cls7 = getattr(_ot, '_class_7mer', None)
    _BACKBONE = getattr(_ot, 'BACKBONE_F', _DEF_BACKBONE)
except Exception as _e:                             # pragma: no cover
    # Loud, not silent: a broken override import must not silently degrade the whole
    # detector to the built-in defaults (AGTACAT/ACCTTGA), which would make every
    # barcode call in the tso-bob funnel quietly wrong.
    import sys as _sys
    _sys.stderr.write('filter_funnel: WARNING could not import override_bob_tso_7mer_illumina ('
                      + type(_e).__name__ + ': ' + str(_e) + '); falling back to '
                      'built-in defaults ' + _DEF_BC_HUMAN + '/' + _DEF_BC_MOUSE
                      + '. tso-bob funnel numbers are NOT trustworthy for this run.\n')
    _ot = None
    BC_HUMAN, BC_MOUSE, rc = _DEF_BC_HUMAN, _DEF_BC_MOUSE, _vrc
    RECOVER_BACKBONE_MIN, RETAIN_BACKBONE_MIN, FUZZY = 2, 5, 2
    _bbmatch, _cls7, _BACKBONE = None, None, _DEF_BACKBONE

if _bbmatch is None:
    def _bbmatch(seq, end):                         # built-in fallback
        k = 0
        for j in range(1, len(_BACKBONE) + 1):
            p = end - j
            if p < 0 or seq[p] != _BACKBONE[-j]:
                break
            k += 1
        return k

if _cls7 is None:
    def _cls7(seven):                               # built-in fallback (Hamming)
        if len(seven) != 7:
            return None
        dh = _ham(seven, BC_HUMAN[:7]); dm = _ham(seven, BC_MOUSE[:7])
        if dh != dm and min(dh, dm) <= FUZZY:
            return 'h' if dh < dm else 'm'
        return None

# --- non-templated G-run after the TSO (template-switch signature) -----------
_BB_ANCHOR = _BACKBONE[-13:] if len(_BACKBONE) >= 13 else _BACKBONE
_GRUN_RE = re.compile(r'G{3,}')
# Some barcode pairs have PER-BARCODE GGG offsets (different spacer lengths, verified
# empirically): GATATGG -> +8 (2-nt spacer), TCGTGAC -> +11 (5-nt spacer). Keyed on
# the actual 7-mer (not the h/m code) so it's robust to the species labeling.
_GRUN_OFFSET_BY_7MER = {'GATATGG': 8, 'TCGTGAC': 11}

# --- generic, config-driven species-barcode caller ---------------------------
# analyze_speciesmix sets these from the run's barcode config (the bobcode/TSO
# sheet): the anchor primer immediately 5' of the barcode + the species barcodes.
# Detection then needs only those two things, so ONE caller handles any chemistry
# (24plex / tso-bob / BOB-J / ...). A barcode absent from _GRUN_OFFSET_BY_7MER ->
# GGG sits right after the barcode (offset 0 = tso-bob).
DETECT_ANCHOR = None            # primer sequence immediately 5' of the barcode
DETECT_BARCODES = {}            # {barcode_seq: 'h'|'m'}
# 'ont_anchor' | 'r2_positional' — see override_24plex_illumina.READ_LAYOUT. In
# positional mode the barcode is at read position 0 with no upstream anchor and the
# read has a fixed orientation, so anchor search and rc scanning are both skipped.
READ_LAYOUT = 'ont_anchor'


def _bcd_positional(seq):
    """ILLUMINA R2: barcode occupies seq[0:7]. Returns (code, dist) or None."""
    seven = seq[:7]
    if len(seven) != 7 or not DETECT_BARCODES:
        return None
    ranked = sorted((_ham(seven, s[:7]), c) for s, c in DETECT_BARCODES.items())
    if not ranked or ranked[0][0] > 1:
        return None
    if len(ranked) > 1 and ranked[1][0] == ranked[0][0]:
        return None                                   # equidistant -> ambiguous
    return ranked[0][1], ranked[0][0]
_DETECT_BC_MM = 1               # max Hamming mismatch for a barcode call
# MORE-SENSITIVE tso-bob detector. When DETECT_BACKBONE is set the
# barcode caller fuzzy-matches the backbone tail at every position (recovering
# ONT-degraded reads the exact-anchor find misses) + classifies the 7-mer + GGG.
# Set only for tso-bob chemistries; 24plex keeps the exact/measure_struct path.
DETECT_BACKBONE = None          # full backbone/anchor for the fuzzy tail-scan
_BB_RECOVER, _BB_RETAIN, _BB_FUZZY = 2, 5, 2

def _find_bc(seq):
    """Best barcode match right after DETECT_ANCHOR over both orientations:
    (strand, pos_after_anchor, barcode_seq, code, dist) or None."""
    best = None
    for strand in (seq, rc(seq)):
        j = strand.find(DETECT_ANCHOR)
        while j >= 0:
            p = j + len(DETECT_ANCHOR)
            for bcseq, code in DETECT_BARCODES.items():
                seg = strand[p:p + len(bcseq)]
                d = _ham(seg, bcseq) if len(seg) == len(bcseq) else 99
                if d <= _DETECT_BC_MM and (best is None or d < best[4]):
                    best = (strand, p, bcseq, code, d)
            j = strand.find(DETECT_ANCHOR, j + 1)
    return best

def _bcd_generic(seq):
    b = _find_bc(seq)
    return None if b is None else (b[3], b[4])


def _bb_tail(seq, end):
    """# of DETECT_BACKBONE chars matching seq ending at `end`, from the 3' end
    (stops at the first mismatch). The fuzzy backbone-tail match that recovers
    ONT-degraded reads."""
    bb = DETECT_BACKBONE
    k = 0
    for j in range(1, len(bb) + 1):
        p = end - j
        if p < 0 or seq[p] != bb[-j]:
            break
        k += 1
    return k


def _cls7_bb(seven):
    """Nearest species code among DETECT_BARCODES 7-mer prefixes: (code, dist) if
    unambiguous within _BB_FUZZY, else (None, 99)."""
    if len(seven) != 7:
        return None, 99
    ds = sorted((_ham(seven, s[:7]), c) for s, c in DETECT_BARCODES.items()
                if len(s) >= 7)
    if not ds:
        return None, 99
    d0, c0 = ds[0]
    d1 = ds[1][0] if len(ds) > 1 else 99
    if d0 <= _BB_FUZZY and d0 < d1:
        return c0, d0
    return None, 99


def _bcd_backbone(seq):
    """Fuzzy TSO-bob detector: at every position, fuzzy backbone-tail match
    + 7-mer species classify + GGG recovery for weak backbones. More sensitive than
    the exact-anchor find (recovers ONT-degraded reads the find() misses)."""
    best = None
    for strand in (seq, rc(seq)):
        n = len(strand)
        for i in range(2, n - 10):
            bbp = _bb_tail(strand, i)
            if bbp < _BB_RECOVER:
                continue
            code, dist = _cls7_bb(strand[i:i + 7])
            if code is None:
                continue
            ggg = strand[i + 7:i + 10].count('G') >= 2
            if bbp < _BB_RETAIN and (not ggg or dist > 1):
                continue
            sc = bbp + (100 if ggg else 0)
            if best is None or sc > best[0]:
                best = (sc, code, dist)
    return None if best is None else (best[1], best[2])

def _grun_generic(seq):
    b = _find_bc(seq)
    if b is None:
        return None
    strand, p, bcseq = b[0], b[1], b[2]
    exp = p + len(bcseq) + _GRUN_OFFSET_BY_7MER.get(bcseq, 0)
    best = 0
    for start in (exp, exp - 1, exp + 1):
        if 0 <= start < len(strand) and strand[start] == 'G':
            k = 0
            while start + k < len(strand) and strand[start + k] == 'G':
                k += 1
            best = max(best, k)
    return best


def _grun_after_bc(seq):
    """Non-templated G-run (template-switch signature): POSITIONAL count of Gs at
    the per-barcode GGG offset after the barcode. Config-driven when DETECT_* is set;
    else the legacy 24plex path."""
    if DETECT_ANCHOR and DETECT_BARCODES:
        return _grun_generic(seq)
    if _ox24 is None:
        return None
    m = _ox24.measure_struct_24plex(seq)
    if m is None or m.get('code') is None:
        return None
    strand = seq if m.get('orient') == 'F' else rc(seq)
    seven = strand[m['bc_pos']:m['bc_pos'] + 7]
    off = min(_GRUN_OFFSET_BY_7MER, key=lambda k: _ham(seven, k))   # nearest known 7-mer
    exp = m['bc_pos'] + 7 + _GRUN_OFFSET_BY_7MER[off]
    best = 0
    for start in (exp, exp - 1, exp + 1):   # anchor at the offset; +/-1 for a single ONT indel
        if 0 <= start < len(strand) and strand[start] == 'G':
            k = 0
            while start + k < len(strand) and strand[start + k] == 'G':
                k += 1
            best = max(best, k)
    return best


_CIG = re.compile(r'(\d+)([MIDN=X])')
_EXP = {'h': 'human', 'm': 'mouse'}

STEP_LABELS = [
    '0. primary-mapped',
    '1. + protein-coding mRNA',
    '2. + exact barcode + single TSO',
    '3. + insert 100-500 bp',
]


def _qaln(c):
    return sum(int(n) for n, o in _CIG.findall(c) if o in 'MI=X')


def _refend(p, c):
    return p + sum(int(n) for n, o in _CIG.findall(c) if o in 'MDN=X') - 1


def _both_bc(seq):
    """True if BOTH species 7-mers occur in the read (a TSO-TSO chimera carries two
    conflicting barcodes).

    Uses the run's configured barcodes (DETECT_BARCODES, set by analyze_speciesmix
    from the guide) and only falls back to the legacy literals 'AGTGGTG'/'AGCAGAC'
    if nothing is configured. A hardcoded pair would leave the chimera filter INERT
    on any run using a different pair (e.g. GATATGG/TCGTGAC or TGGATGA/TGTCCAC):
    the test could never fire and 0 chimeras would be reported regardless of the truth.

    ILLUMINA r2_positional: R2 has a fixed orientation, so the reverse-complement half
    of the search is dropped — it can only add false positives there."""
    codes = list(DETECT_BARCODES) or ['AGTGGTG', 'AGCAGAC']
    if len(codes) < 2:
        return False
    hay = seq if READ_LAYOUT == 'r2_positional' else (seq + '|' + rc(seq))
    return all(c in hay for c in codes)


def _bcd(seq):
    """Species-barcode call: config-driven generic caller (DETECT_*) when set,
    else legacy 24plex (override_24plex). Returns (code, dist) or None."""
    if READ_LAYOUT == 'r2_positional':
        return _bcd_positional(seq)
    if DETECT_BACKBONE and DETECT_BARCODES:
        return _bcd_backbone(seq)              # fuzzy tso-bob backbone scan
    if DETECT_ANCHOR and DETECT_BARCODES:
        return _bcd_generic(seq)
    m=_ox24.measure_struct_24plex(seq)
    return None if m is None else (m['code'], m['dist'])


def _wilson(k, n):
    if not n:
        return (None, None)
    z = 1.96
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (100 * (c - h), 100 * (c + h))


def compute_filter_funnel(bam, gtf=None, rdna=None, term=None):
    """Return the 5-step funnel dict for a STAR BAM, or None if the BAM is
    absent / unreadable. gtf/rdna/term default to the star_taxonomy caches."""
    if not bam or not os.path.exists(bam):
        return None
    if gtf is None:
        gtf = st.parse_gtf(st.GTF)
    if rdna is None:
        rdna = st.load_rdna(st.RDNA_BED)
    if term is None:
        term = st._get_terminal()
    genes, gst, exons, est = gtf
    ridna, rdst = rdna
    term_ivl, term_st = term

    def is_mrna(ch, sx, e):
        if ch.endswith('_MT'):
            return False
        v = ridna.get(ch)
        if v:
            i = bisect.bisect_right(rdst[ch], e)
            for j in range(max(0, i - 1), -1, -1):
                a, b = v[j]
                if b < sx:
                    break
                if a < e and b > sx:
                    return False
        gs = st._overlap(genes, gst, ch, sx, e, collect=True)
        if not gs or not st._overlap(exons, est, ch, sx, e):
            return False
        return any(g[2] == 'protein_coding' for g in gs)

    R = {}
    _grun = []
    try:
        pv = subprocess.Popen(['samtools', 'view', '-F', '0x904', bam],
                              stdout=subprocess.PIPE, text=True)
    except OSError:
        return None
    for line in pv.stdout:
        f = line.split('\t')
        ch = f[2]
        sp = 'human' if ch.startswith('HUMAN_') else 'mouse' if ch.startswith('MOUSE_') else None
        if sp is None:
            continue
        pos = int(f[3])
        e = _refend(pos, f[5])
        bd = _bcd(f[9])
        R[f[0]] = dict(ch=ch, sx=pos - 1, e=e, sp=sp,
                       code=(bd[0] if bd else None),
                       bsp=(_EXP[bd[0]] if bd else None),
                       dist=(bd[1] if bd else None),
                       il=_qaln(f[5]), ex=len(f[9]) - _qaln(f[5]),
                       mrna=is_mrna(ch, pos - 1, e), both=_both_bc(f[9]),
                       term=bool(st._overlap(term_ivl, term_st, ch, pos - 1, e)))
        # ILLUMINA: the G-run is NOT in the aligned SEQ — prep trims the bobcode/UMI/
        # spacer/G-run before alignment, so _grun_after_bc(f[9]) would scan cDNA and
        # produce a near-empty statistic (and "No template switch = 0.0%" when the truth
        # is ~4%). prep measures it at the exact expected position and stores it in the
        # QNAME; read it from there.
        if READ_LAYOUT == 'r2_positional':
            _mq = _ox24.measure_from_qname(f[0]) if _ox24 else None
            if _mq is not None and _mq.get('grun') is not None:
                _grun.append(_mq['grun'])
        else:
            _gl = _grun_after_bc(f[9])
            if _gl is not None:
                _grun.append(_gl)
    pv.wait()

    def acc(ids):
        cl = [r for r in ids if R[r]['code']]
        cor = sum(1 for r in cl if R[r]['bsp'] == R[r]['sp'])
        a = (100 * cor / len(cl)) if cl else None
        lo, hi = _wilson(cor, len(cl))
        return dict(n=len(ids), n_barcoded=len(cl), n_swap=len(cl) - cor,
                    acc=a, ci_lo=lo, ci_hi=hi)

    steps = []
    ids = list(R)
    steps.append(acc(ids))
    ids = [r for r in ids if R[r]['mrna']]
    steps.append(acc(ids))
    ids = [r for r in ids if R[r]['code'] and R[r]['dist'] == 0 and not R[r]['both']]
    steps.append(acc(ids))
    # No '+ excess <=150 bp' filter on Illumina — long P5/P7 adapters legitimately
    # push non-aligned length past 150 bp, so it would reject good reads.
    ids = [r for r in ids if 100 <= R[r]['il'] <= 500]
    steps.append(acc(ids))
    for i, lbl in enumerate(STEP_LABELS):
        steps[i]['label'] = lbl

    final = [r for r in ids if R[r]['code']]

    def split_correct(pred):
        cl = [r for r in final if pred(r)]
        if not cl:
            return None
        cor = sum(1 for r in cl if R[r]['bsp'] == R[r]['sp'])
        return 100 * cor / len(cl)

    return dict(
        steps=steps,
        ms_correct=split_correct(lambda r: R[r]['sp'] == 'mouse'),
        hu_correct=split_correct(lambda r: R[r]['sp'] == 'human'),
        terminal_correct=split_correct(lambda r: R[r]['term']),
        n_final=len(final),
        grun_mean=(sum(_grun) / len(_grun)) if _grun else None,
        grun_median=(statistics.median(_grun) if _grun else None),
        grun_n=len(_grun),
        grun_hist={k: _grun.count(k) for k in set(_grun)},   # {run length: read count}
    )

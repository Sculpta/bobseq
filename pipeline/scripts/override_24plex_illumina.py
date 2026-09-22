#!/usr/bin/env python3
"""24plex variant — species-barcode detection for the V1 (V1-01/V1-02) TSO design.

This chemistry uses a TSO oligo different from the AGTACAT/ACCTTGA TSO-bob
chemistry handled by override_bob_tso_7mer_illumina.py:

  oligo (as provided):  AGACGTGTGCTCTTCCGATCT[7mer]NNNNNN HH rGrGrG /3AmMO/
  as SEEN in reads:     ...CTACACGACG CTCTTCCGATCT [7mer] NNNNNN HH GGG cDNA
                                       └ R1 primer end ┘ └UMI6┘└sp┘└TS┘

  Species 7-mers:  AGTGGTG = human (V1-01),  AGCAGAC = mouse (V1-02).

NOTE (verified on 20k reads of a pilot library): the 5' handle written as
"AGACGTGTG" does NOT appear in the reads (0%); the 7-mer instead sits
immediately 3' of the standard R1 primer end "CTCTTCCGATCT" (present in 94% of
reads), and the GGG template-switch sits at a sharp +8 bp after the 7-mer
(= 6 nt UMI + 2 nt HH spacer). Detection therefore anchors on the R1 end + the
species 7-mer + the +8 GGG, which also disambiguates the barcode/TSO end from
the RT/polyT end (where polyT, not a 7-mer, follows the same R1 anchor).

UMI handling is OUT OF SCOPE here — the UMI (NNNNNN + HH) is located by this
module (umi_span) but not counted/deduplicated; that is done by
umi_dedup_illumina.py / dedup_reads_illumina.py.
"""

import re

_POLYT_RUN = re.compile(r'T{10,}')   # R1 polydT tract (Illumina RT-end evidence)

COMP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
# str.translate instead of a per-character generator + join (~111M generator steps at
# 1M reads). Every byte maps to 'N' first, then ACGTN override -- equivalent to a
# `COMP.get(b, 'N')` fallback. This copy does NOT upper() its input, unlike the
# analyze_speciesmix copy, so lowercase complements to 'N'.
_RC_TAB = {i: 78 for i in range(256)}          # 78 = ord('N')
_RC_TAB.update(str.maketrans('ACGTN', 'TGCAN'))
def rc(s): return s.translate(_RC_TAB)[::-1]

# --- 24plex chemistry constants ---------------------------------------------
R1_ANCHOR = 'CTCTTCCGATCT'      # conserved R1-primer 3' end, immediately 5' of the 7-mer
R1_FULL   = 'CTACACGACGCTCTTCCGATCT'   # full R1 end as observed in reads (for reference)
BC_HUMAN  = 'GATATGG'           # module-default 7-mer #1; 2-nt spacer -> GGG at +8
BC_MOUSE  = 'TCGTGAC'           # module-default 7-mer #2; 5-nt spacer -> GGG at +11
# The species direction is derived from the data (7-mer vs STAR-aligned species), not assumed.
EXP       = {'h': 'human', 'm': 'mouse'}
FUZZY_MM      = 1               # max Hamming mismatch for a 7-mer call
UMI_LEN       = 6               # NNNNNN
SPACER_LEN    = 2               # HH (non-G) -- ONT default only; see GRUN_OFFSET_BY_7MER
GGG_OFFSET    = UMI_LEN + SPACER_LEN     # = 8: GGG starts 8 bp after the 7-mer end
GGG_SLACK     = 2               # tolerate +/-2 bp (ONT indels) around the +8 GGG

# ---------------------------------------------------------------------------
# CONFIGURABLE 7-mers / offsets.  BC_HUMAN/BC_MOUSE above are only a DEFAULT.
# Hardcoding them is a real bug source: any run using a different bobcode pair
# would silently make every 7-mer-dependent filter inert (the chimera / ">=2 TSO"
# metric would read 0). configure() overrides them from the run's guide so nothing
# depends on the module default.
CODES = {BC_HUMAN: 'h', BC_MOUSE: 'm'}          # {7mer: 'h'|'m'}
GRUN_OFFSET_BY_7MER = {BC_HUMAN: 8, BC_MOUSE: 11}   # UMI_LEN + spacer, PER 7-mer
                                                    # (the spacer is 2 nt for one code
                                                    #  and 5 nt for the other -- a single
                                                    #  SPACER_LEN cannot describe both)

# 'ont_anchor'    : ONT/long-read. 7-mer is found by searching for R1_ANCHOR, both strands.
# 'r2_positional' : ILLUMINA R2. The TruSeq2 handle IS the sequencing primer, so it is
#                   NOT in the read: the read BEGINS with the 7-mer at position 0 and has
#                   a FIXED orientation. Searching both strands here would only add false
#                   positives, and searching for the absent anchor finds nothing at all.
READ_LAYOUT = 'ont_anchor'

# --- POLYT MEASURABILITY ----------------------------------------------------
# The R1-polydT 'RT end' is evidence of RT priming, and every check for it assumes
# R1 is long enough to CONTAIN a >=10 nt polyT tract past the 26N UMI. That is a
# property of the RUN's read geometry, not of the library: a provider may sequence
# R1 to just the UMI (a 28-cycle R1 covers only [26N UMI][TT]). Measured on such a
# run, a >=10 nt tract occurs in 0.000% of 600k reads per tube.
#
# Left unguarded this reads as TOTAL STRUCTURAL FAILURE rather than an inapplicable
# filter.
# It must NOT be auto-detected from 'no reads pass', because a genuinely failed
# library is indistinguishable from an unmeasurable one by that test; it is DECLARED
# in the run config (r1_polyt) and cross-checked against the real R1 length by
# prep_reads_illumina, which fails loudly on disagreement.
POLYT_MEASURABLE = True


def configure(species_map=None, grun_offsets=None, layout=None, polyt_measurable=None):
    """Point the detector at this run's bobcodes/offsets/geometry (call once at startup).
    species_map: {7mer: 'human'|'mouse'}; grun_offsets: {7mer: int}; layout: see above;
    polyt_measurable: False when R1 is too short to carry a >=10 nt polydT tract."""
    global CODES, GRUN_OFFSET_BY_7MER, READ_LAYOUT, POLYT_MEASURABLE
    if species_map:
        CODES = {k: ('h' if str(v).lower().startswith('hu') else 'm')
                 for k, v in species_map.items()}
    if grun_offsets:
        GRUN_OFFSET_BY_7MER = dict(grun_offsets)
    if layout:
        READ_LAYOUT = layout
    if polyt_measurable is not None:
        POLYT_MEASURABLE = bool(polyt_measurable)
    return {'CODES': CODES, 'GRUN_OFFSET_BY_7MER': GRUN_OFFSET_BY_7MER,
            'READ_LAYOUT': READ_LAYOUT, 'POLYT_MEASURABLE': POLYT_MEASURABLE}


def _offset_for(seven):
    return GRUN_OFFSET_BY_7MER.get(seven, GGG_OFFSET)


def _ham(a, b):
    return sum(1 for x, y in zip(a, b) if x != y) if len(a) == len(b) else 99


def _match_code_at(seq, p=0):
    """Nearest configured code starting at seq[p], tried longest length first, within
    FUZZY_MM, ties rejected. -> (code_seq, class, dist) or (None, None, 9). Codes are 7 or
    11 nt in the V1 set and the G-run offset counts from the END of the code (measured on
    ONT reads of the V1 set), so nothing may assume 7."""
    if not CODES:
        return None, None, 9
    for L in sorted({len(b) for b in CODES}, reverse=True):
        seg = seq[p:p + L]
        if len(seg) < L:
            continue
        ranked = sorted((_ham(seg, b), b, c) for b, c in CODES.items() if len(b) == L)
        if not ranked or ranked[0][0] > FUZZY_MM:
            continue
        if len(ranked) > 1 and ranked[1][0] == ranked[0][0]:
            return None, None, 9
        return ranked[0][1], ranked[0][2], ranked[0][0]
    return None, None, 9


def _match_7mer(seven):
    """(canonical_7mer, code, dist) for the nearest configured bobcode within FUZZY_MM,
    or (None, None, 9). Ties are rejected (ambiguous), preserving the ONT semantics."""
    if len(seven) != 7 or not CODES:
        return None, None, 9
    ranked = sorted((_ham(seven, b), b, c) for b, c in CODES.items())
    d0, b0, c0 = ranked[0]
    if d0 > FUZZY_MM:
        return None, None, 9
    if len(ranked) > 1 and ranked[1][0] == d0:      # equidistant -> ambiguous
        return None, None, 9
    return b0, c0, d0


def _class_7mer(seven):
    """'h' / 'm' / None for a 7-mer, nearest canonical within FUZZY_MM (tie-guarded)."""
    _b, code, dist = _match_7mer(seven)
    return (code, dist) if code else (None, 9)


def _ggg_confirmed(strand, seven_end):
    """True if the GGG template-switch follows the 6N UMI + spacer. The V1 set uses TWO
    spacers of different length: AGTGGTG(human) has HH -> GGG at +8; TACAAGC(mouse)
    has HHMWMWMW -> GGG at +14. The 6N UMI is always the 6 nt immediately after the
    7-mer (before the spacer), so UMI extraction is unaffected by the spacer length.
    Window covers both offsets (+8 and +14). GGG is confirmation only (the 7-mer
    already fixes the species, Hamming-7 apart)."""
    win = strand[seven_end + 6: seven_end + 18]
    return 'GGG' in win


def _ggg_confirmed_pos(seq, seven):
    """Positional G-run confirmation for Illumina R2: the non-templated Gs begin at
    7 + offset(7mer) (pos 15 for a +8 code, 18 for a +11 code). +/-1 absorbs an indel.
    Requires 2 Gs rather than 3 because the observed run length varies (2-22, mode 5)."""
    s = len(seven) + _offset_for(seven)
    return any(0 <= k and seq[k:k + 2] == 'GG' for k in (s, s - 1, s + 1))


def _measure_positional(seq, mate=None):
    """ILLUMINA R2: the read STARTS with the 7-mer. No anchor search, forward strand
    only. `mate` is R1 and supplies the RT-side (polydT / 26N UMI) evidence, which is
    physically absent from R2."""
    seven, code, dist = _match_code_at(seq, 0)
    if code is None:
        return None
    off = _offset_for(seven)
    return {
        'code': code,
        'species': EXP[code],
        'dist': dist,
        'ggg': _ggg_confirmed_pos(seq, seven),
        'orient': 'F',
        'bc_pos': 0,
        'seven': seven,
        # exact UMI window: the 6 nt immediately 3' of the 7-mer. Unlike the ONT path
        # this excludes the spacer, because the spacer length is known per 7-mer here.
        'umi_span': (len(seven), len(seven) + UMI_LEN),
        'spacer_len': max(0, off - UMI_LEN),
        # TruSeq2 is the sequencing primer -> never present in R2. Report 0 rather than
        # letting a backbone-completeness metric read as "adapter missing/degraded".
        'backbone_bp': 0,
        'polyt': _max_polyt_a(mate) if mate else 0,      # polydT lives on R1, not R2
        'rt_end': has_rt_end(seq, mate),
        'c28': has_rt_end(seq, mate),
    }


def measure_from_qname(qname):
    """ILLUMINA: rebuild the structure dict from the read NAME written by
    prep_reads_illumina.py:  <origid>_<umi1>_<umi2>_<bobcode>_<grun>_<T|F>

    Why the name and not the sequence: prep TRIMS the bobcode/UMI/spacer/G-run off R2
    before alignment, so the BAM SEQ does not contain them. Parsing the name gives the
    FASTQ path and the BAM path identical values and removes any "which sequence am I
    looking at" ambiguity. Returns None if the name is not in prep format."""
    p = qname.split()[0].split('_')
    if len(p) < 6:
        return None
    rt = p[-1]; grun = p[-2]; seven = p[-3]; umi2 = p[-4]; umi1 = p[-5]
    if seven not in CODES or rt not in ('T', 'F') or not grun.isdigit():
        return None
    off = _offset_for(seven)
    return {
        'code': CODES[seven], 'species': EXP[CODES[seven]], 'dist': 0,
        'ggg': int(grun) > 0, 'grun': int(grun), 'orient': 'F', 'bc_pos': 0,
        'seven': seven, 'umi1': umi1, 'umi2': umi2,
        'umi_span': (len(seven), len(seven) + UMI_LEN), 'spacer_len': max(0, off - UMI_LEN),
        # polyT lives on R1 and is NOT in R2: the only evidence here is prep's RT-end
        # flag (a >=10 nt T tract right after the 26N on R1). It is NOT encoded as a
        # fake length: a placeholder 30/0 would be compared to 10 by the funnel and
        # counted by the polyT histogram as a real 30-nt tract.
        'backbone_bp': 0, 'polyt': None,
        'rt_end': rt == 'T', 'c28': rt == 'T',
    }


def measure_struct_24plex(seq, mate=None):
    """Locate the species barcode on the best strand. Returns a dict
    {code, species, dist, ggg, orient, bc_pos, umi_span} or None.

    A candidate = R1 anchor + a canonical 7-mer immediately 3' of it. Candidates
    are scored by GGG-confirmation (strongly preferred) then exact-7-mer; the
    best across both orientations and all anchor hits is returned. The polyT
    RT-end (polyT after the same anchor) yields no 7-mer call and is skipped.
    """
    if READ_LAYOUT == 'r2_positional':
        return _measure_positional(seq, mate)
    best = None  # (score, code, dist, ggg, orient, p)
    for strand, ori in ((seq, 'F'), (rc(seq), 'R')):
        start = 0
        while True:
            j = strand.find(R1_ANCHOR, start)
            if j < 0:
                break
            start = j + 1
            p = j + len(R1_ANCHOR)          # 7-mer start
            cseq, code, dist = _match_code_at(strand, p)
            if code is None:
                continue
            ggg = _ggg_confirmed(strand, p + len(cseq))
            score = (10 if ggg else 0) + (1 if dist == 0 else 0)
            if best is None or score > best[0]:
                best = (score, code, dist, ggg, ori, p)
    if best is None:
        return None
    _, code, dist, ggg, ori, p = best
    strand = seq if ori == 'F' else rc(seq)
    return {
        'code': code,
        'species': EXP[code],
        'dist': dist,
        'ggg': ggg,
        'orient': ori,
        'bc_pos': p,
        'umi_span': (p + 7, p + 7 + UMI_LEN + SPACER_LEN),   # located, not counted (anchor path; 7-mer only)
        'backbone_bp': _r1_tail_match(strand, p),
        'polyt': _max_polyt_a(seq),
        'rt_end': has_rt_end(seq),
        'c28': has_rt_end(seq),     # alias: the "RT-end" element/funnel slot
    }


# --- structural helpers for the barcode_adapter interface --------------------
RETAIN_INSERT_MIN = 50
RETAIN_INSERT_MAX = 10**9          # no upper cap (mirrors TSO)
BACKBONE_F = R1_FULL               # R1 primer end immediately 5' of the barcode
BACKBONE_R = rc(R1_FULL)


def _r1_tail_match(strand, bc_pos):
    """Contiguous R1_FULL bases matched immediately 5' of the barcode position."""
    k = 0
    for j in range(1, len(R1_FULL) + 1):
        q = bc_pos - j
        if q < 0 or strand[q] != R1_FULL[-j]:
            break
        k += 1
    return k


def _max_polyt_a(seq):
    bestT = bestA = rT = rA = 0
    for c in seq:
        rT = rT + 1 if c == 'T' else 0
        rA = rA + 1 if c == 'A' else 0
        bestT = max(bestT, rT); bestA = max(bestA, rA)
    return max(bestT, bestA)


def has_rt_end(seq, mate=None):
    """RT end = R1 anchor immediately followed by a polyT run (the R1_polydT_VN
    primer). Checked in both orientations.

    ILLUMINA r2_positional: the RT end is not on R2 at all — it is the whole of R1
    ([26N UMI][polyT...]). Evidence is a polyT run starting right after the 26N UMI.
    Measured on a 24-plex library: polyT start median 26, length median 30, present in
    99.6% of R1 (vs only 36.5% of R2, where an A/T run instead means 3' polyA
    read-through — i.e. a SHORT insert, not RT priming). Returns False with no mate
    rather than guessing off R2, which would invert the metric's meaning."""
    if READ_LAYOUT == 'r2_positional':
        if not mate:
            return False
        m = _POLYT_RUN.search(mate)
        return bool(m and 18 <= m.start() <= 40)
    for strand in (seq, rc(seq)):
        st = 0
        while True:
            j = strand.find(R1_ANCHOR, st)
            if j < 0:
                break
            st = j + 1
            if strand[j + len(R1_ANCHOR): j + len(R1_ANCHOR) + 12].count('T') >= 9:
                return True
    return False


def retained_call_24plex(m, insert_len):
    """(retained, reason) — funnel: barcode -> insert>=50 -> R1-polydT RT end.
    The GGG template-switch is NOT required: the two 24plex barcodes use
    different-length post-UMI spacers (human HH -> GGG at +8; mouse HHMWMWMW -> GGG
    at +14) and the GGG is reliably visible in only a minority of mouse reads, so
    requiring it would wrongly fail structurally-complete mouse molecules. GGG is
    reported separately as an informational element (has_ggg_element / the GGG row).
    Reasons match the TSO structure-figure drop-reason buckets."""
    if m is None or m.get('code') is None:
        return False, 'no_barcode'
    if insert_len is None or insert_len < RETAIN_INSERT_MIN:
        return False, 'no_insert'
    # Only a gate when R1 can actually carry the tract. When it cannot,
    # the element is absent from the READ, not from the MOLECULE, so gating on it
    # would reject every structurally-complete molecule in the run.
    if POLYT_MEASURABLE and not m.get('rt_end'):
        return False, 'no_c28'          # reuse the RT-end bucket
    return True, 'retained'


def find_24plex_units(seq):
    """All 24plex barcode units = R1 anchor + a canonical species 7-mer, position
    sorted. Hits within 27 bp on a strand collapse to one unit. Returns
    [(pos, code 'h'/'m', orient 'F'/'R', code sequence), ...].

    ILLUMINA r2_positional: there is no anchor to search, so the ONT path returns []
    for EVERY read — and because callers gate on `len(find_24plex_units(seq)) != 1`
    that silently rejects 100% of reads as chimeras. Here the primary unit is the
    7-mer at position 0. A genuine TSO-TSO chimera shows a SECOND 7-mer further in,
    which we accept only when its own G-run sits at that 7-mer's expected offset —
    requiring the G-run keeps a chance 7-mer match from inflating the chimera rate."""
    if READ_LAYOUT == 'r2_positional':
        units = []
        seven0, code0, _d0 = _match_7mer(seq[:7])
        if code0 is None:
            return []
        units.append((0, code0, 'F'))
        for p in range(8, len(seq) - 7):
            sv, cd, _d = _match_code_at(seq, p)
            if cd is None:
                continue
            s = p + len(sv) + _offset_for(sv)
            if seq[s:s + 2] == 'GG' and p - units[-1][0] >= 27:
                units.append((p, cd, 'F'))
        return units
    out = []
    for strand, ori in ((seq, 'F'), (rc(seq), 'R')):
        last = -100
        st = 0
        while True:
            j = strand.find(R1_ANCHOR, st)
            if j < 0:
                break
            st = j + 1
            p = j + len(R1_ANCHOR)
            _cs, code, _d = _match_code_at(strand, p)
            if code is None:
                continue
            if p - last < 27:
                continue
            last = p
            out.append((p, code, ori, _cs))          # (pos, class, orientation, code sequence)
    return sorted(out)


def tso_unit_summary_24plex(seq):
    units = find_24plex_units(seq)
    return {'n_tso': len(units),
            'codes': [c for _p, c, _o, _s in units],
            'orientations': [o for _p, _c, o, _s in units],
            'code_seqs': [s for _p, _c, _o, s in units]}   # pairing by code, not species


def has_ggg_element(seq):
    """% reads with GGG switch element (GGG at +8 after a species 7-mer that
    follows the R1 anchor), either orientation."""
    for strand in (seq, rc(seq)):
        st = 0
        while True:
            j = strand.find(R1_ANCHOR, st)
            if j < 0:
                break
            st = j + 1
            p = j + len(R1_ANCHOR)
            cs, code, _d = _match_code_at(strand, p)
            if code is not None and _ggg_confirmed(strand, p + len(cs)):
                return True
    return False


def classify_read(seq):
    """Convenience: return 'human' / 'mouse' / None for a read."""
    m = measure_struct_24plex(seq)
    return m['species'] if m else None

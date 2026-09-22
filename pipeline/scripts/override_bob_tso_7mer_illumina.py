#!/usr/bin/env python3
"""Override BOB-barcode counts for the TSO-7mer chemistry.

This chemistry puts a 7-nt sample barcode in the TSO oligo (NOT in the
polydT primer as in the BOB-polydT chemistry). All 7 TSOs share the
20-nt Nextera-Read1 backbone GTGACTGGAGTTCAGACGTG followed by a sample-
specific 7-nt barcode and a rGrGrG template-switching tail:

  /5AmMC6/GTGACTGGAGTTCAGACGTG[7-nt barcode]rGrGrG    (5' DBCO-conjugated)

On the read, the backbone appears in either orientation:
  forward:  ...GTGACTGGAGTTCAGACGTG[7mer]GGG-cDNA-polyA-...
  RC:       ...polyA-cDNA-CCC[RC(7mer)]CACGTCTGAACTCCAGTCAC-...

This script scans each per-BC trimmed FASTQ, classifies every read as
'real_h' (ATCGAAA = human-side), 'real_m' (CAGTTGA = mouse-side),
'both' (both backbones detected), or 'neither' (no 7mer match within edit
distance 1 of either canonical barcode at the expected backbone-adjacent
position). Counts and species crosstabs are written back into the per-BC
{bc}_speciesmix_results.json under the keys:

  tso7_classification:       {real_h, real_m, both, neither, backbone_only}
  tso7_counts:               {ATCGAAA, CAGTTGA, both, backbone_only}
  tso7_species:              per-7mer mouse/human crosstab from the alignment pipeline
  tso7_specificity_by_insert_length: per-bin specificity (same bins as the polydT override)
  tso7_genomic_class:        terminal_exon / internal_exon / intron / intergenic counts

Run AFTER analyze_speciesmix_illumina.py (which builds {bc}_per_read.json and the
species/genomic classifications). The report sections read the overlay fields.
"""
import os
import sys
import json
import glob
import re
from collections import Counter, defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)
INDIV_DIR = os.path.join(BASE_DIR, 'individual-analyses')
FASTQ_DIR = os.path.join(BASE_DIR, 'fastq')

COMP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
def rc(s): return ''.join(COMP.get(b, 'N') for b in reversed(s))

# Nextera-R1 backbone (20 bp). Barcode follows in fwd, precedes (RC'd) in RC.
BACKBONE_F = 'GCAGTGGTATCAACGCAG'   # 3'DBCO TSO backbone (18bp as seen in reads; the
                                    # 5' AA of AAGCAGTGGTATCAACGCAG is not sequenced)
BACKBONE_R = rc(BACKBONE_F)

# Canonical 7-mer barcodes (module defaults; overridden from the run's guide).
BC_HUMAN = 'AGTACAT'  # TSO-3'rG3-A (human)
BC_MOUSE = 'ACCTTGA'  # TSO-3'rG3-B (mouse)

# RT-primer end markers (other end of the canonical molecule). A run may use
# TWO RT primers across its barcode groups:
#   * C28 RT primer  (RT_HR ATACTCGTGAC + linker ATCAGCTTCC)
#   * R1_polydT_VN   (R1_short CTACACGACGCTCTTCCGATCT + polyT + VN)
# has_c28_end() returns True if EITHER handle is present, so the "RT-end" row and
# the full-structure gate are correct for both groups. Detected as 11bp halves in
# either orientation (mismatch-tolerant against ONT error / handle truncation).
RT_HR = 'ATACTCGTGAC'        # C28 5' half
C28_LINKER = 'ATCAGCTTCC'    # C28 3' half
R1_HALF_A = 'CTACACGACGC'    # R1_short 5' half (in R1_polydT_VN)
R1_HALF_B = 'TCTTCCGATCT'    # R1_short 3' half (in R1_polydT_VN)
RT_HR_RC = rc(RT_HR)
C28_LINKER_RC = rc(C28_LINKER)
R1_HALF_A_RC = rc(R1_HALF_A)
R1_HALF_B_RC = rc(R1_HALF_B)


def has_c28_end(seq):
    """True if an RT-primer end is detectable: C28 (RT_HR or C28 linker), OR
    R1_short of R1_polydT_VN. Either orientation."""
    return (RT_HR in seq or RT_HR_RC in seq or
            C28_LINKER in seq or C28_LINKER_RC in seq or
            R1_HALF_A in seq or R1_HALF_A_RC in seq or
            R1_HALF_B in seq or R1_HALF_B_RC in seq)


def structure_category(seq, cls):
    """Bucket a read by which ends are present:
      'canonical'  = TSO end (backbone+7mer) AND C28-polydT end
      'tso_only'   = TSO end, no C28/RT end (RT end truncated or TSO-TSO chimera)
      'c28_only'   = C28/RT end, no TSO end (TSO end truncated)
      'neither'    = no detectable adapter
    cls is the classify_read() result for the read."""
    tso = (cls != 'neither')
    c28 = has_c28_end(seq)
    if tso and c28:
        return 'canonical'
    if tso and not c28:
        return 'tso_only'
    if c28 and not tso:
        return 'c28_only'
    return 'neither'


def tso_multiplicity(seq):
    """Classify how many TSO backbones a read carries and, for >=2, whether the
    7-mers encode the same or different species:
      'one'      = exactly one TSO backbone
      'two_same' = >=2 backbones, no two with positively-different species
      'two_diff' = >=2 backbones with at least two different species 7-mers
      None       = no TSO backbone
    """
    classes = []
    i = 0
    while True:
        i = seq.find(BACKBONE_F, i)
        if i < 0:
            break
        end = i + len(BACKBONE_F)
        classes.append(_class_7mer(seq[end:end + 7]) if end + 7 <= len(seq) else None)
        i += 1
    i = 0
    while True:
        i = seq.find(BACKBONE_R, i)
        if i < 0:
            break
        classes.append(_class_7mer(rc(seq[i - 7:i])) if i >= 7 else None)
        i += 1
    n = len(classes)
    if n == 0:
        return None
    if n == 1:
        return 'one'
    species = set(c for c in classes if c in ('h', 'm'))
    return 'two_diff' if len(species) >= 2 else 'two_same'

# Fuzzy 7-mer matching. BC_HUMAN (ATCGAAA) and BC_MOUSE (CAGTTGA) are Hamming-6
# apart, so a 7-mer within FUZZY_BARCODE_MAX_MM (<=2) of one canonical is
# necessarily >=(6 - MAX_MM) from the other (>=4 at MAX_MM=2) — assignment is
# unambiguous. The explicit d_H != d_M guard ("doesn't match both") only ever
# rejects the equidistant case (impossible at MAX_MM<=2; kept for future-proofing).
# At MAX_MM=0 only exact 7-mers match, recovering the strict-only mode.
FUZZY_BARCODE_MAX_MM = 2

def _hamming(a, b):
    if len(a) != len(b):
        return 10**9
    return sum(1 for x, y in zip(a, b) if x != y)

def _within_one(a, b):
    """Legacy <=1mm helper, retained for compatibility (e.g. for callers that
    want the strict 1mm bound). Most barcode-detection code paths should use
    _class_7mer, which honors FUZZY_BARCODE_MAX_MM."""
    return _hamming(a, b) <= 1

def _class_7mer(seven):
    """Return 'h' (human-side BC), 'm' (mouse-side BC), or None.
    Picks the nearest canonical within FUZZY_BARCODE_MAX_MM, with a tie-guard
    (d_h != d_m) so an equidistant 7-mer is never assigned to either."""
    if len(seven) != 7:
        return None
    d_h = _hamming(seven, BC_HUMAN)
    d_m = _hamming(seven, BC_MOUSE)
    if d_h != d_m and min(d_h, d_m) <= FUZZY_BARCODE_MAX_MM:
        return 'h' if d_h < d_m else 'm'
    return None


def classify_read(seq):
    """Scan a read for the TSO backbone in fwd or RC; if found, classify the
    adjacent 7-mer. Returns one of:
      'real_h'         — exactly one orientation found, 7mer = human-side
      'real_m'         — exactly one orientation found, 7mer = mouse-side
      'both'           — both human-side AND mouse-side 7mers found (chimeric)
      'backbone_only'  — backbone found but neither canonical 7mer (within 1 mm)
      'neither'        — no backbone match in either orientation
    """
    found_h = False
    found_m = False
    backbone_any = False
    # Forward orientation: backbone then 7mer
    i = 0
    while True:
        i = seq.find(BACKBONE_F, i)
        if i < 0:
            break
        backbone_any = True
        end = i + len(BACKBONE_F)
        if end + 7 <= len(seq):
            c = _class_7mer(seq[end:end + 7])
            if c == 'h':
                found_h = True
            elif c == 'm':
                found_m = True
        i += 1
    # RC orientation: 7mer (RC'd) before backbone (RC'd)
    i = 0
    while True:
        i = seq.find(BACKBONE_R, i)
        if i < 0:
            break
        backbone_any = True
        if i >= 7:
            c = _class_7mer(rc(seq[i - 7:i]))
            if c == 'h':
                found_h = True
            elif c == 'm':
                found_m = True
        i += 1
    if found_h and found_m:
        return 'both'
    if found_h:
        return 'real_h'
    if found_m:
        return 'real_m'
    if backbone_any:
        return 'backbone_only'
    return 'neither'


def post7mer_window(seq, up=10, dn=10):
    """Return (upstream, downstream) forward-orientation strings flanking the
    7-mer barcode, for reads with exactly ONE clean backbone+7mer match.

    Read layout (forward): [Nextera-R1 backbone][7mer]GGG-cDNA...
      upstream   = the `up` bases 5' of the 7-mer (Nextera backbone tail)
      downstream = the `dn` bases 3' of the 7-mer (the GGG template-switching
                   site followed by the start of the cDNA — this is where the
                   non-templated nucleotide addition shows up)

    Both orientations are handled by searching the read and its reverse
    complement for the forward backbone; a read appears on exactly one strand.
    Returns (None, None) if ambiguous (e.g. chimeric) or the window runs off
    the read end, so the per-position stats aren't polluted.
    """
    cands = []
    for strand in (seq, rc(seq)):
        i = 0
        while True:
            i = strand.find(BACKBONE_F, i)
            if i < 0:
                break
            end = i + len(BACKBONE_F)
            if end - up >= 0 and end + 7 + dn <= len(strand) and \
                    _class_7mer(strand[end:end + 7]) in ('h', 'm'):
                cands.append((strand[end - up:end], strand[end + 7:end + 7 + dn]))
            i += 1
    return cands[0] if len(cands) == 1 else (None, None)


# ============================================================================
# TSO-bob structural retention (QC-filter) — the TSO analog of the polydT
# bob_polydt_retention. Canonical molecule (5'->3'):
#   [Nextera-R1 backbone 20bp][7mer barcode][GGG] .. cDNA insert .. [polyA][polyT-RT C28]
# A read is QC-filtered (kept) only if it has: a canonical 7-mer barcode, the
# GGG template-switch site, an insert >=50 bp, and the C28/polyT RT end.
# Written as results['tso7_retention'] (NOT bob_polydt_retention) so the polydT
# pipeline's read-structure section never fires on TSO data, and vice-versa.
# ============================================================================
RETAIN_INSERT_MIN = 50
RETAIN_INSERT_MAX = 10**9        # no upper cap — long inserts retained
RETAIN_POLYT_MIN = 10            # >=10 bp polyT/A stretch required
RETAIN_BACKBONE_MIN = 5          # >=5 bp of Nextera-R1 adjacent to the 7-mer
BACKBONE_FULL_MIN = 18           # >=18/20 = "full" backbone (else partial)


def _backbone_tail_match(seq, end):
    """Number of Nextera-R1 backbone 3'-bases matching contiguously immediately
    5' of position `end` (i.e., how much backbone survives next to the 7-mer)."""
    k = 0
    for j in range(1, len(BACKBONE_F) + 1):
        p = end - j
        if p < 0 or seq[p] != BACKBONE_F[-j]:
            break
        k += 1
    return k


def _max_polyt_a(seq):
    """Longest run of T or A (the polyT/polyA RT-primer stretch)."""
    bestT = bestA = rT = rA = 0
    for c in seq:
        rT = rT + 1 if c == 'T' else 0
        rA = rA + 1 if c == 'A' else 0
        bestT = max(bestT, rT)
        bestA = max(bestA, rA)
    return max(bestT, bestA)


RECOVER_BACKBONE_MIN = 2   # short backbone accepted when the GGG template-
                           # switch confirms the junction — rescues 5'-truncated TSO reads


def _edist1(seg, ref):
    """Levenshtein distance between seg and ref (small strings)."""
    m, n = len(seg), len(ref)
    dp = list(range(n + 1))
    for a in range(1, m + 1):
        prev = dp[0]; dp[0] = a
        for b in range(1, n + 1):
            cur = dp[b]
            dp[b] = min(dp[b] + 1, dp[b - 1] + 1, prev + (seg[a - 1] != ref[b - 1]))
            prev = cur
    return dp[n]


def _class_barcode(strand, i):
    """Species call at position i. Substitution path first (7-mer, Hamming<=2 via
    _class_7mer); then an indel-tolerant fallback (edit-dist<=1 over a 6/8-mer
    window — the ONT deletion/insertion error mode). Returns (code|None, barcode_len,
    dist) where dist is the Hamming/edit distance of the matched window to the code."""
    code = _class_7mer(strand[i:i + 7])
    if code is not None:
        ref = BC_HUMAN if code == 'h' else BC_MOUSE
        d = sum(1 for a, b in zip(strand[i:i + 7], ref) if a != b)
        return code, 7, d
    for L in (6, 8):
        seg = strand[i:i + L]
        if len(seg) < L:
            continue
        dh = _edist1(seg, BC_HUMAN); dm = _edist1(seg, BC_MOUSE)
        if min(dh, dm) <= 1 and dh != dm:
            return ('h' if dh < dm else 'm'), L, min(dh, dm)
    return None, 7, 9


def measure_struct_tso(seq):
    """TSO-bob structural measurement (the C28-pDt-insert-TSO call). Locates the
    species barcode adjacent to the Nextera-R1/TSO backbone via two acceptance paths:
      (1) >=RETAIN_BACKBONE_MIN bp of backbone immediately 5' of the barcode (strict), OR
      (2) >=RECOVER_BACKBONE_MIN bp of backbone AND a GGG template-switch right after
          the barcode (recovery path: rescues 5'-truncated-backbone TSO reads that
          still carry [barcode]GGG; the GGG requirement keeps the short path specific).
    Barcode matching: the backbone-anchored strict path tolerates Hamming<=2, but the
    short-backbone recovery path requires a TIGHT call (exact / <=1 mismatch / single
    indel) — a 2-mismatch match without a backbone anchor lands on spurious cDNA motifs
    (~random species), so it is only trusted when >=5 bp of backbone fixes the position.
    Returns dict(code, backbone_bp, ggg, polyt, c28) or None."""
    best = None  # (score, bbp, code, i, strand, blen)
    for strand in (seq, rc(seq)):
        n = len(strand)
        for i in range(2, n - 10):
            bbp = _backbone_tail_match(strand, i)
            if bbp < RECOVER_BACKBONE_MIN:
                continue
            code, blen, dist = _class_barcode(strand, i)
            if code is None:
                continue
            ggg = strand[i + blen:i + blen + 3].count('G') >= 2
            if bbp < RETAIN_BACKBONE_MIN and (not ggg or dist > 1):
                continue                       # recovery: GGG-confirmed AND tight (<=1) barcode
            score = bbp + (100 if ggg else 0)  # prefer full-backbone / GGG-confirmed
            if best is None or score > best[0]:
                best = (score, bbp, code, i, strand, blen, dist)
    if best is None:
        return None
    score, bbp, code, i, strand, blen, dist = best
    return {'code': code, 'backbone_bp': bbp, 'dist': dist,   # dist = 7-mer Hamming (0 = exact)
            'ggg': strand[i + blen:i + blen + 3].count('G') >= 2,
            'polyt': _max_polyt_a(seq),
            'c28': has_c28_end(seq)}


def retained_call_tso(m, insert_len):
    """(retained, reason) for the loosened C28-pDt-insert-TSO call. `m` is
    measure_struct_tso() output. Funnel order: barcode -> GGG -> insert (>=50bp,
    no upper cap) -> polyT -> C28."""
    if m is None or m['code'] is None:
        return False, 'no_barcode'         # no canonical 7-mer with >=5bp backbone
    if not m['ggg']:
        return False, 'no_ggg'
    if insert_len < RETAIN_INSERT_MIN:
        return False, 'no_insert'
    # no polyT requirement for the random-hexamer+C28 RT chemistry:
    # these molecules have no polyT (no polydT priming), so the full-structure
    # keep-set is barcode -> GGG -> insert >=50 -> C28 end.
    if not m['c28']:
        return False, 'no_c28'
    return True, 'retained'


def find_tso_units(seq):
    """All TSO units in a read = a canonical 7-mer (<=1 mm) with at least
    RETAIN_BACKBONE_MIN bp of Nextera-R1 backbone immediately 5' of it. Forward
    TSOs are found in the read; reverse (rc) TSOs are found by scanning the reverse
    complement (a forward TSO's rc is 7mer-then-rc-backbone, which has no backbone
    5' of a 7-mer, so it is never double-counted). Hits within 27 bp on the same
    strand collapse to one physical TSO. Returns a position-sorted (5'->3') list of
    (pos_in_read, code 'h'/'m', orientation 'F'/'R')."""
    n = len(seq)
    out = []
    for strand, ori in ((seq, 'F'), (rc(seq), 'R')):
        L = len(strand)
        last = -100
        for i in range(7, L - 9):
            code = _class_7mer(strand[i:i + 7])
            if code is None:
                continue
            if _backbone_tail_match(strand, i) < RETAIN_BACKBONE_MIN:
                continue
            if i - last < 27:
                continue
            last = i
            out.append((i if ori == 'F' else (n - i - 7), code, ori))
    return sorted(out)


def tso_unit_summary(seq):
    """{'n_tso', 'codes', 'orientations'} for the TSO units in a read, position
    sorted (5'->3'). 1 TSO = canonical molecule; >=2 = TSO-TSO chimera. Used by
    the STAR pipeline's TSO-structure multiplicity analysis."""
    units = find_tso_units(seq)
    return {'n_tso': len(units),
            'codes': [c for _p, c, _o in units],
            'orientations': [o for _p, _c, o in units]}


# Insert-length bins (match the polydT override for plot comparability).
INSERT_LENGTH_BINS = [
    ('<300',     lambda l: l < 300),
    ('300-450',  lambda l: 300 <= l < 450),
    ('450-600',  lambda l: 450 <= l < 600),
    ('600-750',  lambda l: 600 <= l < 750),
    ('>750',     lambda l: l >= 750),
]


def compute_specificity(a_h, a_m, t_h, t_m):
    """Chance-corrected piecewise [-1, +1] specificity, exactly the same
    formula as compute_overall_specificity_from_counts in override_bob_polydt.
    Mapping: ATCGAAA (a) -> human (correct = a_h),
             CAGTTGA (t) -> mouse (correct = t_m).
    """
    n_total = a_h + a_m + t_h + t_m
    if n_total == 0:
        return None
    a_total = a_h + a_m
    t_total = t_h + t_m
    correct = a_h + t_m
    base_h = a_h + t_h
    base_pct_h = base_h / n_total
    base_pct_m = 1 - base_pct_h
    if a_total > 0 and t_total > 0:
        p_a = a_total / n_total
        p_t = t_total / n_total
        exp_correct = p_a * base_pct_h + p_t * base_pct_m
    else:
        exp_correct = base_pct_h if a_total > 0 else base_pct_m
    obs_correct = correct / n_total
    if obs_correct >= exp_correct:
        denom = 1 - exp_correct
        return (obs_correct - exp_correct) / denom if denom > 0 else 0.0
    return (obs_correct - exp_correct) / exp_correct if exp_correct > 0 else 0.0


def pct_correct(a_h, a_m, t_h, t_m):
    n = a_h + a_m + t_h + t_m
    if n == 0:
        return None
    return 100.0 * (a_h + t_m) / n


def process_barcode(bc):
    bc_dir = os.path.join(INDIV_DIR, bc)
    results_path = os.path.join(bc_dir, f'{bc}_speciesmix_results.json')
    per_read_path = os.path.join(bc_dir, f'{bc}_per_read.json')
    if not (os.path.exists(results_path) and os.path.exists(per_read_path)):
        print(f'  skip {bc}: missing JSON (run analyze_speciesmix_illumina.py first)')
        return

    fq_dir = os.path.join(FASTQ_DIR, bc)
    fq_paths = (glob.glob(os.path.join(fq_dir, 'noadapter_*.fastq'))
                or glob.glob(os.path.join(fq_dir, 'combined_*.fastq'))
                or glob.glob(os.path.join(fq_dir, '*.fastq')))
    if not fq_paths:
        print(f'  skip {bc}: no fastq')
        return
    fq = fq_paths[0]

    pr = json.load(open(per_read_path))
    rid_to_species = {r['rid']: r.get('species') for r in pr}
    rid_to_aln_len = {r['rid']: (r.get('aln_insert_len') or 0) for r in pr}
    rid_to_gc = {r['rid']: (r.get('genomic_class') or 'unmapped') for r in pr}

    # Classify every read by its TSO 7-mer. Also tally per-position base
    # composition in the 10 bp up/downstream of the 7-mer (the downstream window
    # holds the GGG template-switching site at +1..+3 and the variable
    # non-templated addition at +4/+5) for clean reads.
    NTG_FLANK = 10
    rid_to_class = {}
    cls_counts = Counter()
    ntg_up = [Counter() for _ in range(NTG_FLANK)]  # -10..-1 (backbone tail)
    ntg_dn = [Counter() for _ in range(NTG_FLANK)]  # +1..+10 (GGG + cDNA)
    ntg_n = 0
    # Read-structure category + TSO multiplicity + per-group specificity crosstab.
    struct_cats = Counter()          # canonical/tso_only/c28_only/neither
    tso_mult = Counter()             # one/two_same/two_diff
    # crosstab per structure group: a=ATCGAAA(human-side), t=CAGTTGA(mouse-side)
    struct_spec = {'canonical': Counter(), 'tso_only': Counter()}
    # --- TSO-bob structural retention (loosened C28-pDt-insert-TSO) accumulators ---
    ret_funnel = {'barcode': 0, 'barcode_ggg': 0, 'barcode_ggg_insert': 0,
                  'barcode_ggg_insert_polyt': 0, 'retained': 0}
    ret_drop = Counter()
    ret_bb_bins = {'full': 0, 'partial': 0}      # Nextera-R1 backbone next to 7-mer
    ret_ggg_bins = {'present': 0, 'absent': 0}
    ret_polyt_bins = {'present': 0, 'absent': 0}
    ret_c28_bins = {'present': 0, 'absent': 0}
    ret_spec = {lv: {'a_h': 0, 'a_m': 0, 't_h': 0, 't_m': 0}
                for lv in ('barcode', 'barcode_ggg', 'barcode_ggg_insert',
                           'barcode_ggg_insert_polyt', 'retained')}
    n_total_reads = 0
    with open(fq) as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline().rstrip('\n')
            f.readline(); f.readline()
            rid = h[1:].split()[0]
            c = classify_read(s)
            rid_to_class[rid] = c
            cls_counts[c] += 1
            # Structure category + TSO multiplicity (all reads)
            cat = structure_category(s, c)
            struct_cats[cat] += 1
            mult = tso_multiplicity(s)
            if mult:
                tso_mult[mult] += 1
            # --- TSO-bob structural retention (loosened C28-pDt-insert-TSO) ---
            # barcode found with only PARTIAL Nextera-R1 backbone (>=5bp); adds a
            # polyT/A >=10 requirement; keeps insert 50-500 + GGG + partial C28.
            n_total_reads += 1
            aln_len = rid_to_aln_len.get(rid, 0)
            mc = measure_struct_tso(s)
            retained, reason = retained_call_tso(mc, aln_len)
            ret_drop[reason] += 1
            if mc is not None and mc['code'] is not None:
                ret_funnel['barcode'] += 1
                ret_bb_bins['full' if mc['backbone_bp'] >= BACKBONE_FULL_MIN
                            else 'partial'] += 1
                ret_ggg_bins['present' if mc['ggg'] else 'absent'] += 1
                ret_polyt_bins['present' if mc['polyt'] >= RETAIN_POLYT_MIN
                               else 'absent'] += 1
                ret_c28_bins['present' if mc['c28'] else 'absent'] += 1
                _ins_ok = RETAIN_INSERT_MIN <= aln_len < RETAIN_INSERT_MAX
                _pt_ok = mc['polyt'] >= RETAIN_POLYT_MIN
                if mc['ggg']:
                    ret_funnel['barcode_ggg'] += 1
                    if _ins_ok:
                        ret_funnel['barcode_ggg_insert'] += 1
                        if _pt_ok:
                            ret_funnel['barcode_ggg_insert_polyt'] += 1
                if retained:
                    ret_funnel['retained'] += 1
                sp = rid_to_species.get(rid)
                if sp in ('mouse', 'human'):
                    key = f'{"a" if mc["code"] == "h" else "t"}_{"h" if sp == "human" else "m"}'
                    ret_spec['barcode'][key] += 1
                    if mc['ggg']:
                        ret_spec['barcode_ggg'][key] += 1
                        if _ins_ok:
                            ret_spec['barcode_ggg_insert'][key] += 1
                            if _pt_ok:
                                ret_spec['barcode_ggg_insert_polyt'][key] += 1
                    if retained:
                        ret_spec['retained'][key] += 1
            # Per-group barcode crosstab (canonical vs tso_only): single-barcode
            # reads with a species call; 'both' tallied separately as chimera.
            if cat in ('canonical', 'tso_only'):
                grp = struct_spec[cat]
                grp['n_group'] += 1
                if c == 'both':
                    grp['both'] += 1
                elif c in ('real_h', 'real_m'):
                    sp = rid_to_species.get(rid)
                    if sp in ('mouse', 'human'):
                        if c == 'real_h':
                            grp['a_h' if sp == 'human' else 'a_m'] += 1
                        else:
                            grp['t_h' if sp == 'human' else 't_m'] += 1
            if c in ('real_h', 'real_m'):
                up_s, dn_s = post7mer_window(s, NTG_FLANK, NTG_FLANK)
                if up_s and dn_s and len(up_s) == NTG_FLANK and len(dn_s) == NTG_FLANK:
                    for j, b in enumerate(up_s):
                        ntg_up[j][b] += 1
                    for j, b in enumerate(dn_s):
                        ntg_dn[j][b] += 1
                    ntg_n += 1

    # Build species × 7mer crosstabs on real_h / real_m only
    code_counts = {BC_HUMAN: 0, BC_MOUSE: 0}
    code_species = {BC_HUMAN: {'mouse': 0, 'human': 0},
                    BC_MOUSE: {'mouse': 0, 'human': 0}}
    counts_overall = {'a_h': 0, 'a_m': 0, 't_h': 0, 't_m': 0}
    counts_by_gc = defaultdict(lambda: {'a_h': 0, 'a_m': 0, 't_h': 0, 't_m': 0})
    length_bin_counts = {name: {'a_h': 0, 'a_m': 0, 't_h': 0, 't_m': 0}
                         for name, _ in INSERT_LENGTH_BINS}

    for rid, c in rid_to_class.items():
        if c == 'real_h':
            code = BC_HUMAN
            tag = 'a'
        elif c == 'real_m':
            code = BC_MOUSE
            tag = 't'
        else:
            continue
        code_counts[code] += 1
        sp = rid_to_species.get(rid)
        if sp in ('mouse', 'human'):
            code_species[code][sp] += 1
            key = f'{tag}_{"h" if sp == "human" else "m"}'
            counts_overall[key] += 1
            counts_by_gc[rid_to_gc.get(rid, 'unmapped')][key] += 1
            # Length bin
            ln = rid_to_aln_len.get(rid, 0)
            for name, pred in INSERT_LENGTH_BINS:
                if pred(ln):
                    length_bin_counts[name][key] += 1
                    break

    # Specificity per length bin
    spec_by_len = {}
    for name, _ in INSERT_LENGTH_BINS:
        bk = length_bin_counts[name]
        spec = compute_specificity(bk['a_h'], bk['a_m'], bk['t_h'], bk['t_m'])
        pc = pct_correct(bk['a_h'], bk['a_m'], bk['t_h'], bk['t_m'])
        spec_by_len[name] = {
            'atcgaaa_human': bk['a_h'],
            'atcgaaa_mouse': bk['a_m'],
            'cagttga_human': bk['t_h'],
            'cagttga_mouse': bk['t_m'],
            'n_total': bk['a_h'] + bk['a_m'] + bk['t_h'] + bk['t_m'],
            'overall_specificity': spec,
            'pct_correct': pc,
        }

    # Overall + by-genomic-class specificity
    spec_overall = compute_specificity(**counts_overall)
    pct_overall = pct_correct(**counts_overall)
    spec_by_gc = {}
    pct_by_gc = {}
    for gc_key, c in counts_by_gc.items():
        spec_by_gc[gc_key] = compute_specificity(**c)
        pct_by_gc[gc_key] = pct_correct(**c)

    # Internally-primed: all non-terminal-exon genomic classes combined
    counts_internal = {'a_h': 0, 'a_m': 0, 't_h': 0, 't_m': 0}
    for k in ('internal_exon', 'intron', 'intergenic'):
        for kk in counts_internal:
            counts_internal[kk] += counts_by_gc.get(k, {}).get(kk, 0)
    spec_internal = compute_specificity(**counts_internal)
    pct_internal = pct_correct(**counts_internal)

    # Write overlay fields into the existing results JSON.
    with open(results_path) as f:
        results = json.load(f)
    results['tso7_classification'] = dict(cls_counts)
    results['tso7_counts'] = {
        BC_HUMAN: code_counts[BC_HUMAN],
        BC_MOUSE: code_counts[BC_MOUSE],
        'both': cls_counts.get('both', 0),
        'backbone_only': cls_counts.get('backbone_only', 0),
    }
    results['tso7_species'] = code_species
    # Non-templated G addition at the GGG site (positions +1/+2/+3 after the
    # 7-mer). Read by the cross-barcode summary's G-row when present (replaces
    # the legacy TSO-seq1-variant flanking that doesn't apply to this chemistry).
    results['tso7_nontemplated_g'] = {
        'n': ntg_n,
        'upstream': [dict(c) for c in ntg_up],
        'downstream': [dict(c) for c in ntg_dn],
    }
    # Read-structure breakdown, TSO multiplicity, and per-group barcode crosstab
    # (consumed by the cross-barcode PDF "Read Structure Categories" section).
    results['tso7_structure_categories'] = {
        'canonical': struct_cats.get('canonical', 0),
        'tso_only': struct_cats.get('tso_only', 0),
        'c28_only': struct_cats.get('c28_only', 0),
        'neither': struct_cats.get('neither', 0),
        'n_total': sum(struct_cats.values()),
    }
    results['tso7_tso_multiplicity'] = {
        'one': tso_mult.get('one', 0),
        'two_same': tso_mult.get('two_same', 0),
        'two_diff': tso_mult.get('two_diff', 0),
        'n_with_tso': sum(tso_mult.values()),
    }
    results['tso7_specificity_by_structure'] = {
        'canonical': dict(struct_spec['canonical']),
        'tso_only': dict(struct_spec['tso_only']),
    }
    # TSO-bob structural retention funnel (the TSO analog of bob_polydt_retention)
    results['tso7_retention'] = {
        'n_total': n_total_reads,
        'funnel': ret_funnel,
        'backbone_bins': ret_bb_bins,        # full(>=18) vs partial(5-17) Nextera-R1
        'ggg_bins': ret_ggg_bins,
        'polyt_bins': ret_polyt_bins,
        'c28_bins': ret_c28_bins,
        'sevenmer_bins': {
            'canonical': cls_counts.get('real_h', 0) + cls_counts.get('real_m', 0),
            'off': cls_counts.get('backbone_only', 0),
            'chimera': cls_counts.get('both', 0),
            'absent': cls_counts.get('neither', 0),
        },
        'drop_reason': dict(ret_drop),
        'spec_by_level': ret_spec,
        'thresholds': {'insert_min': RETAIN_INSERT_MIN,
                       'insert_max': RETAIN_INSERT_MAX,
                       'polyt_min': RETAIN_POLYT_MIN,
                       'backbone_min': RETAIN_BACKBONE_MIN},
    }
    # No-UMI PCR-duplicate proxy for the diversity table. The structural
    # insert-extraction dedup doesn't fire for this chemistry, so estimate
    # duplicates by collapsing reads on (species | transcript | aligned insert
    # length) — long-read cDNA molecules sharing transcript AND exact aligned
    # span are very likely the same amplified molecule.
    _dupk = Counter()
    for r in pr:
        sp = r.get('species')
        if sp not in ('mouse', 'human'):
            continue
        tx = r.get('transcript')
        ln = r.get('aln_insert_len') or 0
        if not tx or ln <= 0:
            continue
        _dupk[(sp, tx, ln)] += 1
    _dup_considered = sum(_dupk.values())
    _dup_reads = sum(c for c in _dupk.values() if c > 1)
    _dup_groups = sum(1 for c in _dupk.values() if c > 1)
    results['tso7_duplicates'] = {
        'dup_reads': _dup_reads,
        'dup_groups': _dup_groups,
        'n_considered': _dup_considered,
    }
    results['tso7_overall_counts'] = counts_overall
    results['tso7_overall_specificity'] = spec_overall
    results['tso7_overall_pct_correct'] = pct_overall
    results['tso7_terminal_specificity'] = spec_by_gc.get('terminal_exon')
    results['tso7_terminal_pct_correct'] = pct_by_gc.get('terminal_exon')
    results['tso7_internal_specificity'] = spec_internal
    results['tso7_internal_pct_correct'] = pct_internal
    results['tso7_specificity_by_insert_length'] = spec_by_len
    # Genomic-class breakdown of species-called reads (used by the PDF
    # generator's terminal-exon row).
    gc_breakdown = Counter()
    for r in pr:
        if r.get('species') in ('mouse', 'human'):
            gc_breakdown[r.get('genomic_class') or 'unmapped'] += 1

    # ----- Cross-barcode-PDF compatibility aliases ----------------------
    # The existing generate_cross_barcode_summary() reads bob_polydt_* fields
    # gated by `is_bob_polydt` (= any BC has bob_polydt_classification). To
    # render the PDF with TSO 7-mer data, we write the TSO counts into the
    # same field names but keyed on the TSO codes (ATCGAAA/CAGTTGA). The
    # rendering code reads codes from a constant that analyze_speciesmix_illumina.py
    # sets to ('ATCGAAA', 'CAGTTGA') when tso7_classification is present.
    results['bob_polydt_classification'] = dict(cls_counts)
    results['bob_polydt_counts'] = {
        BC_HUMAN: code_counts[BC_HUMAN],
        BC_MOUSE: code_counts[BC_MOUSE],
        'c28_contam': 0,
        'chimeric': 0,
        'both': cls_counts.get('both', 0),
    }
    # Canonical / mRNA / mRNA+terminal-exon counts — for the TSO chemistry, we
    # don't have a stringent "canonical structure" filter. Use the species-
    # called reads with a 7-mer hit as the closest analog of all four levels.
    base_counts = {
        BC_HUMAN: code_counts[BC_HUMAN],
        BC_MOUSE: code_counts[BC_MOUSE],
        'total': code_counts[BC_HUMAN] + code_counts[BC_MOUSE],
    }
    base_species = code_species
    results['bob_polydt_canonical_counts'] = base_counts
    results['bob_polydt_canonical_species'] = base_species
    results['bob_polydt_canonical_mrna_counts'] = base_counts
    results['bob_polydt_canonical_mrna_species'] = base_species
    results['bob_polydt_mrna_counts'] = base_counts
    results['bob_polydt_mrna_species'] = base_species
    # ----- Main 3C table + master-summary compatibility ------------------
    # The canonical generate_cross_barcode_summary 3C tables and the master-
    # summary "% barcode correct" / "Overall (specificity)" rows read the core
    # crosstab field `bob_barcode_species` (+ `bob_barcode_counts`). The BOB
    # 4-mer extractor leaves these EMPTY for 7-mer TSO chemistry, so populate
    # them here from the all-reads 7-mer crosstab so those rows render for TSO
    # runs (computed AFTER add_biotype_species, so this is the final value).
    results['bob_barcode_species'] = code_species
    results['bob_barcode_counts'] = {
        BC_HUMAN: code_counts[BC_HUMAN],
        BC_MOUSE: code_counts[BC_MOUSE],
    }
    # mRNA + terminal-exon — slice by terminal-exon genomic_class on real 7-mer reads
    term_counts_tso = {BC_HUMAN: 0, BC_MOUSE: 0}
    term_species_tso = {BC_HUMAN: {'mouse': 0, 'human': 0},
                        BC_MOUSE: {'mouse': 0, 'human': 0}}
    for rid, c in rid_to_class.items():
        code = BC_HUMAN if c == 'real_h' else BC_MOUSE if c == 'real_m' else None
        if code is None:
            continue
        sp = rid_to_species.get(rid)
        if sp not in ('mouse', 'human'):
            continue
        if rid_to_gc.get(rid) == 'terminal_exon':
            term_counts_tso[code] += 1
            term_species_tso[code][sp] += 1
    results['bob_polydt_mrna_terminal_counts'] = {
        BC_HUMAN: term_counts_tso[BC_HUMAN],
        BC_MOUSE: term_counts_tso[BC_MOUSE],
        'total': term_counts_tso[BC_HUMAN] + term_counts_tso[BC_MOUSE],
    }
    results['bob_polydt_mrna_terminal_species'] = term_species_tso
    # Specificity-by-insert-length — copy with renamed inner keys so the
    # rendering code's lookup of atcgct_human / atcgct_mouse / tcgata_human
    # / tcgata_mouse still works (it expects those exact keys).
    spec_by_len_alias = {}
    for name, _ in INSERT_LENGTH_BINS:
        bk = length_bin_counts[name]
        spec_by_len_alias[name] = {
            'atcgct_human': bk['a_h'],
            'atcgct_mouse': bk['a_m'],
            'tcgata_human': bk['t_h'],
            'tcgata_mouse': bk['t_m'],
            'n_total': bk['a_h'] + bk['a_m'] + bk['t_h'] + bk['t_m'],
            'overall_specificity': compute_specificity(
                bk['a_h'], bk['a_m'], bk['t_h'], bk['t_m']),
        }
    results['bob_polydt_specificity_by_insert_length'] = spec_by_len_alias
    # Genomic-class breakdown
    n_gc_total = sum(gc_breakdown.values())
    results['bob_polydt_mrna_genomic_class_breakdown'] = {
        'n_total': n_gc_total,
        'terminal_exon': gc_breakdown.get('terminal_exon', 0),
        'internal_exon': gc_breakdown.get('internal_exon', 0),
        'intron': gc_breakdown.get('intron', 0),
        'intergenic': gc_breakdown.get('intergenic', 0),
        'unmapped': gc_breakdown.get('unmapped', 0),
    }
    results['tso7_mrna_genomic_class_breakdown'] = dict(gc_breakdown)
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)

    n_h = code_counts[BC_HUMAN]
    n_m = code_counts[BC_MOUSE]
    n_real = n_h + n_m
    pct_h = 100 * n_h / n_real if n_real else 0
    n_back_only = cls_counts.get('backbone_only', 0)
    n_neither = cls_counts.get('neither', 0)
    spec_str = f'{spec_overall:+.3f}' if spec_overall is not None else 'NA'
    pct_str = f'{pct_overall:.1f}%' if pct_overall is not None else 'NA'
    print(f'  {bc}: {n_real:>6,} TSO-7mer ({n_h:,} H + {n_m:,} M, '
          f'{pct_h:.0f}%H)  backbone_only={n_back_only:,}  no_TSO={n_neither:,}  '
          f'spec={spec_str}  %correct={pct_str}')


def main():
    bcs = sorted(d for d in os.listdir(INDIV_DIR)
                 if d.startswith('barcode') and os.path.isdir(os.path.join(INDIV_DIR, d)))
    for bc in bcs:
        process_barcode(bc)


if __name__ == '__main__':
    main()

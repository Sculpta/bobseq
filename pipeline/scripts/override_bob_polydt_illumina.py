#!/usr/bin/env python3
"""Override BOB-barcode counts in per-barcode JSONs to use 6mer BOB-polydT
codes (ATCGCT / TCGATA) instead of the standard BOB-C5 4mer codes.

Why: analyze_speciesmix_illumina.py's BOB detection looks for ATGTC[4bp]ACATG context,
which doesn't match the 6mer BOB-polydT chemistry (ATGTC[6bp]TTTTT, no
BOB_HR suffix). The downstream report fields (bob_barcode_counts,
bob_barcode_species) end up empty for this experiment.

This script: scans the trimmed FASTQ per barcode, detects bob-pdt-1 (ATCGCT)
and bob-pdt-2 (TCGATA) in their expected polyT-flanked context (fwd + RC),
classifies each read against the C28-RT contamination signature
(ATACTCGTGACATCAGCTTCC), and ALSO checks for the full canonical read
structure (TSO + insert + polydT + 6mer + bob-c in the right orientation
and order). Pulls the species call from the per_read JSON, then writes the
6mer counts back to bob_barcode_counts / bob_barcode_species and the
canonical-only counts to bob_polydt_canonical_* on each
{bc}_speciesmix_results.json.

Run *after* analyze_speciesmix_illumina.py has produced the per-read JSON and the
species / biotype classifications.

Diagnostic fields added:
  bob_polydt_classification:           {real_p1, real_p2, chimeric, c28_only, neither, both}
  bob_polydt_counts:                   {ATCGCT, TCGATA, c28_contam, chimeric, both}
  bob_polydt_canonical_counts:         {ATCGCT, TCGATA, total}    # full structure only
  bob_polydt_canonical_species:        per-6mer mouse/human crosstab on canonical subset
  bob_polydt_canonical_mrna_counts:    {ATCGCT, TCGATA, total}    # canonical + mRNA only
                                       (excludes MT genes + rRNA + ribosomal proteins)
  bob_polydt_canonical_mrna_species:   per-6mer mouse/human crosstab on canonical+mRNA subset
"""
import os
import sys
import json
import glob
import re
from collections import Counter

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)
INDIV_DIR = os.path.join(BASE_DIR, 'individual-analyses')
FASTQ_DIR = os.path.join(BASE_DIR, 'fastq')

COMP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
def rc(s): return ''.join(COMP.get(b, 'N') for b in reversed(s))

# 6mer detection ANCHOR (11bp exact): ATGTC + 6mer  (last 5bp of bob-c + the 6mer)
# We then verify a polyT-like region (≥3 T's in next 5 bases, fuzzy-tolerant
# for ONT homopolymer errors) — rather than a strict TTTTT exact match.
P1_ANCHOR = 'ATGTC' + 'ATCGCT'      # 11bp fwd anchor
P2_ANCHOR = 'ATGTC' + 'TCGATA'
P1_ANCHOR_RC = rc(P1_ANCHOR)        # 11bp rc anchor (preceded by polyA in RC)
P2_ANCHOR_RC = rc(P2_ANCHOR)

# Fuzzy polyT requirement: in the 5 bases adjacent to the anchor, require at
# least MIN_T_IN_WINDOW T's (or A's on the RC side). A strict 5/5 would allow no
# errors; this allows up to 2 substitution errors per 5 bases.
POLYT_WINDOW = 5
MIN_T_IN_WINDOW = 3

# C28-RT contamination signature: RT_HR + C28 linker (21bp)
C28_FWD = 'ATACTCGTGACATCAGCTTCC'
C28_RC = rc(C28_FWD)
# Full bob-c primer (20bp) — required for canonical structure
BOBC_FWD = 'ATGTAGTCCGTACGCATGTC'
BOBC_RC = rc(BOBC_FWD)
# TSO primer (18bp) — diagnostics only; not required by the
# canonical filter (read truncation drops TSO too often).
TSO_FWD = 'GCAGTGGTATCAACGCAG'
TSO_RC = rc(TSO_FWD)


def has_anchor_with_fuzzy_polyT(seq, anchor_fwd, anchor_rc):
    """True if `seq` contains the 11bp ATGTC[6mer] anchor (fwd or RC) followed
    (or preceded, in RC) by a polyT-like region with at least MIN_T_IN_WINDOW
    T's (or A's, in RC) within POLYT_WINDOW adjacent bases. Tolerates up to
    POLYT_WINDOW - MIN_T_IN_WINDOW substitution errors in the immediate polyT
    flank — addresses the ONT homopolymer error mode, where a strict
    TTTTT exact match would drop many real BOB-polydT reads.
    """
    # Forward orientation: anchor + polyT downstream
    i = 0
    while True:
        i = seq.find(anchor_fwd, i)
        if i < 0:
            break
        end = i + len(anchor_fwd)
        window = seq[end:end + POLYT_WINDOW]
        if window.count('T') >= MIN_T_IN_WINDOW:
            return True
        i += 1
    # Reverse-complement orientation: anchor preceded by polyA upstream
    i = 0
    while True:
        i = seq.find(anchor_rc, i)
        if i < 0:
            break
        window = seq[max(0, i - POLYT_WINDOW):i]
        if window.count('A') >= MIN_T_IN_WINDOW:
            return True
        i += 1
    return False


def has_full_bobc_primer(seq):
    """True if the full 20bp bob-c primer is present in either orientation."""
    return (BOBC_FWD in seq) or (BOBC_RC in seq)


# ---------------------------------------------------------------------------
# Structure API consumed by barcode_adapter._polydt_adapter (P1_6MER / measure_structure /
# retained_call / RETAIN_*), i.e. what selecting barcode_location='BOBcode' (the
# DEFAULT) needs. Definitions follow the canonical structure stated
# above: bob-c primer (full 20 bp) + 6mer + fuzzy polyT (>= MIN_T_IN_WINDOW T in
# POLYT_WINDOW) + an aligned insert. No TSO or C28 requirement.
# ---------------------------------------------------------------------------
P1_6MER = P1_ANCHOR[5:]          # 'ATCGCT'
P2_6MER = P2_ANCHOR[5:]          # 'TCGATA'
RETAIN_POLYT_MIN = MIN_T_IN_WINDOW   # T's in the POLYT_WINDOW right after the 6mer
RETAIN_INSERT_MIN = 20               # aligned insert (bp); reads below are primer/adapter only
RETAIN_INSERT_MAX = 10**9            # no upper cap


def _polyt_window_count(seq, anchor_fwd, anchor_rc):
    """Best polyT evidence for one anchor: max T count in the POLYT_WINDOW after a
    forward anchor hit, or A count before an RC anchor hit. -1 = anchor absent."""
    best = -1
    i = 0
    while True:
        i = seq.find(anchor_fwd, i)
        if i < 0:
            break
        end = i + len(anchor_fwd)
        best = max(best, seq[end:end + POLYT_WINDOW].count('T'))
        i += 1
    i = 0
    while True:
        i = seq.find(anchor_rc, i)
        if i < 0:
            break
        best = max(best, seq[max(0, i - POLYT_WINDOW):i].count('A'))
        i += 1
    return best


def measure_structure(seq, mate=None):
    """Per-read structure for the BOB-polydT chemistry.
    Returns None when neither 6mer anchor is present; else a dict with
      code      'p1' | 'p2' | 'both' ('both' = two species anchors = chimera; it is
                deliberately NOT a key of the adapter's EXP so it is never scored)
      polyt     T count in the POLYT_WINDOW after the anchor (0..POLYT_WINDOW)
      bobc_full full 20 bp bob-c primer present (either orientation)
      c28, tso  presence of the C28 / TSO sequences (informational only)
      bc_len    6
    The insert is NOT in this dict: for polydT it is the STAR-aligned length, which
    the adapter receives separately as `aln_len`."""
    t1 = _polyt_window_count(seq, P1_ANCHOR, P1_ANCHOR_RC)
    t2 = _polyt_window_count(seq, P2_ANCHOR, P2_ANCHOR_RC)
    if t1 < 0 and t2 < 0:
        return None
    if t1 >= 0 and t2 >= 0:
        code, polyt = 'both', max(t1, t2)
    elif t1 >= 0:
        code, polyt = 'p1', t1
    else:
        code, polyt = 'p2', t2
    return {'code': code, 'polyt': polyt, 'bobc_full': has_full_bobc_primer(seq),
            'c28': (C28_FWD in seq) or (C28_RC in seq),
            'tso': (TSO_FWD in seq) or (TSO_RC in seq), 'bc_len': 6}


def retained_call(m, aln_len=0):
    """(kept, reason) for a measure_structure() result. Canonical = single species
    anchor + full bob-c primer + fuzzy polyT + aligned insert >= RETAIN_INSERT_MIN."""
    if not m or m.get('code') not in ('p1', 'p2'):
        return False, 'no_code'
    if not m.get('bobc_full'):
        return False, 'no_bobc_primer'
    if m.get('polyt', 0) < RETAIN_POLYT_MIN:
        return False, 'polyt_short'
    if not (RETAIN_INSERT_MIN <= (aln_len or 0) < RETAIN_INSERT_MAX):
        return False, 'insert_short'
    return True, 'retained'


# Insert-length bins for the per-BC specificity-vs-length table.
# First bin is everything <300 (<150 and 150-300 are merged — short reads are too
# sparse to be useful on their own; the merged <300 bin still resolves
# short-read effects while doubling the per-bin n).
INSERT_LENGTH_BINS = [
    ('<300',     lambda l: l < 300),
    ('300-450',  lambda l: 300 <= l < 450),
    ('450-600',  lambda l: 450 <= l < 600),
    ('600-750',  lambda l: 600 <= l < 750),
    ('>750',     lambda l: l >= 750),
]


def compute_overall_specificity_from_counts(a_h, a_m, t_h, t_m, swapped=False):
    """Pooled [-1, +1] specificity score using the same formula as the cross-
    barcode PDF. Returns None if no reads.

    Default expected mapping:    ATCGCT->human, TCGATA->mouse  (correct = a_h + t_m)
    Swapped (per-BC override):    ATCGCT->mouse, TCGATA->human  (correct = a_m + t_h)

    The "base composition" terms (base_pct_h / base_pct_m) describe what
    fraction of the read pool actually is human/mouse — these are independent
    of which 6mer SHOULD have mapped where. Only the `correct` count and
    `exp_correct` flip.
    """
    n_total = a_h + a_m + t_h + t_m
    if n_total == 0:
        return None
    a_total = a_h + a_m
    t_total = t_h + t_m
    if swapped:
        correct = a_m + t_h  # ATCGCT->mouse + TCGATA->human is correct
    else:
        correct = a_h + t_m  # ATCGCT->human + TCGATA->mouse is correct
    base_h = a_h + t_h
    base_pct_h = base_h / n_total
    base_pct_m = 1 - base_pct_h
    if a_total > 0 and t_total > 0:
        p_a = a_total / n_total
        p_t = t_total / n_total
        if swapped:
            # ATCGCT (a) expected to be mouse, TCGATA (t) expected to be human
            exp_correct = p_a * base_pct_m + p_t * base_pct_h
        else:
            exp_correct = p_a * base_pct_h + p_t * base_pct_m
    else:
        if swapped:
            exp_correct = base_pct_m if a_total > 0 else base_pct_h
        else:
            exp_correct = base_pct_h if a_total > 0 else base_pct_m
    obs_correct = correct / n_total
    if obs_correct >= exp_correct:
        denom = 1 - exp_correct
        return (obs_correct - exp_correct) / denom if denom > 0 else 0.0
    denom = exp_correct
    return (obs_correct - exp_correct) / denom if denom > 0 else 0.0


def parse_gene_from_transcript(transcript):
    """Return (species_prefix, gene_name) or (None, None) if not parseable."""
    if not transcript:
        return None, None
    parts = transcript.split('|')
    if len(parts) < 8:
        return None, None
    if parts[0].startswith('MOUSE_'):
        sp = 'mouse'
    elif parts[0].startswith('HUMAN_'):
        sp = 'human'
    else:
        return None, None
    return sp, parts[5]

# mRNA filter — define what counts as a "real mRNA" read.
# Required:  transcript biotype = protein_coding
# Excluded:  mitochondrial genes (gene name starts MT-/mt-/Mt-)
#            ribosomal protein subunits (cytoplasmic RPL/RPS/RPLP and
#            mitochondrial MRPL/MRPS)
# Rationale: bulk RNA-seq is dominated by rRNA + mt-rRNA + ribosomal protein
# mRNAs which behave as repetitive / over-abundant. Excluding them gives a
# cleaner readout of regulated mRNA expression.
MT_RE = re.compile(r'^[Mm][Tt]-')
# Match ribosomal protein subunits — cytoplasmic (RPL/RPS/RPLP) and
# mitochondrial (MRPL/MRPS) — including paralog variants like RPL26L1,
# RPS27L, RPL36AL.
#
# Approach: require the gene name to start with the (M?)RPL or (M?)RPS prefix
# and then contain only digits/letters with NO 'K' character. The 'K' rule
# excludes kinases (RPS6KA1, RPS6KB1, RPS6KC1) which are NOT ribosomal
# subunits — they're S6 kinases and should NOT be filtered out.
#
# Pattern is case-insensitive (mouse symbols use mixed case: Rpl3, Rps17,
# Rplp1, etc.)
RIBO_RE = re.compile(r'^[Mm]?[Rr][Pp][LlSs][^Kk]*$')

def is_mrna_read(transcript):
    """Given the per-read transcript field (pipe-delimited), decide whether
    the read aligns to a 'real' mRNA (protein_coding, non-MT, non-ribosomal-protein)."""
    if not transcript:
        return False
    parts = transcript.split('|')
    if len(parts) < 8:
        return False
    gene_name = parts[5]
    biotype = parts[7]
    if biotype != 'protein_coding':
        return False
    if MT_RE.match(gene_name):
        return False
    if RIBO_RE.match(gene_name):
        return False
    return True


def classify_read(seq):
    has_p1 = has_anchor_with_fuzzy_polyT(seq, P1_ANCHOR, P1_ANCHOR_RC)
    has_p2 = has_anchor_with_fuzzy_polyT(seq, P2_ANCHOR, P2_ANCHOR_RC)
    has_c28 = (C28_FWD in seq) or (C28_RC in seq)
    if has_p1 and has_p2:
        return 'both'
    if has_p1 and has_c28:
        return 'chimeric_p1'
    if has_p2 and has_c28:
        return 'chimeric_p2'
    if has_p1:
        return 'real_p1'
    if has_p2:
        return 'real_p2'
    if has_c28:
        return 'c28_only'
    return 'neither'


# Canonical structure for BOB-polydT chemistry.
# The TSO is not required: ONT read truncation drops the TSO end on roughly 60%
# of otherwise-good BOB-polydT reads. The practical canonical structure is:
#
#   bob-c primer (full 20bp) + 6mer barcode + polyT (with errors allowed) + insert
#
# That is, *no* TSO requirement. The polyT uses fuzzy detection (≥3 of
# next 5 bases are T) — see has_anchor_with_fuzzy_polyT above.
# "Insert" is implied by the read having a mouse/human alignment (species
# call), since unmapped reads have no insert call.
#
# Implementation: a read is canonical iff
#   1) classify_read returned 'real_p1' or 'real_p2'  (6mer + fuzzy polyT)
#   2) the full 20bp bob-c primer is present in the sequence (fwd or RC)
#   3) the read has a mouse/human species call from the alignment pipeline


# Per-BC species-map overrides — parsed once from the analysisguide so the
# specificity formula uses the right "expected mapping" for each BC. The
# convention: if the BC appears here, ATCGCT->mouse and TCGATA->human is the
# correct mapping (rather than the default ATCGCT->human / TCGATA->mouse).
def _load_swapped_bc_set():
    swapped = set()
    mrg = os.path.join(BASE_DIR, 'machine_readable_guide')
    if not os.path.isdir(mrg):
        return swapped
    for fname in os.listdir(mrg):
        if not (fname.endswith('_analysisguide.txt') or fname.endswith('_guide.txt')):
            continue
        try:
            with open(os.path.join(mrg, fname)) as f:
                in_cfg = False
                for line in f:
                    s = line.strip()
                    if s == 'BEGIN_CONFIG':
                        in_cfg = True
                        continue
                    if s == 'END_CONFIG':
                        in_cfg = False
                        continue
                    if not in_cfg or '\t' not in s:
                        continue
                    k, _, v = s.partition('\t')
                    if k.strip() != 'barcode_species_overrides':
                        continue
                    for entry in v.split(','):
                        entry = entry.strip()
                        if ':' not in entry:
                            continue
                        bc_name, pairs = entry.split(':', 1)
                        ov = {}
                        for pair in pairs.split(';'):
                            if '=' in pair:
                                code, sp = pair.strip().split('=', 1)
                                ov[code.strip()] = sp.strip()
                        if ov.get('ATCGCT') == 'mouse' and ov.get('TCGATA') == 'human':
                            swapped.add(bc_name.strip())
        except Exception:
            pass
    return swapped


SWAPPED_BCS = _load_swapped_bc_set()
if SWAPPED_BCS:
    print(f'  per-BC species-map overrides (ATCGCT=mouse,TCGATA=human): '
          f'{", ".join(sorted(SWAPPED_BCS))}')


def process_barcode(bc):
    bc_dir = os.path.join(INDIV_DIR, bc)
    results_path = os.path.join(bc_dir, f'{bc}_speciesmix_results.json')
    per_read_path = os.path.join(bc_dir, f'{bc}_per_read.json')
    if not os.path.exists(results_path) or not os.path.exists(per_read_path):
        print(f'  skip {bc}: missing JSON')
        return

    fq_dir = os.path.join(FASTQ_DIR, bc)
    fq_paths = glob.glob(os.path.join(fq_dir, 'noadapter_*.fastq')) \
        or glob.glob(os.path.join(fq_dir, 'combined_*.fastq')) \
        or glob.glob(os.path.join(fq_dir, '*.fastq'))
    if not fq_paths:
        print(f'  skip {bc}: no fastq')
        return
    fq = fq_paths[0]

    pr = json.load(open(per_read_path))
    rid_to_species = {r['rid']: r['species'] for r in pr if r.get('species') in ('mouse', 'human')}
    # Per-rid mRNA gate (canonical + mRNA filter): True iff this read aligns
    # to a protein_coding transcript that's not MT and not ribosomal-protein.
    rid_is_mrna = {r['rid']: is_mrna_read(r.get('transcript')) for r in pr}
    # Per-rid terminal-exon gate: pipeline classifies each mapped read's
    # genomic_class as one of {terminal_exon, internal_exon, intron,
    # intergenic, unmapped}. Used by the "mRNA + terminal exon" filter.
    rid_is_terminal = {r['rid']: r.get('genomic_class') == 'terminal_exon' for r in pr}
    # Per-rid insert length (alignment) and gene name (for top-mRNA tables).
    rid_insert_len = {r['rid']: (r.get('aln_insert_len') or 0) for r in pr}
    rid_gene = {}  # rid -> gene_name (or None)
    for r in pr:
        sp_pref, gname = parse_gene_from_transcript(r.get('transcript'))
        rid_gene[r['rid']] = gname
    # Per-rid full-bob-c gate is filled below during the FASTQ scan.
    rid_has_full_bobc = {}

    rid_to_class = {}
    classification_counts = {'real_p1': 0, 'real_p2': 0, 'both': 0,
                             'chimeric_p1': 0, 'chimeric_p2': 0,
                             'c28_only': 0, 'neither': 0}
    with open(fq) as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline().rstrip('\n')
            f.readline()
            f.readline()
            rid = h[1:].split()[0]
            cls = classify_read(s)
            rid_to_class[rid] = cls
            classification_counts[cls] = classification_counts.get(cls, 0) + 1
            # Cache full bob-c presence for the canonical filter — saves us a
            # second pass over the FASTQ later.
            rid_has_full_bobc[rid] = has_full_bobc_primer(s)

    # Build 6mer counts and species crosstab from real_p1 / real_p2 reads only.
    bob_counts = {'ATCGCT': 0, 'TCGATA': 0}
    bob_species = {'ATCGCT': {'mouse': 0, 'human': 0},
                   'TCGATA': {'mouse': 0, 'human': 0}}
    canon_counts = {'ATCGCT': 0, 'TCGATA': 0}
    canon_species = {'ATCGCT': {'mouse': 0, 'human': 0},
                     'TCGATA': {'mouse': 0, 'human': 0}}
    canon_mrna_counts = {'ATCGCT': 0, 'TCGATA': 0}
    canon_mrna_species = {'ATCGCT': {'mouse': 0, 'human': 0},
                          'TCGATA': {'mouse': 0, 'human': 0}}
    # mRNA-only filter (no canonical-structure requirement). Allows
    # truncated reads that have a 6mer + an mRNA alignment, even if the
    # other adapter components are missing.
    mrna_counts = {'ATCGCT': 0, 'TCGATA': 0}
    mrna_species = {'ATCGCT': {'mouse': 0, 'human': 0},
                    'TCGATA': {'mouse': 0, 'human': 0}}
    # mRNA + terminal-exon filter (also no canonical requirement). The
    # terminal-exon genomic_class is the cleanest 3'-end signal that the
    # read came from a properly polyadenylated mature mRNA.
    mrna_term_counts = {'ATCGCT': 0, 'TCGATA': 0}
    mrna_term_species = {'ATCGCT': {'mouse': 0, 'human': 0},
                         'TCGATA': {'mouse': 0, 'human': 0}}
    # Per-bin insert-length specificity tracker.
    # For each insert-length bin, accumulate (a_h, a_m, t_h, t_m) counts on
    # canonical-filtered reads. Specificity per bin computed at the end.
    length_bin_counts = {name: {'a_h': 0, 'a_m': 0, 't_h': 0, 't_m': 0}
                         for name, _ in INSERT_LENGTH_BINS}
    # Per-species top-10 mRNA tables, computed on ANY mRNA-aligned read
    # (NOT gated by canonical structure / 6mer presence). This makes the
    # tables populated even for C28 single-species controls and other BCs
    # where BOB-polydT chemistry is absent.
    mrna_genes_human = Counter()
    mrna_genes_mouse = Counter()
    # Genomic-class breakdown for mRNA-aligned reads — diagnostic for
    # whether polyT priming landed at the polyA tail (terminal_exon) vs
    # internal sites (internal_exon, intron) or alignment artifacts
    # (intergenic, unmapped).
    mrna_genomic_class = Counter()
    for rid, cls in rid_to_class.items():
        if cls == 'real_p1':
            code = 'ATCGCT'
        elif cls == 'real_p2':
            code = 'TCGATA'
        else:
            continue
        bob_counts[code] += 1
        sp = rid_to_species.get(rid)
        if sp in ('mouse', 'human'):
            bob_species[code][sp] += 1
        # Canonical filter: real_p1/real_p2 + full 20bp bob-c
        # primer + species call (= insert aligned). NO TSO requirement.
        if rid_has_full_bobc.get(rid) and sp in ('mouse', 'human'):
            canon_counts[code] += 1
            canon_species[code][sp] += 1
            # Tally insert-length bin (canonical reads only — these are the
            # ones with both species and structure information).
            ln = rid_insert_len.get(rid, 0)
            for name, pred in INSERT_LENGTH_BINS:
                if pred(ln):
                    bk = length_bin_counts[name]
                    if code == 'ATCGCT':
                        bk['a_h' if sp == 'human' else 'a_m'] += 1
                    else:
                        bk['t_h' if sp == 'human' else 't_m'] += 1
                    break
            # Canonical + mRNA-only filter (exclude MT + ribosomal-protein +
            # any non-protein_coding biotype).
            if rid_is_mrna.get(rid):
                canon_mrna_counts[code] += 1
                canon_mrna_species[code][sp] += 1
        # mRNA-only filter (no canonical-structure requirement)
        if rid_is_mrna.get(rid):
            mrna_counts[code] += 1
            if sp in ('mouse', 'human'):
                mrna_species[code][sp] += 1
            # mRNA + terminal-exon filter (still no canonical requirement)
            if rid_is_terminal.get(rid):
                mrna_term_counts[code] += 1
                if sp in ('mouse', 'human'):
                    mrna_term_species[code][sp] += 1

    # Independent (non-canonical) tallies — computed once per per-read entry
    # so they work for ANY chemistry, including C28-only controls that don't
    # have BOB-polydT reads at all.
    for r in pr:
        rid = r['rid']
        if not rid_is_mrna.get(rid):
            continue
        sp = rid_to_species.get(rid)
        if sp not in ('mouse', 'human'):
            continue
        # Per-species top-10 mRNA tally (any read that aligns to a real mRNA)
        gname = rid_gene.get(rid)
        if gname:
            if sp == 'human':
                mrna_genes_human[gname] += 1
            else:
                mrna_genes_mouse[gname] += 1
        # Genomic-class breakdown of mRNA-aligned reads
        gc = r.get('genomic_class') or 'unmapped'
        mrna_genomic_class[gc] += 1

    # Update the per-barcode results JSON
    with open(results_path) as f:
        results = json.load(f)
    results['bob_barcode_counts'] = bob_counts
    results['bob_barcode_species'] = bob_species
    # Clear BOB-C5-style canonical/dedup overrides — those tables expect
    # ATGTC[4mer]ACATG context with BOB_HR suffix that doesn't apply here.
    # We add our own BOB-polydT-specific canonical fields below.
    results['bob_species_clean'] = {}
    results['bob_species_canonical_clean'] = {}
    results['bob_species_dedup'] = {}
    results['bob_polydt_classification'] = classification_counts
    results['bob_polydt_counts'] = {
        'ATCGCT': bob_counts['ATCGCT'],
        'TCGATA': bob_counts['TCGATA'],
        'c28_contam': classification_counts['c28_only'],
        'chimeric': classification_counts['chimeric_p1'] + classification_counts['chimeric_p2'],
        'both': classification_counts['both'],
    }
    # canonical-structure-only counts (TSO + insert + polyT + 6mer + bob-c).
    results['bob_polydt_canonical_counts'] = {
        'ATCGCT': canon_counts['ATCGCT'],
        'TCGATA': canon_counts['TCGATA'],
        'total': canon_counts['ATCGCT'] + canon_counts['TCGATA'],
    }
    results['bob_polydt_canonical_species'] = canon_species
    # canonical + mRNA-only counts (additionally requires aligned to
    # protein_coding transcript that is not MT and not a ribosomal protein).
    results['bob_polydt_canonical_mrna_counts'] = {
        'ATCGCT': canon_mrna_counts['ATCGCT'],
        'TCGATA': canon_mrna_counts['TCGATA'],
        'total': canon_mrna_counts['ATCGCT'] + canon_mrna_counts['TCGATA'],
    }
    results['bob_polydt_canonical_mrna_species'] = canon_mrna_species
    # mRNA-only counts (allows non-canonical structure — only requires
    # the read carries a 6mer code AND aligns to an mRNA transcript).
    results['bob_polydt_mrna_counts'] = {
        'ATCGCT': mrna_counts['ATCGCT'],
        'TCGATA': mrna_counts['TCGATA'],
        'total': mrna_counts['ATCGCT'] + mrna_counts['TCGATA'],
    }
    results['bob_polydt_mrna_species'] = mrna_species
    # mRNA + terminal-exon counts (allows non-canonical structure;
    # restricts to reads with terminal_exon genomic_class — the cleanest
    # 3'-end signal of polyadenylated mature mRNA).
    results['bob_polydt_mrna_terminal_counts'] = {
        'ATCGCT': mrna_term_counts['ATCGCT'],
        'TCGATA': mrna_term_counts['TCGATA'],
        'total': mrna_term_counts['ATCGCT'] + mrna_term_counts['TCGATA'],
    }
    results['bob_polydt_mrna_terminal_species'] = mrna_term_species
    # per-BC specificity by insert length (canonical reads, binned).
    spec_by_len = {}
    bc_swapped = bc in SWAPPED_BCS
    for name, _ in INSERT_LENGTH_BINS:
        bk = length_bin_counts[name]
        spec = compute_overall_specificity_from_counts(
            bk['a_h'], bk['a_m'], bk['t_h'], bk['t_m'], swapped=bc_swapped)
        spec_by_len[name] = {
            'atcgct_human': bk['a_h'],
            'atcgct_mouse': bk['a_m'],
            'tcgata_human': bk['t_h'],
            'tcgata_mouse': bk['t_m'],
            'n_total': bk['a_h'] + bk['a_m'] + bk['t_h'] + bk['t_m'],
            'overall_specificity': spec,
        }
    results['bob_polydt_specificity_by_insert_length'] = spec_by_len
    # per-species top-10 mRNA gene tables.
    # Computed on any read aligning to an mRNA (protein_coding, non-MT, non-
    # ribosomal-protein) with a mouse/human species call — NOT gated by
    # canonical / 6mer presence. This populates C28 single-species controls
    # that lack BOB-polydT chemistry too.
    results['bob_polydt_top_mrna_human'] = mrna_genes_human.most_common(10)
    results['bob_polydt_top_mrna_mouse'] = mrna_genes_mouse.most_common(10)
    # genomic-class breakdown of mRNA-aligned reads.
    # Diagnostic for polyT priming fidelity — should be dominated by
    # terminal_exon (the 3' end of polyA mRNA). Lower terminal_exon % =
    # more internal priming, intron retention, or annotation gaps.
    n_mrna_total = sum(mrna_genomic_class.values())
    results['bob_polydt_mrna_genomic_class_breakdown'] = {
        'n_total': n_mrna_total,
        'terminal_exon': mrna_genomic_class.get('terminal_exon', 0),
        'internal_exon': mrna_genomic_class.get('internal_exon', 0),
        'intron': mrna_genomic_class.get('intron', 0),
        'intergenic': mrna_genomic_class.get('intergenic', 0),
        'unmapped': mrna_genomic_class.get('unmapped', 0),
    }
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)

    n_real = bob_counts['ATCGCT'] + bob_counts['TCGATA']
    n_canon = canon_counts['ATCGCT'] + canon_counts['TCGATA']
    n_canon_mrna = canon_mrna_counts['ATCGCT'] + canon_mrna_counts['TCGATA']
    n_mrna = mrna_counts['ATCGCT'] + mrna_counts['TCGATA']
    n_mrna_term = mrna_term_counts['ATCGCT'] + mrna_term_counts['TCGATA']
    pct_canon = (100 * n_canon / n_real) if n_real else 0
    print(f'  {bc}: {n_real:,} BOB-polydT | '
          f'canonical: {n_canon:,} ({pct_canon:.0f}%) | '
          f'canon+mRNA: {n_canon_mrna:,} | '
          f'mRNA: {n_mrna:,} | '
          f'mRNA+terminal: {n_mrna_term:,} | '
          f'C28-contam={classification_counts["c28_only"]:,}')


def main():
    bcs = sorted(d for d in os.listdir(INDIV_DIR)
                 if d.startswith('barcode') and os.path.isdir(os.path.join(INDIV_DIR, d)))
    for bc in bcs:
        process_barcode(bc)
    print('\nDone. Run regen_cross_barcode_pdf_illumina.py to refresh the PDF.')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Authoritative bobcode / TSO registry, loaded from two reference sheets (located via
BOBCODE_REF_DIR or the default candidate directory below):

  24plex_bobcode_set_v1.tsv   -- the V1 24plex set (ID, Name, Barcode, Anchor, ...)
  bobcode_tso_oligo_list.tsv  -- the mixed stock list (tso-bob / BOB-J / polydT TSOs)

These two sheets ARE the bobcode configuration. This module turns them into a
lookup that maps a barcode SEQUENCE observed in reads back to its NAME / ID /
architecture / (intended) species, so the pipeline can auto-identify which
bobcodes an experiment used and print their names + sequences in the PDF.

Exposes:
  ARCH_SIGNATURES  -- ordered [(arch_label, search_motif)], longest/most-specific
                      first; the motif is a robust internal slice of the anchor
                      that survives ONT error and disambiguates the chemistry.
  BARCODES         -- {seq: {name, id, arch, species, anchor, ggg_offset, source}}
  build_lookup(mm) -- {variant_seq: canonical_seq} for exact + <=mm-mismatch match
  identify(seq)    -- (record, mismatches) best hit or (None, None)
  lengths_for(arch)-- sorted barcode lengths registered for an architecture
"""
import os, re, csv

# ---- locate the two authoritative sheets -----------------------------------
_CANDIDATE_DIRS = [
    os.environ.get('BOBCODE_REF_DIR', ''),
    os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))), 'Species Mixing Projects', 'all_read_summary'),
]

def _find_sheet(fname):
    for d in _CANDIDATE_DIRS:
        if d and os.path.exists(os.path.join(d, fname)):
            return os.path.join(d, fname)
    return None

# ---- architecture anchor signatures ----------------------------------------
# Order matters: the 24plex-D anchor (TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT) contains
# BOTH a Nextera-like head AND the TruSeq tail, so we test the TruSeq tail FIRST
# and only then the SMARTer / Nextera-only signatures.
ARCH_SIGNATURES = [
    ('24plex',            'GCTCTTCCGATCT'),      # TruSeq tail (A/B/C 7-mer + D extended)
    ('tso-bob (SMARTer)', 'GCAGTGGTATCAACGCAG'), # SMARTer core
    ('BOB-J (Nextera-R1)', 'GTGACTGGAGTTCAGACGT'),# Nextera Read1
]

BARCODES = {}   # seq -> record
_LOADED = False

def _add(seq, name, bid, arch, species=None, anchor=None, ggg=None, source=None):
    seq = (seq or '').upper().strip()
    if not seq or not re.fullmatch(r'[ACGT]{6,14}', seq):
        return
    # first writer wins, but fill in a species if a later row supplies one
    if seq in BARCODES:
        if species and not BARCODES[seq].get('species'):
            BARCODES[seq]['species'] = species
        return
    BARCODES[seq] = dict(seq=seq, name=name, id=bid, arch=arch,
                         species=species, anchor=anchor, ggg_offset=ggg,
                         source=source)

def _load_24plex(path):
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            bc = (row.get('Barcode') or '').strip().upper()
            if not re.fullmatch(r'[ACGT]{6,14}', bc):
                continue
            sett = (row.get('Set') or '').strip()
            arch = '24plex-D' if sett.startswith('D') else '24plex'
            try:
                ggg = int(row.get('GGG_offset_from_7mer_end') or 0)
            except ValueError:
                ggg = None
            _add(bc, row.get('Name', '').strip(), row.get('ID', '').strip(),
                 arch, anchor=(row.get('Anchor') or '').strip(), ggg=ggg,
                 source='24plex_v1')

# tso-bob / BOB-J barcodes are embedded in the oligo Sequence between the anchor
# and the rGrGrG tail; the Notes column states them EXPLICITLY ("human 7mer
# AGTACAT", "8bp barcode ATCGCTCC", "10bp barcode AAGTACGTTA"). Only an explicit
# barcode declaration counts -- a bare A/C/G/T run in Notes is usually a design
# element (capture tail "GGGCAT", "seq2 leader GCTAATCATTGC"), NOT a sample code.
_NOTE_PATTERNS = [
    re.compile(r'\b(?:human|mouse)\s+\d*\s*mer\s+([ACGT]{6,12})\b', re.I),  # "human 7mer AGTACAT"
    re.compile(r'\b\d+\s*bp\s+barcode\s+([ACGT]{6,12})\b', re.I),           # "8bp barcode ATCGCTCC"
    re.compile(r'\b\d+\s*mer\s+([ACGT]{6,12})\b', re.I),                    # "7mer ATCGAAA"
    re.compile(r'\bbarcode\s+([ACGT]{6,12})\b', re.I),                      # "barcode GCATGA"
]
_NOTE_SP = re.compile(r'(human|mouse)', re.I)

def _note_barcode(notes):
    for pat in _NOTE_PATTERNS:
        m = pat.search(notes or '')
        if m:
            return m.group(1).upper()
    return None

def _arch_for_oligo(name, seq, notes):
    s = (seq or '').upper()
    if 'GTGACTGGAGTTCAGACGT' in s or 'nextera' in (notes or '').lower():
        return 'BOB-J (Nextera-R1)'
    if 'GCAGTGGTATCAACGCAG' in s:
        return 'tso-bob (SMARTer)'
    return None

def _load_mixed(path):
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            name = (row.get('Name') or '').strip()
            seq = (row.get('Sequence') or '').strip()
            notes = (row.get('Notes') or '')
            cat = (row.get('Category') or '').strip()
            arch = _arch_for_oligo(name, seq, notes)
            if arch is None:
                continue
            sp = None
            m = _NOTE_SP.search(notes)
            if m:
                sp = m.group(1).lower()
            # only an EXPLICIT barcode declaration in Notes is trusted
            bc = _note_barcode(notes)
            if bc:
                _add(bc, name, name, arch, species=sp, source='oligo_list')

def load(force=False):
    global _LOADED
    if _LOADED and not force:
        return
    BARCODES.clear()
    p24 = _find_sheet('24plex_bobcode_set_v1.tsv')
    pmx = _find_sheet('bobcode_tso_oligo_list.tsv')
    # A sheet that cannot be found is an ERROR, not an empty registry. Otherwise a
    # mis-staged assets dir (or an unset BOBCODE_REF_DIR) loads nothing, every
    # identify() returns (None, None), and "no bobcodes in this data" looks like a
    # finding. The NEB index registry guards against the same failure class.
    if not p24 and not pmx:
        raise SystemExit(
            'bobcode_reference: neither 24plex_bobcode_set_v1.tsv nor '
            'bobcode_tso_oligo_list.tsv found. Searched: '
            + ', '.join(_CANDIDATE_DIRS) + '. Set BOBCODE_REF_DIR to the assets dir.')
    if p24:
        _load_24plex(p24)
    if pmx:
        _load_mixed(pmx)
    import sys as _sys
    _sys.stderr.write('bobcode_reference: loaded '
                      + ', '.join(x for x in (p24, pmx) if x)
                      + f' ({len(BARCODES)} barcodes)\n')
    _LOADED = True

# ---- fuzzy lookup ----------------------------------------------------------
_BASES = 'ACGT'

def build_lookup(max_mm=1):
    """{variant_seq: canonical_seq} covering exact + up to max_mm substitutions.
    Exact matches take precedence over 1-mismatch collisions."""
    load()
    exact = {s: s for s in BARCODES}
    if max_mm < 1:
        return exact
    fuzzy = {}
    for s in BARCODES:
        for i in range(len(s)):
            for b in _BASES:
                if b == s[i]:
                    continue
                v = s[:i] + b + s[i + 1:]
                if v not in exact:      # never shadow an exact barcode
                    fuzzy.setdefault(v, s)
    fuzzy.update(exact)
    return fuzzy

def lengths_for(arch=None):
    load()
    return sorted({len(s) for s, r in BARCODES.items()
                   if arch is None or r['arch'] == arch})

def _ham(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)

def identify(seq, max_mm=1):
    load()
    seq = seq.upper()
    if seq in BARCODES:
        return BARCODES[seq], 0
    best, bd = None, 99
    for s, r in BARCODES.items():
        if len(s) != len(seq):
            continue
        d = _ham(seq, s)
        if d <= max_mm and d < bd:
            best, bd = r, d
    return (best, bd) if best else (None, None)

if __name__ == '__main__':
    load()
    print(f"loaded {len(BARCODES)} barcodes from reference sheets")
    by_arch = {}
    for s, r in sorted(BARCODES.items(), key=lambda kv: (kv[1]['arch'], kv[1]['id'] or '')):
        by_arch.setdefault(r['arch'], []).append(r)
    for arch, recs in by_arch.items():
        print(f"\n[{arch}]  ({len(recs)} barcodes, lengths {lengths_for(arch)})")
        for r in recs:
            sp = f"  species={r['species']}" if r['species'] else ''
            print(f"  {r['seq']:12s}  {r['id'] or '':10s}  {r['name']:16s}{sp}")

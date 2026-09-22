#!/usr/bin/env python3
"""Chemistry adapter — gives star_taxonomy_illumina.py a uniform barcode/structure
interface for BOTH species-mixing chemistries, selected by
EXPERIMENT_CONFIG['barcode_location']:

  * 'BOBcode'  -> BOB-polydT: 6mer between the bob-c primer (…ATGTC) and polyT
                  (override_bob_polydt.measure_structure; codes p1/p2).
  * 'TSO'      -> TSO-bob: 7mer inside the TSO, after the Nextera-R1 backbone,
                  before GGG (override_bob_tso_7mer.measure_struct_tso; codes h/m).

The species barcode SEQUENCES live in the override modules (detection logic);
the code->species mapping (the correctness oracle) is taken from config
(primary_species_map = {seq: 'human'/'mouse'}). Imports of the override modules
are LAZY so this file works in an experiment dir that only ships one of them.

get_adapter(barcode_location, species_map) -> Adapter, where Adapter exposes:
  .measure(seq)  -> dict with 'code' (key of EXP/CODE_SEQ) + structural fields, or None
  .EXP           -> {code: 'human'/'mouse'}            (expected species per barcode)
  .CODE_SEQ      -> {code: display sequence}
  .FUNNEL        -> [(key, label, gate(m, aln_len)->bool), ...] cumulative funnel;
                    the last level == "Full structure" (the retained set)
  .LEN_HIST      -> [(measure-field, label), ...] numeric fields to histogram
  .label         -> short chemistry label for captions
"""


class Adapter:
    def __init__(self, measure, EXP, CODE_SEQ, FUNNEL, LEN_HIST, label, ELEMENTS=None,
                 ret_key='bob_polydt_retention', no_code_reason='no_code',
                 retained_fn=None, multiplicity=None):
        self.measure = measure
        # multiplicity(seq) -> {'n_tso','codes','orientations'} or None. Set for
        # chemistries where a read can carry >1 species barcode (TSO-TSO chimeras);
        # star_taxonomy aggregates it into the TSO-structure multiplicity stats.
        self.multiplicity = multiplicity
        self.EXP = EXP
        self.CODE_SEQ = CODE_SEQ
        self.FUNNEL = FUNNEL
        self.LEN_HIST = LEN_HIST
        self.label = label
        # ELEMENTS: [(key, label, detector(seq)->bool), ...] structural-element
        # detectors run over ALL reads for the Master-Summary "% reads with …"
        # rows. Chemistry-specific so the rows are correct for each construct.
        self.ELEMENTS = ELEMENTS or []
        # Read-structure retention (§1): JSON key the structure-figure module
        # reads, the drop-reason for reads with no barcode, and the per-read
        # retained_call(m, aln_len) -> (kept, reason) (reason strings match the
        # structure module's drop-reason buckets).
        self.RET_KEY = ret_key
        self.NO_CODE_REASON = no_code_reason
        self._retained_fn = retained_fn

    def retained_reason(self, m, aln_len):
        """(kept, reason) for a measure() result; reason ∈ the structure-figure
        drop-reason buckets. aln_len = STAR aligned length (used as the insert)."""
        if self._retained_fn is None:
            return False, self.NO_CODE_REASON
        return self._retained_fn(m, aln_len)

    def funnel_level(self, m, aln_len):
        """Index of the deepest cumulative funnel level the read reaches
        (0..len(FUNNEL)-1). retained iff == len(FUNNEL)-1."""
        lvl = -1
        for i, (_k, _lbl, gate) in enumerate(self.FUNNEL):
            if gate(m, aln_len):
                lvl = i
            else:
                break
        return lvl


# Canonical 24plex 7-mer allowlist — a FALLBACK only. New runs should declare the
# chemistry explicitly via the guide field `tso_arch` (24plex | tso-bob | polydt);
# this list only rescues legacy guides written before tso_arch existed. Single source
# of truth so the chemistry selector can't drift across call sites.
KNOWN_24PLEX_7MERS = frozenset({'AGTGGTG', 'AGCAGAC', 'TACAAGC', 'ACTCACA',
                                'GATATGG', 'TCGTGAC'})


def is_24plex(species_map=None, tso_arch=None):
    """True for the 24plex R1-anchor TSO chemistry. Prefers the explicit guide field
    `tso_arch`; falls back to detecting a known 24plex 7-mer in the species map for
    legacy guides that predate tso_arch."""
    if tso_arch:
        return str(tso_arch).strip().lower() in ('24plex', '24-plex')
    return any(m in (species_map or {}) for m in KNOWN_24PLEX_7MERS)


def get_adapter(barcode_location, species_map, tso_arch=None):
    species_map = species_map or {}
    if (barcode_location or 'BOBcode') == 'TSO':
        # 24plex (R1 anchor + 7-mer + UMI + GGG) vs classic AGTACAT/ACCTTGA TSO-bob,
        # selected by is_24plex (explicit guide tso_arch, else the known-7mer allowlist).
        if is_24plex(species_map, tso_arch):
            return _24plex_adapter(species_map)
        return _tso_adapter(species_map)
    return _polydt_adapter(species_map)


def _24plex_adapter(species_map):
    import override_24plex_illumina as ox
    # Include ONLY a side whose 7-mer this run actually declares.
    # analyze_speciesmix sets ox.BC_HUMAN/BC_MOUSE from the run's species map but
    # guards each with `if _hu:` / `if _ms:`, so a SINGLE-SPECIES tube (one code in
    # the map) leaves the other side at its module default. Building the code
    # table from both would inject a PHANTOM 7-mer that is not in the library: it
    # would appear in the crosstab with zero reads, and detect_species_map would fail
    # the run with "too few reads (n=0)" on a barcode that was never there. Two-code
    # runs declare both sides, so both are included.
    _sides = [('h', ox.BC_HUMAN), ('m', ox.BC_MOUSE)]
    code_seq = {k: v for k, v in _sides if v in species_map}
    if not code_seq:                       # legacy guide with no/unknown map: keep both
        code_seq = {k: v for k, v in _sides}
    exp = {k: species_map.get(v, 'human' if k == 'h' else 'mouse')
           for k, v in code_seq.items()}
    funnel = [
        ('code', 'Has 7mer barcode', lambda m, a: True),
        ('code_insert', '+ insert >=50 bp',
         lambda m, a: a is not None and ox.RETAIN_INSERT_MIN <= a < ox.RETAIN_INSERT_MAX),
        # DELEGATE to retained_call_24plex rather than re-testing rt_end here, so the
        # retention rule has one definition: a separate copy would not know about read
        # geometry and would report ~0% 'Full structure' on a run whose R1 is simply too
        # short to show a polydT tract.
        ('retained',
         ('+ R1-polydT RT end = Full structure' if getattr(ox, 'POLYT_MEASURABLE', True)
          else '+ insert = Full structure (R1-polydT not in the read at this length)'),
         lambda m, a: ox.retained_call_24plex(m, a)[0]),
    ]
    # NOTE: GGG template-switch removed as a funnel/retention requirement — the two
    # 24plex barcodes use different-length post-UMI spacers (human HH -> GGG at +8;
    # mouse HHMWMWMW -> GGG at +14), so GGG is reliably visible only in a minority of
    # mouse reads. It remains reported as an informational element (the GGG row).
    F, R = ox.BACKBONE_F, ox.BACKBONE_R
    elements = [
        ('backbone', '% reads with R1 backbone', lambda s: F in s or R in s),
        ('ggg', '% reads with GGG switch', ox.has_ggg_element),
        ('polyt', '% reads with polyT', lambda s: ('T' * 10 in s) or ('A' * 10 in s)),
        ('rt_end', '% reads with RT-end (R1-polydT)', ox.has_rt_end),
    ]
    return Adapter(ox.measure_struct_24plex, exp, code_seq, funnel,
                   [('backbone_bp', 'R1 backbone (bp)'), ('polyt', 'polyT/A run (bp)')],
                   'TSO-24plex 7mer', ELEMENTS=elements,
                   ret_key='tso7_retention', no_code_reason='no_barcode',
                   retained_fn=lambda m, a: ox.retained_call_24plex(m, a),
                   multiplicity=ox.tso_unit_summary_24plex)


def _polydt_adapter(species_map):
    import override_bob_polydt_illumina as ob
    code_seq = {'p1': ob.P1_6MER, 'p2': ob.P2_6MER}
    exp = {'p1': species_map.get(ob.P1_6MER, 'human'),
           'p2': species_map.get(ob.P2_6MER, 'mouse')}
    # Gates follow override_bob_polydt's canonical definition (bob-c primer +
    # 6mer + fuzzy polyT + aligned insert; no TSO/C28 requirement). The insert for
    # polydT is the STAR-aligned length `a`, not a field of m.
    funnel = [
        ('code', 'Has code (bobcode)', lambda m, a: True),
        ('code_bobc', '+ full bob-c primer (20 bp)',
         lambda m, a: bool(m.get('bobc_full'))),
        ('code_bobc_polyt', f'+ polyT (>={ob.RETAIN_POLYT_MIN} of {ob.POLYT_WINDOW} nt T)',
         lambda m, a: m.get('polyt', 0) >= ob.RETAIN_POLYT_MIN),
        ('retained', f'+ insert >={ob.RETAIN_INSERT_MIN} bp (STAR) = Full structure',
         lambda m, a: ob.retained_call(m, a)[0]),
    ]
    elements = [
        ('bobc', '% reads with bob-c primer (20 bp)', ob.has_full_bobc_primer),
        ('polyt', '% reads with 6mer anchor + polyT',
         lambda s: ob.classify_read(s) not in ('neither', 'c28_only')),
        ('c28', '% reads with C28', lambda s: (ob.C28_FWD in s) or (ob.C28_RC in s)),
        ('tso', '% reads with TSO', lambda s: (ob.TSO_FWD in s) or (ob.TSO_RC in s)),
    ]
    return Adapter(ob.measure_structure, exp, code_seq, funnel,
                   [('polyt', f'polyT (T in {ob.POLYT_WINDOW}-nt window)')], 'BOB-polydT 6mer',
                   ELEMENTS=elements,
                   ret_key='star_retention', no_code_reason='no_code',
                   retained_fn=lambda m, a: ob.retained_call(m, a))
    # polydT keeps its analyze_barcode-derived "% reads with …" rows (no ELEMENTS),
    # and its §1 still reads the override-written bob_polydt_retention (which has
    # c28_bins etc.); star_taxonomy's retention goes to the unused 'star_retention'
    # key so it doesn't clobber the override. (TSO §1 uses the STAR retention.)


def _tso_adapter(species_map):
    import override_bob_tso_7mer_illumina as ot
    code_seq = {'h': ot.BC_HUMAN, 'm': ot.BC_MOUSE}
    exp = {'h': species_map.get(ot.BC_HUMAN, 'human'),
           'm': species_map.get(ot.BC_MOUSE, 'mouse')}
    # insert length for TSO = the STAR aligned length (passed as `a`).
    funnel = [
        ('code', 'Has 7mer barcode', lambda m, a: True),
        ('code_ggg', '+ GGG template-switch', lambda m, a: bool(m['ggg'])),
        ('code_ggg_insert', '+ insert ≥50 bp',
         lambda m, a: a is not None and ot.RETAIN_INSERT_MIN <= a < ot.RETAIN_INSERT_MAX),
        ('retained', '+ C28 end = Full structure', lambda m, a: bool(m['c28'])),
    ]
    F, R = ot.BACKBONE_F, ot.BACKBONE_R

    def _has_ggg(s):
        i = s.find(F)
        while i >= 0:
            if s[i + len(F) + 7:i + len(F) + 10].count('G') >= 2:
                return True
            i = s.find(F, i + 1)
        i = s.find(R)
        while i >= 0:
            if s[i - 10:i - 7].count('C') >= 2:   # rc(GGG)=CCC, 5' of the rc-backbone
                return True
            i = s.find(R, i + 1)
        return False

    elements = [
        ('backbone', '% reads with TSO backbone', lambda s: F in s or R in s),
        ('ggg', '% reads with GGG switch', _has_ggg),
        ('polyt', '% reads with polyT', lambda s: ('T' * 10 in s) or ('A' * 10 in s)),
        ('c28', '% reads with RT-end (C28 or R1)', ot.has_c28_end),
    ]
    return Adapter(ot.measure_struct_tso, exp, code_seq, funnel,
                   [('backbone_bp', 'TSO backbone (bp)'), ('polyt', 'polyT/A run (bp)')],
                   'TSO-bob 7mer', ELEMENTS=elements,
                   ret_key='tso7_retention', no_code_reason='no_barcode',
                   retained_fn=lambda m, a: ot.retained_call_tso(m, a),
                   multiplicity=ot.tso_unit_summary)

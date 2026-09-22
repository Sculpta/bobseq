#!/usr/bin/env python3
"""
Analysis for RT+TSO constructs (barcodes 01-04) — species mix experiment.
Read structure: 5'-[RTprimer_HR2: ATACTCGTGAC][ATCAGCTTCC][polyT]...[cDNA]...[nontemplated nt][TSO]-3'

TSO variants (4bp barcode at end):
  GCAGTGGTATCAACGCAGAGTACA  → TSO barcode: TACA (→ mouse)
  GCAGTGGTATCAACGCAGAGGGCA  → TSO barcode: GGCA (→ human)

In nanopore reads, orientation is random so TSO can appear as fwd or RC.
"""

import sys, os, json, re
from datetime import datetime
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def load_sequence_list(path):
    """Parse sequence_list_illumina.txt (or dataset-specific YYMMDD_sequence_list_illumina.txt).
    Returns tuple: (flat_dict, sections_dict).
      flat_dict: name -> sequence (backward-compatible as GROUND_TRUTH)
      sections_dict: section_name -> list of (name, seq, description, usage)
    """
    flat = {}
    sections = {}
    current_section = 'Uncategorized'
    sections[current_section] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line or line.startswith('#'):
                continue
            # Section header
            if line.startswith('[') and line.endswith(']'):
                current_section = line[1:-1]
                sections.setdefault(current_section, [])
                continue
            # Tab-separated entry: Name<tab>Sequence[<tab>Description[<tab>Usage]]
            if '\t' in line:
                parts = line.split('\t')
                name = parts[0].strip()
                seq = parts[1].strip() if len(parts) > 1 else ''
                desc = parts[2].strip() if len(parts) > 2 else ''
                usage = parts[3].strip() if len(parts) > 3 else ''
                if name and seq:
                    flat[name] = seq
                    sections[current_section].append((name, seq, desc, usage))
    # Remove empty Uncategorized section if unused
    if not sections.get('Uncategorized'):
        sections.pop('Uncategorized', None)
    return flat, sections


def find_sequence_list_path():
    """Find sequence_list_illumina.txt. Search order:
    1. Parent of script dir (SpeciesMixing/ or experiment root)
    2. Script dir itself
    3. machine_readable_guide/ sibling for dataset-specific YYMMDD_sequence_list_illumina.txt
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(script_dir)
    # Parent dir (e.g. SpeciesMixing/ when running from template, or experiment root)
    p = os.path.join(parent_dir, 'sequence_list_illumina.txt')
    if os.path.exists(p):
        return p
    # Script dir itself
    p = os.path.join(script_dir, 'sequence_list_illumina.txt')
    if os.path.exists(p):
        return p
    # Dataset-specific: look for YYMMDD_sequence_list_illumina.txt in machine_readable_guide/
    mrg_dir = os.path.join(parent_dir, 'machine_readable_guide')
    if os.path.isdir(mrg_dir):
        for f in sorted(os.listdir(mrg_dir), reverse=True):
            if f.endswith('_sequence_list_illumina.txt'):
                return os.path.join(mrg_dir, f)
    raise FileNotFoundError(f"No sequence_list_illumina.txt found in {parent_dir} or machine_readable_guide/")


def load_experiment_config(guide_path):
    """Parse experiment config from analysis guide (between BEGIN_CONFIG/END_CONFIG).
    Returns dict with parsed values. All experiment-specific configuration that would
    otherwise be hardcoded in the script is read from here.

    Expected keys:
      barcodes             comma-separated barcode names (e.g. barcode05,barcode06)
      fastq_dir            name of FASTQ directory (e.g. fastq, fastq_pass)
      experiment_title     title for cross-barcode summary PDF
      construct_description  one-line construct description for reports (optional; auto-generated if absent)
      barcode_location     'BOBcode' or 'TSO' — where the species barcode is located (default: BOBcode)
      bob_species_map      comma-separated CODE=species pairs (e.g. ATCG=human,CCGA=mouse)
      tso_species_map      comma-separated CODE=species pairs (e.g. TACA=mouse,GGCA=human)
      rt_primers_used      comma-separated RT primer names (e.g. PolydT_5P_C28,C28-RT-N9)
      condition_labels     comma-separated barcode=label pairs (e.g. barcode05=1x)
      transcriptome_ref    path to combined_transcriptome.fa
      genome_ref           path to combined_genome.fa
      mouse_gtf            path to mouse GTF
      human_gtf            path to human GTF
    """
    raw = {}
    in_block = False
    with open(guide_path) as f:
        for line in f:
            line = line.rstrip()
            if line == 'BEGIN_CONFIG':
                in_block = True
                continue
            if line == 'END_CONFIG':
                break
            if in_block and '\t' in line:
                key, val = line.split('\t', 1)
                raw[key.strip()] = val.strip()

    config = {}
    config['barcodes'] = [b.strip() for b in raw.get('barcodes', '').split(',') if b.strip()]
    config['fastq_dir'] = raw.get('fastq_dir', 'fastq')
    config['experiment_title'] = raw.get('experiment_title', 'Species Mix Analysis')
    config['construct_description'] = raw.get('construct_description', '')
    config['transcriptome_ref'] = raw.get('transcriptome_ref', '')
    config['genome_ref'] = raw.get('genome_ref', '')
    config['mouse_gtf'] = raw.get('mouse_gtf', '')
    config['human_gtf'] = raw.get('human_gtf', '')

    # Parse bob_species_map: "ATCG=human,CCGA=mouse" -> {'ATCG': 'human', 'CCGA': 'mouse'}
    config['bob_species_map'] = {}
    for pair in raw.get('bob_species_map', '').split(','):
        if '=' in pair:
            code, species = pair.strip().split('=', 1)
            config['bob_species_map'][code.strip()] = species.strip()

    # Parse barcode_location: 'BOBcode' or 'TSO' — where the species barcode is located
    config['barcode_location'] = raw.get('barcode_location', 'BOBcode')
    # Explicit TSO architecture (24plex | tso-bob | polydt). Preferred over inferring
    # the chemistry from the barcode 7-mers; blank falls back to that inference.
    config['tso_arch'] = raw.get('tso_arch', '')

    # Parse tso_species_map: same format as bob_species_map
    config['tso_species_map'] = {}
    for pair in raw.get('tso_species_map', '').split(','):
        if '=' in pair:
            code, species = pair.strip().split('=', 1)
            config['tso_species_map'][code.strip()] = species.strip()

    # Parse grun_offsets: "GATATGG=8,TCGTGAC=11" -> {barcode_seq: int} (per-barcode
    # GGG position after the barcode; drives filter_funnel's generic G-run caller).
    config['grun_offsets'] = {}
    for pair in raw.get('grun_offsets', '').split(','):
        if '=' in pair:
            s, o = pair.split('=', 1)
            try:
                config['grun_offsets'][s.strip()] = int(o.strip())
            except ValueError:
                pass

    # Parse rt_primers_used: comma-separated list of RT primer canonical names
    config['rt_primers_used'] = [p.strip() for p in raw.get('rt_primers_used', '').split(',') if p.strip()]
    config['rt_end_label'] = raw.get('rt_end_label', 'C28')   # PDF RT-end label (R1 / C28)

    # Derived: primary_species_map is whichever map matches barcode_location
    if config['barcode_location'] == 'TSO':
        if not config['tso_species_map']:
            raise ValueError("barcode_location=TSO but tso_species_map is empty in experiment config")
        config['primary_species_map'] = config['tso_species_map']
        config['secondary_species_map'] = config['bob_species_map']
    else:
        config['primary_species_map'] = config['bob_species_map']
        config['secondary_species_map'] = config['tso_species_map']

    # Parse condition_labels: "barcode05=1x,barcode06=0.5x" -> {'barcode05': '1x', ...}
    config['condition_labels'] = {}
    for pair in raw.get('condition_labels', '').split(','):
        if '=' in pair:
            bc, label = pair.strip().split('=', 1)
            config['condition_labels'][bc.strip()] = label.strip()

    # Parse barcode_species_overrides: per-BC species map overrides. Format:
    #   barcode17:ATCGCT=mouse;TCGATA=human,barcode18:ATCGCT=mouse;TCGATA=human,...
    # Used when subset of BCs have swapped/flipped chemistry codes.
    # Each entry maps a barcode name to a dict that REPLACES the default
    # primary_species_map for that BC only.
    config['barcode_species_overrides'] = {}
    raw_overrides = raw.get('barcode_species_overrides', '')
    if raw_overrides:
        for bc_entry in raw_overrides.split(','):
            bc_entry = bc_entry.strip()
            if ':' not in bc_entry:
                continue
            bc, pairs = bc_entry.split(':', 1)
            bc = bc.strip()
            override_map = {}
            for pair in pairs.split(';'):
                if '=' in pair:
                    code, species = pair.strip().split('=', 1)
                    override_map[code.strip()] = species.strip()
            if override_map:
                config['barcode_species_overrides'][bc] = override_map

    # MAPQ filter (default 0 = no filter). Applied at SAM-parse time in both
    # transcriptome and genome alignments. Set in BEGIN_CONFIG with `min_mapq	N`.
    try:
        config['min_mapq'] = int(raw.get('min_mapq', '0') or 0)
    except ValueError:
        config['min_mapq'] = 0

    # ILLUMINA platform fields. This loader is an explicit WHITELIST — any guide field
    # not copied here is silently dropped, which is why these two must be listed: without
    # them EXPERIMENT_CONFIG['read_layout'] is None, IS_ILLUMINA is False, and every
    # Illumina behaviour (positional barcode detection, dropped/greyed rows) quietly
    # stays off while the report still renders.
    config['platform'] = raw.get('platform', '')
    config['read_layout'] = raw.get('read_layout', '')
    # Geometry-dependent elements this run's READ LENGTH cannot observe,
    # as a comma-separated list. prep_reads derives the same set from the measured
    # reads and refuses to run unless the two agree, so by the time the analysis sees
    # this field it has already been checked against the data.
    config['unmeasurable'] = raw.get('unmeasurable', '')

    return config


def build_construct_description(config):
    """Auto-generate construct description string from barcode_location + rt_primers_used.
    If construct_description is already set in config, returns it as-is (manual override).
    """
    override = config.get('construct_description', '')
    if override:
        return override

    loc = config.get('barcode_location', 'BOBcode')
    rt_primers = config.get('rt_primers_used', [])
    has_c28 = any('C28' in p for p in rt_primers)
    has_polydt = any('polydT' in p.lower() or 'polydt' in p.lower() for p in rt_primers)

    if loc == 'BOBcode':
        parts = ["5'-[BOB-c primer][4bp BOBcode][BOB_HR][RT_HR]"]
        if has_c28:
            parts.append("[C28 linker]")
        if has_polydt:
            parts.append("[polyT]")
        parts.append("...[cDNA]...[NTs][TSO primer]-3'")
        return ''.join(parts)
    else:
        parts = ["5'-[TSO primer][4bp TSO barcode][NTs]...[cDNA]..."]
        if has_polydt:
            parts.append("[polyT]")
        if has_c28:
            parts.append("[C28 linker]")
        parts.append("[RT_HR]-3'")
        return ''.join(parts)


def build_construct_diagram(config, results):
    """Build elements_diag list for the construct structure diagram.
    Returns list of (label, size_bp, color) tuples in 5'→3' order.
    """
    loc = config.get('barcode_location', 'BOBcode')
    cDNA_len = int(results.get('len_mean', 200))
    rt_primers = config.get('rt_primers_used', [])
    has_c28 = any('C28' in p for p in rt_primers)
    has_polydt = any('polydT' in p.lower() or 'polydt' in p.lower() for p in rt_primers)

    if loc == 'BOBcode':
        elems = [
            ('BOB-c\nprimer', 20, '#2CA02C'),
            ('BOB\ncode', 4, '#FFC000'),
            ('BOB_HR', 10, '#98DF8A'),
            ('RT_HR', 11, '#ED7D31'),
        ]
        if has_c28:
            elems.append(('C28\nlinker', 10, '#A5A5A5'))
        if has_polydt:
            elems.append(('polyT', 30, '#FF6B6B'))
        elems += [
            ('cDNA', cDNA_len, '#9B59B6'),
            ('NTs', 5, '#CCCCCC'),
            ('TSO\nprimer', 18, '#E74C3C'),
        ]
        return elems
    else:
        elems = [
            ('TSO\nprimer', 18, '#E74C3C'),
            ('TSO\nbarcode', 4, '#FFC000'),
            ('NTs', 5, '#CCCCCC'),
            ('cDNA', cDNA_len, '#9B59B6'),
        ]
        if has_polydt:
            elems.append(('polyT', 30, '#FF6B6B'))
        if has_c28:
            elems.append(('C28\nlinker', 10, '#A5A5A5'))
        elems.append(('RT_HR', 11, '#ED7D31'))
        return elems


def create_dataset_sequence_list(master_path, output_dir, config, ground_truth, rc_func):
    """Copy master sequence_list_illumina.txt to output_dir/guide/ as YYMMDD_sequence_list_illumina.txt,
    add timestamp header and compute derivative sequences section.
    Returns path to the created file.
    """
    date_prefix = datetime.now().strftime('%y%m%d')
    guide_dir = os.path.join(output_dir, 'machine_readable_guide')
    os.makedirs(guide_dir, exist_ok=True)
    output_path = os.path.join(guide_dir, f'{date_prefix}_sequence_list_illumina.txt')

    # Read master content
    with open(master_path) as f:
        master_content = f.read()

    # IDEMPOTENCY: if the "master" already carries a [Derivative Sequences]
    # section (i.e. it's a previously-generated dataset copy), strip it before
    # regenerating. Otherwise reverse-complement derivatives compound every run
    # (_RC -> _RC_RC -> _RC_RC_RC ...), bloating the searchable list and slowing
    # the per-read scan ~Nx. Always regenerate derivatives from the base list.
    _cut = master_content.find('[Derivative Sequences]')
    if _cut != -1:
        master_content = master_content[:_cut].rstrip() + '\n'

    # Compute derivative sequences.
    #
    # IMPORTANT: derive from the BASE entries parsed out of the (already
    # derivative-stripped) master_content — NOT from the passed-in `ground_truth`.
    # `ground_truth` is loaded from the dataset copy, which carries the PREVIOUS
    # run's derivatives; deriving from it re-appended `_RC` on every regeneration,
    # compounding `_RC -> _RC_RC -> _RC_RC_RC ...` (and duplicating the 4bp
    # BOB_barcode_/TSO_barcode codes once per compounded copy). Parsing the base list
    # here makes derivative generation idempotent regardless of how many times it runs.
    base = {}
    for _bl in master_content.splitlines():
        _bl = _bl.rstrip()
        if not _bl or _bl.startswith('#') or (_bl.startswith('[') and _bl.endswith(']')):
            continue
        if '\t' in _bl:
            _bp = _bl.split('\t')
            _bn = _bp[0].strip()
            _bs = _bp[1].strip() if len(_bp) > 1 else ''
            if _bn and _bs and not (_bn.endswith('_RC') or '_RC_' in _bn):
                base[_bn] = _bs

    derivatives = []

    # 4bp BOB barcode codes
    for gt_name, gt_seq in base.items():
        if gt_name.startswith('BOB-C5-') and len(gt_seq) == 34 and '[' not in gt_seq:
            code = gt_seq[20:24]
            derivatives.append((f'BOB_barcode_{code}', code,
                f'4bp BOB barcode from {gt_name} positions 20-23', 'bob_barcode_extraction'))

    # 4bp TSO barcode codes
    for gt_name in ('TSO-seq1-variant', 'TSO-seq2-variant'):
        full = base.get(gt_name, '')
        if full and len(full) >= 24:
            code = full[20:24]
            derivatives.append((f'TSO_barcode_{code}', code,
                f'4bp TSO barcode from {gt_name} positions 20-23', 'tso_barcode_extraction'))

    # Reverse complements of all pure-DNA base entries (base excludes anything
    # already named *_RC, so _RC can never compound)
    for name, seq in base.items():
        if '[' in seq or ']' in seq:
            continue  # skip pattern descriptors
        if any(c not in 'ATGCN' for c in seq.upper()):
            continue  # skip non-DNA
        rc_seq = rc_func(seq)
        derivatives.append((f'{name}_RC', rc_seq,
            f'Reverse complement of {name}', 'detection_rc'))

    # Composite sequences
    rt_hr = base.get('RT_HR', '')
    c28_link = base.get('C28_linker', '')
    if rt_hr and c28_link:
        rt_c28_prefix = rt_hr + c28_link + 'T'
        derivatives.append(('RT_C28_PREFIX', rt_c28_prefix,
            f'RT_HR + C28_linker + T ({len(rt_c28_prefix)}bp) for C28 subtype detection', 'rt_subtype_detection'))
        rt_c28_6t = rt_c28_prefix + 'TTTTT'
        derivatives.append(('RT_C28_6T', rt_c28_6t,
            f'RT_C28_PREFIX + 5T ({len(rt_c28_6t)}bp) for polyT C28 detection', 'rt_subtype_detection'))

    # Write dataset-specific file
    with open(output_path, 'w') as f:
        f.write(f'# Dataset-specific sequence list copied from {os.path.basename(master_path)} '
                f'on {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n')
        f.write(master_content)
        if not master_content.endswith('\n'):
            f.write('\n')
        f.write('\n[Derivative Sequences]\n')
        for name, seq, desc, usage in derivatives:
            f.write(f'{name}\t{seq}\t{desc}\t{usage}\n')

    return output_path


def find_guide_path():
    """Find the experiment-specific analysis guide.
    Search order:
    1. machine_readable_guide/ sibling of scripts/ (in experiment dir)
    2. Parent dir (SpeciesMixing/) for SpeciesMix_guide.txt (template mode)
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(script_dir)
    # Experiment dir: machine_readable_guide/ contains *_analysisguide.txt or *_guide.txt
    mrg_dir = os.path.join(parent_dir, 'machine_readable_guide')
    if os.path.isdir(mrg_dir):
        for f in sorted(os.listdir(mrg_dir)):
            if f.endswith('_analysisguide.txt') or f.endswith('_guide.txt'):
                return os.path.join(mrg_dir, f)
    # Template mode: SpeciesMix_guide.txt in parent dir
    for f in sorted(os.listdir(parent_dir)):
        if f.endswith('_guide.txt') or f.endswith('_analysisguide.txt'):
            return os.path.join(parent_dir, f)
    raise FileNotFoundError(f"No experiment guide found in {mrg_dir} or {parent_dir}")


# Load ground truth and experiment config from analysis guide
GUIDE_PATH = find_guide_path()
SEQ_LIST_PATH = find_sequence_list_path()
GROUND_TRUTH, SEQUENCE_SECTIONS = load_sequence_list(SEQ_LIST_PATH)
EXPERIMENT_CONFIG = load_experiment_config(GUIDE_PATH)

# --- ILLUMINA switch --------------------------------------------------------
# True when the run is Illumina (registry: platform=illumina_pe150 / read_layout=
# r2_positional). Sections and metrics that measure a property of an ONT read (a
# read spanning the WHOLE molecule, ONT adapters/native barcode, molecule-length
# inserts) are either dropped or greyed off this flag. Kept as one predicate so
# every such decision is greppable and reversible.
class _IllRemoved(Exception):
    """Raised to skip a figure/section REMOVED on Illumina; the enclosing try/except
    logs the reason so the removal is visible in the run log."""


IS_ILLUMINA = (str(EXPERIMENT_CONFIG.get('platform', '')).lower().startswith('illumina')
               or EXPERIMENT_CONFIG.get('read_layout') == 'r2_positional')

# --- METHODS PROVENANCE ------------------------------------------------------
# The Methods facts (aligner, sequencer, duplicate definition) are not fixed strings:
# a report can render perfectly while describing a different pipeline, and corrected
# strings drift again the next time the pipeline changes. They are DERIVED from the
# same values that drive the analysis, so the Methods section cannot disagree with
# what actually ran.
_PLATFORM_RAW = str(EXPERIMENT_CONFIG.get('platform', '') or '').strip()

METHODS_ALIGNER = 'STAR 2.7.11b' if IS_ILLUMINA else 'STARlong 2.7.11b'
METHODS_ALIGNER_BIN = 'STAR' if IS_ILLUMINA else 'STARlong'
METHODS_ALIGNER_TUNING = (
    'short-read parameters (local alignment, mismatch ratio <=0.04, score and '
    'match-length >=0.66 of read length, multimapping retained)' if IS_ILLUMINA else
    'long-read-tuned parameters (local alignment, relaxed score/match-length and '
    'mismatch ratios, multimapping retained)')
METHODS_SEQUENCER = (
    (f'Illumina paired-end (platform={_PLATFORM_RAW})' if _PLATFORM_RAW
     else 'Illumina paired-end') if IS_ILLUMINA
    else (EXPERIMENT_CONFIG.get('sequencer') or
          'Oxford Nanopore (flow cell chemistry and basecaller not recorded in the run config; add `sequencer` to runs/&lt;run&gt;.json)'))

# Read the ACTUAL dedup constants rather than restating them. If that import fails
# the value is left visibly unresolved on purpose: a plausible-looking wrong string
# is exactly what this block exists to prevent, so an obvious placeholder is safer.
if IS_ILLUMINA:
    try:
        import umi_dedup_illumina as _ud_meta
        _ud_meta.configure_umi(EXPERIMENT_CONFIG.get('rt_primers_used'))
        if _ud_meta.UMI_SOURCE == '6N+HH':
            # random priming, paired-end: 6N + HH = 8 nt, the benchmark's definition; the read name
            # carries the whole degenerate stretch, dedup uses its first 8 nt; the key carries the fragment length
            _umi_desc = 'UMI = 6N + HH (8 nt; the first 8 nt of the degenerate stretch after the bobcode, 7 for the 11-mer BOB25D)'
        elif _ud_meta.UMI_SOURCE == '6N+spacer':
            _umi_desc = ('UMI = the 6N plus the degenerate spacer between bobcode and G-run, 8 / 11 / 14 nt '
                         'by code set (7 nt for the 11-mer BOB25D)')
        else:
            _umi_desc = f'{_ud_meta.UMI_SOURCE} UMI ({"26 nt R1 + 6 nt R2" if _ud_meta.UMI1_LEN else "6 nt TSO only"})'
        METHODS_DUPLICATE = (f'{_umi_desc}, clustered at edit distance {_ud_meta.UMI_MAX_ED}; '
                             f'deduplicated before analysis')
        DEDUP_KEY_EXTRA = (', fragment length (both mates aligned paired-end; 0 when the mate is unmapped)'
                           if _ud_meta.UMI_SOURCE in ('6N+HH', '6N+spacer') else '')
    except Exception as _e_meta:                      # noqa: BLE001
        METHODS_DUPLICATE = f'UNRESOLVED (umi_dedup_illumina unavailable: {_e_meta})'
        DEDUP_KEY_EXTRA = ''
else:
    METHODS_DUPLICATE = 'First + last 20bp of insert'
    DEDUP_KEY_EXTRA = ''

# Auto-generate construct description if not set in config
if not EXPERIMENT_CONFIG.get('construct_description'):
    EXPERIMENT_CONFIG['construct_description'] = build_construct_description(EXPERIMENT_CONFIG)

# 24plex detector: drive override_24plex's detected 7-mers from the config so a run
# can use a different species barcode (e.g. a third, human-only 7-mer) without
# editing the module. Human/mouse = first of each in the primary species map.
# Chemistry is decided by barcode_adapter.is_24plex (single source of truth): the guide
# field `tso_arch` if present, else the known-7mer allowlist.
import re as _re
import barcode_adapter_illumina as _ba
_pm = EXPERIMENT_CONFIG.get('primary_species_map', {})
_hu = next((k for k, v in _pm.items() if v == 'human'), None)
_ms = next((k for k, v in _pm.items() if v == 'mouse'), None)
_is24 = (EXPERIMENT_CONFIG.get('barcode_location') == 'TSO' and _ba.is_24plex(
    EXPERIMENT_CONFIG.get('tso_species_map'), EXPERIMENT_CONFIG.get('tso_arch')))
if EXPERIMENT_CONFIG.get('barcode_location') == 'TSO':
    # 24plex detector 7-mers
    import override_24plex_illumina as _ox24
    if _hu:
        _ox24.BC_HUMAN = _hu
    if _ms:
        _ox24.BC_MOUSE = _ms
    # tso-bob detector 7-mers
    try:
        import override_bob_tso_7mer_illumina as _ot
        if _hu:
            _ot.BC_HUMAN = _hu
        if _ms:
            _ot.BC_MOUSE = _ms
    except Exception:
        pass
    # filter_funnel's generic, config-driven caller: anchor = 3' end of the TSO/R1
    # primer (immediately 5' of the barcode); barcodes + GGG offsets from the config.
    try:
        import filter_funnel_illumina as _ff
        if _hu:
            _ff.BC_HUMAN = _hu
        if _ms:
            _ff.BC_MOUSE = _ms
        _mp = _re.search(r'\[(?:TSO|R1) ([ACGT]{8,})\]',
                         EXPERIMENT_CONFIG.get('construct_description', ''))
        _anchor = (_mp.group(1) if _mp else
                   ('AGACGTGTGCTCTTCCGATCT' if _is24 else 'AAGCAGTGGTATCAACGCAG'))
        _ff.DETECT_ANCHOR = _anchor[-13:]
        # use the MORE-SENSITIVE backbone scan for tso-bob (non-24plex);
        # 24plex keeps the exact/measure_struct path (its barcode accuracy comes from
        # the full-molecule funnel, not filter_funnel's _bcd).
        _ff.DETECT_BACKBONE = None if _is24 else _anchor
        _ff.DETECT_BARCODES = {k: ('h' if v == 'human' else 'm')
                               for k, v in _pm.items() if v in ('human', 'mouse')}
        if EXPERIMENT_CONFIG.get('grun_offsets'):
            _ff._GRUN_OFFSET_BY_7MER = dict(EXPERIMENT_CONFIG['grun_offsets'])
        else:
            # Loud, not silent. Without this the module default {'GATATGG':8,'TCGTGAC':11}
            # stands and the G-run is measured INSIDE the UMI for any other bobcode pair.
            print("  WARNING: no 'grun_offsets' in the guide — G-run offsets fall back to "
                  f"module defaults {getattr(_ff, '_GRUN_OFFSET_BY_7MER', {})}; add the field "
                  "to the registry section or G-run/insert boundaries will be wrong.")
        # ILLUMINA: point the 24plex detector at this run's bobcodes/offsets AND switch it
        # to positional R2 mode. Without this the anchor 'CTCTTCCGATCT' is searched for in
        # a read where it is the sequencing primer (absent) -> 0% of reads get a barcode.
        try:
            import override_24plex_illumina as _ox24
            _layout = EXPERIMENT_CONFIG.get('read_layout') or (
                'r2_positional' if str(EXPERIMENT_CONFIG.get('platform', '')).startswith('illumina')
                else 'ont_anchor')
            _ff.READ_LAYOUT = _layout       # positional _bcd + orientation-correct _both_bc
            _unmeas = {x.strip() for x in
                       str(EXPERIMENT_CONFIG.get('unmeasurable', '') or '').split(',') if x.strip()}
            _polyt_ok = 'polyt' not in _unmeas
            _cfg = _ox24.configure(species_map=_pm,
                                   grun_offsets=EXPERIMENT_CONFIG.get('grun_offsets'),
                                   layout=_layout,
                                   polyt_measurable=_polyt_ok)
            print(f"  24plex detector: layout={_cfg['READ_LAYOUT']} codes={_cfg['CODES']} "
                  f"offsets={_cfg['GRUN_OFFSET_BY_7MER']} "
                  f"polyt_measurable={_cfg['POLYT_MEASURABLE']}")
        except Exception as _e:
            print(f"  (24plex detector configure skipped: {_e})")
        # full-molecule funnel's barcode->insert distance step: same per-barcode offsets
        try:
            import full_mol_funnel_illumina as _fmf
            if EXPERIMENT_CONFIG.get('grun_offsets'):
                _fmf.GRUN_OFFSET = dict(EXPERIMENT_CONFIG['grun_offsets'])
        except Exception:
            pass
        print(f"  detector (config-driven): anchor='{_ff.DETECT_ANCHOR}' "
              f"barcodes={_ff.DETECT_BARCODES} offsets={_ff._GRUN_OFFSET_BY_7MER}")
    except Exception as _e:
        print(f"  (filter_funnel config injection skipped: {_e})")

print(f"Sequence list loaded from: {os.path.basename(SEQ_LIST_PATH)}")
print(f"Experiment config loaded from: {os.path.basename(GUIDE_PATH)}")
print(f"  Barcodes: {EXPERIMENT_CONFIG['barcodes']}")
print(f"  Barcode location: {EXPERIMENT_CONFIG['barcode_location']}")
print(f"  Primary species map: {EXPERIMENT_CONFIG['primary_species_map']}")
print(f"  Construct: {EXPERIMENT_CONFIG['construct_description']}")
print(f"  FASTQ dir: {EXPERIMENT_CONFIG['fastq_dir']}")

# str.translate (C-level, whole string) instead of a per-character Python generator
# + str.join (~83M generator steps at 1M reads). Table maps EVERY byte to 'N'
# first, then overrides ACGTN -- equivalent to a `comp.get(b, 'N')` fallback.
# NB: this version upper()s its input (the override_24plex copy does NOT).
_RC_TAB = {i: 78 for i in range(256)}          # 78 = ord('N')
_RC_TAB.update(str.maketrans('ACGTN', 'TGCAN'))

def rc(seq):
    return seq.upper().translate(_RC_TAB)[::-1]


def _fuzzy_find_first(read, primer, max_mm):
    """Leftmost start where `primer` matches `read` with <= max_mm substitutions,
    or -1. EXACTLY equivalent to sliding every position and counting mismatches,
    but without the Python-level scan.

    PIGEONHOLE: split the primer into max_mm+1 blocks. A window differing in
    <= max_mm positions cannot have a mismatch in every block, so at least ONE
    block matches EXACTLY -- and str.find locates exact blocks in C. Only those
    few candidate starts are verified. This is a shortcut, NOT an approximation:
    every position a full scan would accept is still reachable as a candidate.

    A full scan slides an 18 bp primer over ~130 positions x 2 strands for EVERY
    read (~291M len() calls at 1M reads) and, when the primer is absent, fails at
    every position -- worst case on every single read.
    """
    L, n = len(primer), len(read)
    if n < L:
        return -1
    if max_mm <= 0:
        return read.find(primer)
    nb = max_mm + 1
    cands = set()
    for bi in range(nb):
        s, e = bi * L // nb, (bi + 1) * L // nb
        blk = primer[s:e]
        st = 0
        while True:
            i = read.find(blk, st)
            if i < 0:
                break
            si = i - s
            if 0 <= si <= n - L:
                cands.add(si)
            st = i + 1
    for si in sorted(cands):                 # leftmost, matching a full scan's order
        mm = 0
        for j in range(L):
            if read[si + j] != primer[j]:
                mm += 1
                if mm > max_mm:
                    break
        else:
            return si
    return -1

def read_fastq(path):
    with open(path) as f:
        while True:
            h = f.readline().rstrip()
            if not h: break
            seq = f.readline().rstrip().upper()
            f.readline(); f.readline()
            yield h, seq

# --- Derive all constants from sequence list ---
# Core sequences (looked up from sequence_list_illumina.txt, with fallbacks for safety)
TSO_PRIMER = GROUND_TRUTH.get('TSO_primer', 'GCAGTGGTATCAACGCAG')
TSO_PRIMER_RC = rc(TSO_PRIMER)
RT_PRIMER = GROUND_TRUTH.get('RT_HR', 'ATACTCGTGAC')
C28_LINKER = GROUND_TRUTH.get('C28_linker', 'ATCAGCTTCC')
POLYT_10 = GROUND_TRUTH.get('polyT_10', 'TTTTTTTTTT')
POLYA_10 = GROUND_TRUTH.get('polyA_10', 'AAAAAAAAAA')
C28_PCR = GROUND_TRUTH.get('C28-PCR-primer', 'ATACTCGTGACATCAGCTTCC')

# Full TSO variants (derived from TSO-seq*-variant entries, trimming terminal base)
TSO_VARIANTS = {}
for _gt_name in ('TSO-seq1-variant', 'TSO-seq2-variant'):
    _full = GROUND_TRUTH.get(_gt_name, '')
    if _full and len(_full) >= 24:
        _trimmed = _full[:-1]  # strip terminal T (25bp -> 24bp)
        _code = _trimmed[20:24]  # 4bp barcode at positions 20-23
        TSO_VARIANTS[_code] = _trimmed

# TSO barcode: 4bp at positions 20-23 after the 18bp TSO primer + 2bp AG linker
TSO_BARCODE_OFFSET = 20  # primer(18) + linker AG(2) = offset to barcode start

# BOB barcode detection: 4bp code between Bob-c prefix (20bp) and BOB_HR suffix (10bp)
BOB_PREFIX = GROUND_TRUTH.get('Bob-c_primer', 'ATGTAGTCCGTACGCATGTC')
BOB_SUFFIX = GROUND_TRUTH.get('BOB_HR', 'ACATGGTAGC')
BOB_PREFIX_RC = rc(BOB_PREFIX)
BOB_SUFFIX_RC = rc(BOB_SUFFIX)

# Derive expected BOB barcodes from ground truth BOB-C5-* entries
BOB_EXPECTED = {}
for gt_name, gt_seq in GROUND_TRUTH.items():
    if gt_name.startswith('BOB-C5-') and len(gt_seq) == 34:
        code = gt_seq[20:24]  # 4bp barcode at positions 20-23
        BOB_EXPECTED[code] = gt_name.replace('BOB-', '')

# Additional sequences to detect
READ1F = GROUND_TRUTH.get('Read1-f', 'GACTGGAGTTCAGACGTGT')

# Build SEQUENCE_REFERENCE from ground truth
# Each entry: (label, sequence) — searched fwd + RC in each read
# Skip entries with brackets (pattern descriptors like BOB_code)
SEQUENCE_REFERENCE = []
for gt_name, gt_seq in GROUND_TRUTH.items():
    if '[' in gt_seq or ']' in gt_seq:
        continue  # skip pattern descriptors
    label = f"{gt_name} ({len(gt_seq)}bp)"
    SEQUENCE_REFERENCE.append((label, gt_seq))

# Validate: print ground truth summary at import time
print(f"Ground truth loaded: {len(GROUND_TRUTH)} entries, "
      f"{len(SEQUENCE_REFERENCE)} searchable sequences, "
      f"BOB barcodes: {BOB_EXPECTED}")

def _config_bobcode_extras():
    """The run's own configured bobcodes (from its guide's species map) as scan
    `extra` records, enriched with reference-sheet names where the sequence is
    catalogued. This lets auto-detect name a bobcode that is genuinely in use but
    absent from the reference sheets."""
    try:
        import bobcode_reference_illumina as _bref
        _bref.load()
        _B = _bref.BARCODES
    except Exception:
        _B = {}
    smap = EXPERIMENT_CONFIG.get('primary_species_map') or {}
    extra = {}
    for seq, sp in smap.items():
        s = str(seq).upper()
        if not s.isalpha():
            continue
        rec = _B.get(s)
        extra[s] = dict(seq=s, name=(rec['name'] if rec else s),
                        id=(rec['id'] if rec else ''),
                        arch=(rec['arch'] if rec else 'config'),
                        species=sp, source='config')
    return extra


# ---- Empty-products metric (empty template-switch: polyA/polyT runs into the bobcode) ----
# "Empty products %" = reads where the RT-primer polyA/polyT abuts the bobcode with essentially
# no cDNA insert, as a percent of reads carrying >=1 barcode sequence (the headline denominator;
# the % of ALL reads is reported alongside). insert = (bobcode 7-mer -> first polyA/polyT run) -
# structural overhead (UMI + spacer + GGG, group-specific); a read is empty if insert < THRESHOLD.
# 24plex/V1 only — returns None for other architectures (rows render blank).
_EMPTY_ANCHOR = 'CTCTTCCGATCT'
_EMPTY_OVERHEAD = {'A': 11, 'B': 14, 'C': 17}   # UMI6 + HH2 + spacer(0/3/6) + GGG3 after the 7-mer
_EMPTY_V1_GROUP = {
    'AGTGGTG': 'A', 'AGCAGAC': 'A', 'TGCCGTT': 'A', 'GATATGG': 'A', 'TGGATGA': 'A',
    'GAATCTG': 'A', 'ACTCACA': 'A', 'CTCTTAC': 'A',
    'ACACCAT': 'B', 'GCGAATG': 'B', 'CTAACCT': 'B', 'CTGCGAT': 'B', 'CAAGCAA': 'B',
    'GTTGCTT': 'B', 'TCGTGAC': 'B', 'TGTCCAC': 'B',
    'TACAAGC': 'C', 'ATACTCC': 'C', 'CCAAGGA': 'C', 'AGCTCTA': 'C', 'CCTGTCT': 'C',
    'GAGTAGT': 'C', 'TAAGACG': 'C', 'GTCCATC': 'C',
}
_EMPTY_HOMO = _re.compile(r'(A{12,}|T{12,})')

def _empty_codes():
    """{code: overhead} for the artefact rows. Overhead = UMI + spacer + GGG = the run's
    per-code grun_offset + 3, for EVERY code the run config declares (any length); the
    hardcoded V1 group table is only the fallback for a config without offsets, so a run
    with non-V1 codes does not silently render these rows blank."""
    cfg = globals().get('EXPERIMENT_CONFIG') or {}
    offs = cfg.get('grun_offsets') or {}
    codes = {c: int(o) + 3 for c, o in offs.items()}
    if not codes:
        codes = {m: _EMPTY_OVERHEAD[g] for m, g in _EMPTY_V1_GROUP.items()}
    return codes


def _empty_detect(s):
    """(strand, anchor index, code) using the SAME barcode rule as every other row: anchor +
    declared code within 1 mismatch, longest code length first, ties rejected, forward
    strand then rc, first hit."""
    codes = _empty_codes()
    Ls = sorted({len(c) for c in codes}, reverse=True)
    for t in (s, _empty_rc(s)):
        j = t.find(_EMPTY_ANCHOR)
        while j >= 0:
            p0 = j + len(_EMPTY_ANCHOR)
            for L in Ls:
                seg = t[p0:p0 + L]
                if len(seg) < L:
                    continue
                ranked = sorted((sum(1 for x, y in zip(seg, c) if x != y), c) for c in codes if len(c) == L)
                if ranked and ranked[0][0] <= 1 and not (len(ranked) > 1 and ranked[1][0] == ranked[0][0]):
                    return (t, j, ranked[0][1])
            j = t.find(_EMPTY_ANCHOR, j + 1)
    return None



def _empty_rc(s):
    return s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def compute_empty_products(fastq_path, threshold=20):
    """Empty template-switch % for a 24plex/V1 barcode's reads. Returns
    {'threshold','n_total','n_bc','n_empty','pct_bc','pct_all'} or None if no V1 bobcodes
    are present (non-24plex chemistries)."""
    if not fastq_path or not os.path.exists(fastq_path):
        return None
    codes = _empty_codes()
    n_total = n_bc = n_empty = 0
    with open(fastq_path) as fh:
        for i, line in enumerate(fh):
            if i % 4 != 1:
                continue
            n_total += 1
            s = line.rstrip('\n')
            found = _empty_detect(s)
            if not found:
                continue
            n_bc += 1
            t, j, mer = found
            p = j + len(_EMPTY_ANCHOR) + len(mer)
            m = _EMPTY_HOMO.search(t, p, p + 500)
            if m and ((m.start() - p) - codes[mer]) < threshold:
                n_empty += 1
    if n_bc == 0:
        return None
    return {'threshold': threshold, 'n_total': n_total, 'n_bc': n_bc, 'n_empty': n_empty,
            'pct_bc': 100 * n_empty / n_bc, 'pct_all': 100 * n_empty / n_total if n_total else 0.0}


_EMPTY_ADAPTER = 'GATCGGAAGAGC'   # TruSeq adapter core (downstream boundary for poly(dT)-less)


def compute_polydt_less(fastq_path, threshold=20):
    """Poly(dT)-less artifact % for a 24plex/V1 barcode's reads: a bobcode is present, there is NO
    polyA/polyT run anywhere, and a TruSeq adapter begins within <threshold bp of the bobcode GGG
    (a short A-rich stub, no polydT/RT-primer arm — the ~170 bp band). Distinct from the empty
    template-switch products (which carry the polyA). Returns
    {'threshold','n_total','n_bc','n_art','pct_all','pct_bc'} or None (non-24plex chemistries)."""
    if not fastq_path or not os.path.exists(fastq_path):
        return None
    codes = _empty_codes()
    n_total = n_bc = n_art = 0
    with open(fastq_path) as fh:
        for i, line in enumerate(fh):
            if i % 4 != 1:
                continue
            n_total += 1
            s = line.rstrip('\n')
            found = _empty_detect(s)
            if not found:
                continue
            n_bc += 1
            t, j, mer = found
            if _EMPTY_HOMO.search(t):            # any polyA/polyT run -> NOT poly(dT)-less
                continue
            ins_start = j + len(_EMPTY_ANCHOR) + len(mer) + codes[mer]
            k = t.find(_EMPTY_ADAPTER, ins_start, ins_start + threshold + len(_EMPTY_ADAPTER))
            if k >= 0 and (k - ins_start) < threshold:
                n_art += 1
    if n_bc == 0:
        return None
    return {'threshold': threshold, 'n_total': n_total, 'n_bc': n_bc, 'n_art': n_art,
            'pct_all': 100 * n_art / n_total if n_total else 0.0, 'pct_bc': 100 * n_art / n_bc}


def analyze_barcode(bc_dir, fastq_file, output_dir=None):
    fastq_path = os.path.join(bc_dir, fastq_file)
    bc_name = os.path.basename(bc_dir.rstrip('/'))
    if output_dir is None:
        output_dir = bc_dir

    print(f"\n{'='*60}")
    print(f"RT+TSO Analysis: {bc_name}")
    print(f"{'='*60}")

    reads = []
    read_ids = []  # FASTQ header IDs (first word after @)
    for header, seq in read_fastq(fastq_path):
        reads.append(seq)
        # Extract read ID from header (first whitespace-delimited token, strip @)
        rid = header.lstrip('@').split()[0]
        read_ids.append(rid)
    total = len(reads)
    print(f"Total reads: {total:,}")

    # --- Auto-detect the bobcodes in use, straight from the reads (config-free) ---
    # Scans a read sample for the anchor signatures (24plex TruSeq / tso-bob
    # SMARTer / BOB-J Nextera) and matches the downstream barcode window to the
    # reference registry, so we know WHICH named bobcodes this sample used and in
    # what proportion -- reported at the top of the cross-barcode PDF. Cheap
    # (sampled str.find + O(1) dict lookups) and never fatal.
    detected_bobcodes = None
    try:
        import scan_bobcodes_illumina as _scan
        detected_bobcodes = _scan.scan(reads, n=30000, extra=_config_bobcode_extras())
    except Exception as _e:
        print(f"  (bobcode auto-detect skipped: {_e})")
    if detected_bobcodes and detected_bobcodes.get('architecture'):
        _db = detected_bobcodes
        _bl = ", ".join(
            f"{b['id'] or b['name']} {b['seq']}"
            + (f"/{b['species']}" if b['species'] else "")
            + f" ({b['frac']*100:.0f}%)" for b in _db['barcodes'][:4])
        print(f"  auto-detected: {_db['architecture']} | {_bl}")

    # --- BOB barcode extraction ---
    bob_barcode_counts = Counter()
    bob_barcode_reads = defaultdict(list)  # barcode -> list of read indices
    read_bob_barcode = {}  # read_id -> bob_barcode

    # --- TSO barcode extraction ---
    tso_barcode_counts = Counter()
    tso_barcode_reads = defaultdict(list)  # barcode -> list of read indices
    reads_with_tso = 0
    reads_with_rt = 0
    reads_with_c28_pcr = 0     # C28-PCR-primer: RT_HR + linker (21bp)
    reads_with_rt_c28_6t = 0   # RT_HR + C28 linker + 6T
    reads_with_rt_c28_6n = 0   # RT_HR + C28 linker + T + 6N (not all T)
    reads_with_polyt = 0
    rt_inferred = Counter()    # 4-category RT primer type: PolydT-5P, PolydT_5P_C28, C28-RT-N9, RT-N9
    reads_complete = 0  # has RT + polyT + TSO
    reads_with_read1f = 0
    reads_with_tso_and_read1f = 0  # has both TSO and Read1f
    reads_with_tso_or_read1f = 0   # has either TSO or Read1f
    reads_without_tso_or_read1f = 0  # has neither TSO nor Read1f
    reads_with_bob_c = 0     # BOB-c prefix (20bp)
    reads_with_bob_hr = 0    # BOB_HR suffix (10bp)
    reads_with_c28_both_ends = 0  # C28 linker on BOTH ends of read (2xRT self-primed)
    LIGATED_SEQ = GROUND_TRUTH.get('BOB_HR-RT_HR', 'ACATGGTAGCATACTCGTGAC')  # BOB_HR + RT_HR ligation junction
    LIGATED_SEQ_RC = rc('ACATGGTAGCATACTCGTGAC')
    reads_with_ligated = 0
    reads_multi_bob_c = 0    # 2+ BOB-c primers
    reads_multi_ligated = 0  # 2+ ligation junctions

    # For cDNA extraction
    insert_lengths = []
    inserts_fasta = []
    read_tso_barcode = {}  # read_id -> tso_barcode

    # Read structure classification
    RT_PRIMER_RC = rc(RT_PRIMER)
    read1f_rc = rc(READ1F)
    structure_counts = Counter()  # e.g. "RT > polydt > insert > TSO" -> count
    structure_insert_lengths = defaultdict(list)  # struct_str -> list of insert lengths
    read_structures = []  # per-read structure string
    reads_multi_rt = 0   # reads with 2+ distinct RT primer hits
    reads_multi_tso = 0  # reads with 2+ distinct TSO primer hits
    reads_multi_read1f = 0  # reads with 2+ distinct Read1f hits

    # RT subtype detection sequences
    RT_C28_PREFIX = RT_PRIMER + C28_LINKER + 'T'  # RT_HR + C28 linker + T (22bp)
    RT_C28_PREFIX_RC = rc(RT_C28_PREFIX)
    RT_C28_6T = RT_C28_PREFIX + 'TTTTT'  # + 5 more T (total 6T after linker+T)
    RT_C28_6T_RC = rc(RT_C28_6T)
    POLYDT_NOC28 = RT_PRIMER + 'TTTTT'      # 16bp: RT_HR + 5T for PolydT-5P detection
    POLYDT_NOC28_RC = rc(POLYDT_NOC28)

    # PRECOMPILED per-read patterns (loop-invariant). Building them INSIDE the
    # per-read loop as re.finditer(re.escape(X), read) would run re.escape and probe
    # the re cache on EVERY read: ~14M re.escape calls and 12M finditer calls at
    # 1M reads. They depend only on module constants and the two RCs above.
    _rx_bob_fwd   = re.compile(re.escape(BOB_PREFIX) + r'(.{4})' + re.escape(BOB_SUFFIX))
    _rx_bob_rc    = re.compile(re.escape(BOB_SUFFIX_RC) + r'(.{4})' + re.escape(BOB_PREFIX_RC))
    _rx_rt        = re.compile(re.escape(RT_PRIMER))
    _rx_rt_rc     = re.compile(re.escape(RT_PRIMER_RC))
    _rx_polyt     = re.compile(r'T{10,}')
    _rx_polya     = re.compile(r'A{10,}')
    _rx_polya_end = re.compile(r'A{10,}$')
    _rx_tso       = re.compile(re.escape(TSO_PRIMER))
    _rx_tso_rc    = re.compile(re.escape(TSO_PRIMER_RC))
    _rx_r1f       = re.compile(re.escape(READ1F))
    _rx_r1f_rc    = re.compile(re.escape(read1f_rc))
    _rx_bobsuf    = re.compile(re.escape(BOB_SUFFIX))
    _rx_bobsuf_rc = re.compile(re.escape(BOB_SUFFIX_RC))
    _rx_bobpre    = re.compile(re.escape(BOB_PREFIX))
    _rx_bobpre_rc = re.compile(re.escape(BOB_PREFIX_RC))

    for i, read in enumerate(reads):
        has_rt = RT_PRIMER in read or RT_PRIMER_RC in read
        has_polyt = POLYT_10 in read or POLYA_10 in read
        has_read1f = READ1F in read or read1f_rc in read
        has_bob_c = BOB_PREFIX in read or BOB_PREFIX_RC in read
        has_bob_hr = BOB_SUFFIX in read or BOB_SUFFIX_RC in read
        has_ligated = LIGATED_SEQ in read or LIGATED_SEQ_RC in read
        if has_rt: reads_with_rt += 1
        if has_polyt: reads_with_polyt += 1
        if has_read1f: reads_with_read1f += 1
        if has_bob_c: reads_with_bob_c += 1
        if has_bob_hr: reads_with_bob_hr += 1
        if has_ligated: reads_with_ligated += 1
        # Count 2+ occurrences of BOB-c and ligation junction
        bob_c_fwd = read.count(BOB_PREFIX) + read.count(BOB_PREFIX_RC)
        if bob_c_fwd >= 2: reads_multi_bob_c += 1
        lig_fwd = read.count(LIGATED_SEQ) + read.count(LIGATED_SEQ_RC)
        if lig_fwd >= 2: reads_multi_ligated += 1

        # C28-PCR-primer: RT_HR + linker (21bp)
        c28_pcr = C28_PCR
        if c28_pcr in read or rc(c28_pcr) in read:
            reads_with_c28_pcr += 1
        # C28 linker on BOTH ends — indicator of 2xRT self-primed products.
        # Search the 10bp C28 linker fwd in the 5' region and its RC in the 3' region.
        _c28_linker = 'ATCAGCTTCC'
        _c28_linker_rc = rc(_c28_linker)
        _head = read[:80]; _tail = read[-80:]
        _c28_fwd_5p = _c28_linker in _head
        _c28_rc_3p = _c28_linker_rc in _tail
        _c28_fwd_3p = _c28_linker in _tail  # reverse-oriented read
        _c28_rc_5p = _c28_linker_rc in _head
        # Either orientation counts as "both ends"
        if (_c28_fwd_5p and _c28_rc_3p) or (_c28_fwd_3p and _c28_rc_5p):
            reads_with_c28_both_ends += 1
        # RT subtype: C28_6T (has 6T after C28 prefix)
        if RT_C28_6T in read or RT_C28_6T_RC in read:
            reads_with_rt_c28_6t += 1
        # RT subtype: C28_6N (has C28 prefix + 6 bases that are NOT all T/A)
        elif RT_C28_PREFIX in read or RT_C28_PREFIX_RC in read:
            # Check forward
            found_6n = False
            for search_seq, is_rc in [(RT_C28_PREFIX, False), (RT_C28_PREFIX_RC, True)]:
                idx = read.find(search_seq)
                if idx >= 0:
                    if is_rc:
                        # For RC, the 6N would be BEFORE the prefix
                        start = idx - 6
                        if start >= 0:
                            sixmer = read[start:idx]
                            # In RC context, all-T in original = all-A in RC
                            if sixmer != 'AAAAAA':
                                found_6n = True
                    else:
                        start = idx + len(RT_C28_PREFIX)
                        if start + 6 <= len(read):
                            sixmer = read[start:start+6]
                            if sixmer != 'TTTTTT':
                                found_6n = True
            if found_6n:
                reads_with_rt_c28_6n += 1

        # Inferred RT primer type (4-category decision tree)
        if has_rt:
            has_c28_in_read = c28_pcr in read or rc(c28_pcr) in read
            if has_c28_in_read:
                # Check for >=10 Ts after C28 region (positional)
                polyt_after_c28 = False
                c28_pcr_rc = rc(c28_pcr)
                fwd_idx = read.find(c28_pcr)
                if fwd_idx >= 0:
                    after = read[fwd_idx + len(c28_pcr):]
                    if _rx_polyt.match(after):
                        polyt_after_c28 = True
                else:
                    rc_idx = read.find(c28_pcr_rc)
                    if rc_idx >= 0:
                        # In RC orientation, poly-T becomes poly-A upstream
                        before = read[:rc_idx]
                        if _rx_polya_end.search(before):
                            polyt_after_c28 = True
                if polyt_after_c28:
                    rt_inferred['PolydT_5P_C28'] += 1
                else:
                    rt_inferred['C28-RT-N9'] += 1
            else:
                # No C28: check for poly-T after RT_HR (fast substring check)
                if POLYDT_NOC28 in read or POLYDT_NOC28_RC in read:
                    rt_inferred['PolydT-5P'] += 1
                else:
                    rt_inferred['RT-N9'] += 1

        # Find TSO primer (18bp) with fuzzy matching (1 mismatch allowed)
        tso_bc = None
        tso_pos = None
        tso_max_mm = len(TSO_PRIMER) // 10  # 1 mismatch for 18bp

        # Forward: find TSO primer, then extract 4bp barcode at offset 20.
        # NB: only the FIRST hit is considered -- if the barcode window runs off the
        # read end, later positions are not tried, so tso_bc stays None
        # (and the RC branch below still gets its turn).
        si = _fuzzy_find_first(read, TSO_PRIMER, tso_max_mm)
        if si >= 0:
            bc_start = si + TSO_BARCODE_OFFSET
            if bc_start + 4 <= len(read):
                tso_bc = read[bc_start:bc_start + 4]
                tso_pos = si

        # RC: find TSO_PRIMER_RC, extract barcode upstream of it
        if tso_bc is None:
            si = _fuzzy_find_first(read, TSO_PRIMER_RC, tso_max_mm)
            if si >= 0:
                # RC primer found; barcode is 4bp before offset upstream
                bc_end = si - (TSO_BARCODE_OFFSET - len(TSO_PRIMER))
                bc_start = bc_end - 4
                if bc_start >= 0:
                    tso_bc = rc(read[bc_start:bc_end])
                    tso_pos = si

        if tso_bc:
            reads_with_tso += 1
            tso_barcode_counts[tso_bc] += 1
            tso_barcode_reads[tso_bc].append(i)
            read_tso_barcode[read_ids[i]] = tso_bc

        if has_rt and has_polyt and tso_bc:
            reads_complete += 1
        has_tso = tso_bc is not None
        if has_tso and has_read1f:
            reads_with_tso_and_read1f += 1
        if has_tso or has_read1f:
            reads_with_tso_or_read1f += 1
        if not has_tso and not has_read1f:
            reads_without_tso_or_read1f += 1

        # --- BOB barcode extraction ---
        # Forward: BOB_PREFIX + [4bp] + BOB_SUFFIX
        bob_bc = None
        m = _rx_bob_fwd.search(read)
        if m:
            bob_bc = m.group(1)
        else:
            # RC: RC(BOB_SUFFIX) + [4bp_RC] + RC(BOB_PREFIX)
            m = _rx_bob_rc.search(read)
            if m:
                bob_bc = rc(m.group(1))
        if bob_bc:
            bob_barcode_counts[bob_bc] += 1
            bob_barcode_reads[bob_bc].append(i)
            read_bob_barcode[read_ids[i]] = bob_bc

        # --- Read structure classification ---
        # Find positions (start, end) of all components to determine order
        # and infer insert (cDNA) from gaps between components
        components = []  # list of (type, start, end)

        # RT primer (find all occurrences, fwd and RC)
        for mt in _rx_rt.finditer(read):
            components.append(('RT_HR', mt.start(), mt.end()))
        for mt in _rx_rt_rc.finditer(read):
            components.append(('RT_HR', mt.start(), mt.end()))

        # PolyT (≥10 consecutive T or A)
        for mt in _rx_polyt.finditer(read):
            components.append(('polydt', mt.start(), mt.end()))
        for mt in _rx_polya.finditer(read):
            components.append(('polydt', mt.start(), mt.end()))

        # TSO (using primer, not just barcode context)
        for mt in _rx_tso.finditer(read):
            components.append(('TSO', mt.start(), mt.end()))
        for mt in _rx_tso_rc.finditer(read):
            components.append(('TSO', mt.start(), mt.end()))

        # Read1f (Illumina adapter)
        for mt in _rx_r1f.finditer(read):
            components.append(('Read1f', mt.start(), mt.end()))
        for mt in _rx_r1f_rc.finditer(read):
            components.append(('Read1f', mt.start(), mt.end()))

        # BOB_HR (10bp suffix)
        for mt in _rx_bobsuf.finditer(read):
            components.append(('BOB_HR', mt.start(), mt.end()))
        for mt in _rx_bobsuf_rc.finditer(read):
            components.append(('BOB_HR', mt.start(), mt.end()))

        # BOB-c (20bp prefix)
        for mt in _rx_bobpre.finditer(read):
            components.append(('bob-c', mt.start(), mt.end()))
        for mt in _rx_bobpre_rc.finditer(read):
            components.append(('bob-c', mt.start(), mt.end()))

        # Sort by position and deduplicate adjacent same-type
        components.sort(key=lambda x: x[1])
        deduped = []
        for comp_type, start, end in components:
            if not deduped or deduped[-1][0] != comp_type:
                deduped.append((comp_type, start, end))
            else:
                # Extend end if same type and overlapping/adjacent
                deduped[-1] = (comp_type, deduped[-1][1], max(deduped[-1][2], end))

        # Count reads with multiple RT, TSO, or Read1f components (after dedup)
        rt_count_deduped = sum(1 for c in deduped if c[0] == 'RT_HR')
        tso_count_deduped = sum(1 for c in deduped if c[0] == 'TSO')
        read1f_count_deduped = sum(1 for c in deduped if c[0] == 'Read1f')
        if rt_count_deduped >= 2:
            reads_multi_rt += 1
        if tso_count_deduped >= 2:
            reads_multi_tso += 1
        if read1f_count_deduped >= 2:
            reads_multi_read1f += 1

        # Insert "insert" labels for gaps ≥ 30bp between components
        # and measure the insert (cDNA) length from the gap
        MIN_INSERT_GAP = 30
        with_inserts = []
        insert_gap_len = 0  # largest gap = primary insert
        for j, (comp_type, start, end) in enumerate(deduped):
            if j > 0:
                prev_end = deduped[j-1][2]
                gap = start - prev_end
                if gap >= MIN_INSERT_GAP:
                    with_inserts.append('insert')
                    if gap > insert_gap_len:
                        insert_gap_len = gap
            with_inserts.append(comp_type)

        # Remove polydt from structure display (keep bob-c)
        filtered_inserts = [c for c in with_inserts if c not in ('polydt',)]
        # Collapse adjacent 'insert' entries that may result from removal
        collapsed = []
        for c in filtered_inserts:
            if c == 'insert' and collapsed and collapsed[-1] == 'insert':
                continue
            collapsed.append(c)
        raw_struct = ' > '.join(collapsed) if collapsed else 'none'

        # Normalize orientation: prefer BOB_HR first when present, else alphabetically first
        rev_struct = ' > '.join(reversed(collapsed)) if collapsed else 'none'
        if 'BOB_HR' in collapsed:
            # Pick whichever puts BOB_HR first
            if raw_struct.startswith('BOB_HR'):
                struct_str = raw_struct
            elif rev_struct.startswith('BOB_HR'):
                struct_str = rev_struct
            else:
                struct_str = min(raw_struct, rev_struct)
        else:
            struct_str = min(raw_struct, rev_struct)

        structure_counts[struct_str] += 1
        read_structures.append(struct_str)

        # Track insert length per structure type
        if insert_gap_len >= MIN_INSERT_GAP:
            structure_insert_lengths[struct_str].append(insert_gap_len)

        # --- cDNA extraction ---
        # Pattern A: ...polyT...[cDNA]...[TSO_RC]
        # Pattern B: [TSO_fwd]...[cDNA]...polyA...
        insert_seq = None

        polyt_pos = read.find(POLYT_10)
        polya_pos = read.find(POLYA_10)

        # Find any TSO anchor (AGTACA or AGGGCA or TSO primer)
        for tso_full in TSO_VARIANTS.values():
            tso_rc_seq = rc(tso_full)

            # Pattern: ...polyT...[cDNA]...[TSO_full_RC]
            if polyt_pos != -1:
                anchor_pos = read.find(tso_rc_seq, polyt_pos)
                if anchor_pos != -1 and anchor_pos > polyt_pos:
                    end_polyt = polyt_pos
                    while end_polyt < len(read) and read[end_polyt] == 'T':
                        end_polyt += 1
                    if anchor_pos > end_polyt and (anchor_pos - end_polyt) >= 30:
                        insert_seq = read[end_polyt:anchor_pos]
                        break

            # Pattern: [TSO_full_fwd]...[cDNA]...polyA
            if polya_pos != -1:
                anchor_pos = read.find(tso_full)
                if anchor_pos != -1 and anchor_pos < polya_pos:
                    start = anchor_pos + len(tso_full)
                    pa_start = polya_pos
                    while pa_start > start and read[pa_start-1] == 'A':
                        pa_start -= 1
                    if (pa_start - start) >= 30:
                        insert_seq = read[start:pa_start]
                        break

        # Fallback: use TSO primer as anchor
        if insert_seq is None:
            tso_primer_rc_pos = read.find(TSO_PRIMER_RC)
            if polyt_pos != -1 and tso_primer_rc_pos != -1 and tso_primer_rc_pos > polyt_pos:
                end_polyt = polyt_pos
                while end_polyt < len(read) and read[end_polyt] == 'T':
                    end_polyt += 1
                # TSO RC region starts ~6bp before the primer RC (barcode+AG)
                anchor = max(tso_primer_rc_pos - 6, end_polyt)
                if anchor > end_polyt and (anchor - end_polyt) >= 30:
                    insert_seq = read[end_polyt:anchor]

            tso_fwd_pos = read.find(TSO_PRIMER)
            if insert_seq is None and polya_pos != -1 and tso_fwd_pos != -1 and tso_fwd_pos < polya_pos:
                start = tso_fwd_pos + len(TSO_PRIMER) + 6  # skip primer + AG + 4bp barcode
                pa_start = polya_pos
                while pa_start > start and read[pa_start-1] == 'A':
                    pa_start -= 1
                if (pa_start - start) >= 30:
                    insert_seq = read[start:pa_start]

        if insert_seq and len(insert_seq) >= 30:
            insert_lengths.append(len(insert_seq))
            inserts_fasta.append((read_ids[i], insert_seq))

    print(f"\nStructure detection:")
    print(f"  RT primer (ATACTCGTGAC): {reads_with_rt:,} ({reads_with_rt/total*100:.1f}%)")
    print(f"  PolyT/A (≥10bp):         {reads_with_polyt:,} ({reads_with_polyt/total*100:.1f}%)")
    print(f"  TSO detected:            {reads_with_tso:,} ({reads_with_tso/total*100:.1f}%)")
    print(f"  Read1f detected:         {reads_with_read1f:,} ({reads_with_read1f/total*100:.1f}%)")
    print(f"  Complete (RT+polyT+TSO):  {reads_complete:,} ({reads_complete/total*100:.1f}%)")

    print(f"\nRead structure classification:")
    import statistics as _stats
    for struct, count in structure_counts.most_common(10):
        lens = structure_insert_lengths.get(struct, [])
        if lens:
            med = _stats.median(lens)
            print(f"  {struct}: {count:,} ({count/total*100:.1f}%) — median insert {med:.0f}bp")
        else:
            print(f"  {struct}: {count:,} ({count/total*100:.1f}%)")

    # Biological grouping of structures
    bio_groups = Counter()
    bio_group_insert_lengths = defaultdict(list)  # bio_group -> list of insert lengths
    for idx_s, struct in enumerate(read_structures):
        # Count components (ignoring 'insert' which is inferred)
        parts = [p for p in struct.split(' > ') if p != 'insert']
        n_rt = parts.count('RT_HR')
        n_tso = parts.count('TSO')
        has_pd = 'polydt' in parts

        if n_rt >= 1 and n_tso >= 1:
            grp = 'RT + TSO (intended)'
        elif n_rt >= 2:
            grp = 'RT x2 (self-primed)'
        elif n_rt == 1 and has_pd:
            grp = 'RT + polydt only (truncated)'
        elif n_rt == 1:
            grp = 'RT only (truncated)'
        elif n_tso >= 1:
            grp = 'TSO only (no RT)'
        elif has_pd:
            grp = 'polydt only'
        else:
            grp = 'none detected'
        bio_groups[grp] += 1

        # Track insert lengths for this bio group
        lens = structure_insert_lengths.get(struct, [])
        # Each struct's insert lengths list corresponds to reads with that struct in order
        # We need per-read tracking; use structure_insert_lengths which was filled per read
        # Actually, structure_insert_lengths[struct] has one entry per read with that struct that had an insert
        # We can't directly index, so instead compute from the struct's stats

    # Compute bio group insert stats by aggregating from structure-level data
    for struct, lens in structure_insert_lengths.items():
        parts = [p for p in struct.split(' > ') if p != 'insert']
        n_rt = parts.count('RT_HR')
        n_tso = parts.count('TSO')
        has_pd = 'polydt' in parts
        if n_rt >= 1 and n_tso >= 1:
            grp = 'RT + TSO (intended)'
        elif n_rt >= 2:
            grp = 'RT x2 (self-primed)'
        elif n_rt == 1 and has_pd:
            grp = 'RT + polydt only (truncated)'
        elif n_rt == 1:
            grp = 'RT only (truncated)'
        elif n_tso >= 1:
            grp = 'TSO only (no RT)'
        elif has_pd:
            grp = 'polydt only'
        else:
            grp = 'none detected'
        bio_group_insert_lengths[grp].extend(lens)

    print(f"\nBiological grouping:")
    for grp, count in bio_groups.most_common():
        lens = bio_group_insert_lengths.get(grp, [])
        if lens:
            med = _stats.median(lens)
            print(f"  {grp}: {count:,} ({count/total*100:.1f}%) — median insert {med:.0f}bp")
        else:
            print(f"  {grp}: {count:,} ({count/total*100:.1f}%)")

    print(f"\nTSO 4bp barcodes:")
    for bc, count in tso_barcode_counts.most_common():
        print(f"  {bc}: {count:,} ({count/total*100:.1f}%)")

    print(f"\ncDNA inserts extracted: {len(insert_lengths):,}")

    # Save inserts FASTA
    fasta_path = os.path.join(output_dir, "cdna_inserts.fasta")
    with open(fasta_path, 'w') as f:
        for rid, seq in inserts_fasta:
            f.write(f">{rid}\n{seq}\n")

    # Length stats
    if insert_lengths:
        import statistics
        len_mean = statistics.mean(insert_lengths)
        len_median = statistics.median(insert_lengths)
        len_min = min(insert_lengths)
        len_max = max(insert_lengths)
    else:
        len_mean = len_median = len_min = len_max = 0

    print(f"  Mean length: {len_mean:.0f}bp, Median: {len_median:.0f}bp")

    # --- Duplicate analysis ---
    insert_seqs = [seq for _, seq in inserts_fasta]
    n_inserts = len(insert_seqs)

    # Exact insert duplicates
    insert_exact_counts = Counter(insert_seqs)
    exact_dup_reads = sum(c for c in insert_exact_counts.values() if c > 1)
    exact_dup_groups = sum(1 for c in insert_exact_counts.values() if c > 1)

    # Start+end fingerprint duplicates (20bp each side)
    insert_fps = Counter()
    for seq in insert_seqs:
        if len(seq) >= 40:
            fp = seq[:20] + '|' + seq[-20:]
        else:
            fp = seq
        insert_fps[fp] += 1
    startend_dup_reads = sum(c for c in insert_fps.values() if c > 1)
    startend_dup_groups = sum(1 for c in insert_fps.values() if c > 1)

    print(f"\nDuplicate inserts:")
    print(f"  Exact duplicates: {exact_dup_reads} reads in {exact_dup_groups} groups")
    print(f"  Start+end(20bp) duplicates: {startend_dup_reads} reads in {startend_dup_groups} groups")

    # Unique insert count (after dedup by start+end fingerprint)
    unique_inserts = len(insert_fps)

    # --- Sequence reference detection: exact match (fwd + RC) ---
    # ILLUMINA: SKIPPED. The "Sequence Reference Detection" PDF section is gated
    # `if not IS_ILLUMINA` (it searches ONT-era full-molecule primer/adapter
    # sequences that a 150 bp trimmed mate cannot contain), so on Illumina NOTHING
    # consumes seq_ref_counts -- and the scan would cost 114 reference sequences
    # x every read = ~107M substring searches, ~27% of the whole run at 1M reads.
    # Left EMPTY rather than zero-filled: an explicit 0 would read as "searched and
    # not found", which is a different claim from "not searched".
    seq_ref_counts = {}
    if not IS_ILLUMINA:
        for label, seq in SEQUENCE_REFERENCE:
            seq_rc = rc(seq)
            count = sum(1 for read in reads if seq in read or seq_rc in read)
            seq_ref_counts[label] = count
        print(f"\nSequence reference detection:")
        for label, count in seq_ref_counts.items():
            print(f"  {label}: {count} ({count/total*100:.1f}%)")

    # --- Nucleotide frequency flanking primers ---
    def compute_flank_freq(reads, primer_seq, flank=5):
        """Compute per-position nucleotide frequency around a primer.
        Returns dict with 'upstream' and 'downstream' keys, each a list of dicts.
        Positions are relative to the primer: upstream[-5..-1], downstream[+1..+5].
        All orientations normalized to forward primer direction.
        """
        primer_rc = rc(primer_seq)
        plen = len(primer_seq)
        # upstream[i] = Counter for position -(flank-i) relative to primer start
        upstream = [Counter() for _ in range(flank)]
        downstream = [Counter() for _ in range(flank)]
        n_found = 0
        for read in reads:
            # Try forward
            idx = read.find(primer_seq)
            if idx >= 0:
                n_found += 1
                for j in range(flank):
                    pos = idx - flank + j
                    if 0 <= pos < len(read):
                        upstream[j][read[pos]] += 1
                for j in range(flank):
                    pos = idx + plen + j
                    if pos < len(read):
                        downstream[j][read[pos]] += 1
                continue
            # Try RC
            idx = read.find(primer_rc)
            if idx >= 0:
                n_found += 1
                # RC found: the region after primer_rc in read = upstream of primer in fwd
                # We need to RC the flanking context to get fwd orientation
                # In fwd orientation: [upstream][primer][downstream]
                # In read (RC): [RC(downstream)][RC(primer)][RC(upstream)]
                # So: downstream of primer (fwd) = RC of bases BEFORE primer_rc in read
                # upstream of primer (fwd) = RC of bases AFTER primer_rc in read
                for j in range(flank):
                    pos = idx + plen + (flank - 1 - j)
                    if pos < len(read):
                        base = {'A':'T','T':'A','G':'C','C':'G'}.get(read[pos], 'N')
                        upstream[j][base] += 1
                for j in range(flank):
                    pos = idx - 1 - j
                    if 0 <= pos < len(read):
                        base = {'A':'T','T':'A','G':'C','C':'G'}.get(read[pos], 'N')
                        downstream[j][base] += 1
        return {
            'n_reads': n_found,
            'upstream': [dict(c) for c in upstream],
            'downstream': [dict(c) for c in downstream],
        }

    flank_size = 10
    flank_data = {}
    flank_primers = [
        ('RT_HR', RT_PRIMER),
        ('TSO primer', TSO_PRIMER),
        ('Read1-f', READ1F),
        ('adaptor-1-short', GROUND_TRUTH.get('adaptor-1-short', 'GCTCTTCCGATC')),
        ('TSO-seq1-variant', GROUND_TRUTH.get('TSO-seq1-variant', 'GCAGTGGTATCAACGCAGAGTACAT')),
    ]
    for primer_name, primer_seq in flank_primers:
        flank_data[primer_name] = compute_flank_freq(reads, primer_seq, flank=flank_size)
        print(f"\nFlanking frequency for {primer_name} ({flank_data[primer_name]['n_reads']} reads):")
        for j in range(flank_size):
            up = flank_data[primer_name]['upstream'][j]
            total_up = sum(up.values())
            if total_up > 0:
                pcts = {b: f"{up.get(b,0)/total_up*100:.0f}%" for b in 'ATGC'}
                print(f"  pos -{flank_size-j}: {pcts}")
        for j in range(flank_size):
            dn = flank_data[primer_name]['downstream'][j]
            total_dn = sum(dn.values())
            if total_dn > 0:
                pcts = {b: f"{dn.get(b,0)/total_dn*100:.0f}%" for b in 'ATGC'}
                print(f"  pos +{j+1}: {pcts}")

    # Save results
    results = {
        'barcode_name': os.path.basename(bc_dir.rstrip('/')),
        'fastq_file': fastq_file,
        # ILLUMINA: SKIPPED (None, not 0 -- see seq_ref_counts note above). Both rows
        # ('Empty products % (of all reads)', 'Poly(dT)-less artifacts % (of all
        # reads)') are in _DROP_ROWS for Illumina: they key on the ONT full-molecule
        # read, and a 150 bp mate sees one end only. Each re-scans the whole FASTQ.
        'empty_products': None if IS_ILLUMINA else compute_empty_products(fastq_path),
        'polydt_less': None if IS_ILLUMINA else compute_polydt_less(fastq_path),
        'construct_type': 'rt_tso_only',
        'total_reads': total,
        'reads_with_rt': reads_with_rt,
        'reads_with_c28_pcr': reads_with_c28_pcr,
        'reads_with_c28_both_ends': reads_with_c28_both_ends,
        'reads_with_rt_c28_6t': reads_with_rt_c28_6t,
        'reads_with_rt_c28_6n': reads_with_rt_c28_6n,
        'reads_with_polyt': reads_with_polyt,
        'rt_inferred': dict(rt_inferred),
        'reads_with_tso': reads_with_tso,
        'reads_with_read1f': reads_with_read1f,
        'reads_with_tso_and_read1f': reads_with_tso_and_read1f,
        'reads_with_tso_or_read1f': reads_with_tso_or_read1f,
        'reads_without_tso_or_read1f': reads_without_tso_or_read1f,
        'reads_with_bob_c': reads_with_bob_c,
        'reads_with_bob_hr': reads_with_bob_hr,
        'reads_with_ligated': reads_with_ligated,
        'reads_multi_bob_c': reads_multi_bob_c,
        'reads_multi_ligated': reads_multi_ligated,
        'reads_complete': reads_complete,
        'reads_multi_rt': reads_multi_rt,
        'reads_multi_tso': reads_multi_tso,
        'reads_multi_read1f': reads_multi_read1f,
        'tso_barcode_counts': dict(tso_barcode_counts),
        'bob_barcode_counts': dict(bob_barcode_counts),
        'read_bob_barcodes': dict(read_bob_barcode),
        'insert_count': n_inserts,
        'len_mean': len_mean,
        'len_median': len_median,
        'len_min': len_min,
        'len_max': len_max,
        'read_tso_barcodes': dict(read_tso_barcode),
        'structure_counts': dict(structure_counts),
        'structure_insert_stats': {
            s: {'n': len(lens), 'median': round(_stats.median(lens)), 'mean': round(_stats.mean(lens))}
            for s, lens in structure_insert_lengths.items() if lens
        },
        'bio_groups': dict(bio_groups),
        'bio_group_insert_stats': {
            g: {'n': len(lens), 'median': round(_stats.median(lens)), 'mean': round(_stats.mean(lens))}
            for g, lens in bio_group_insert_lengths.items() if lens
        },
        'duplicates': {
            'exact_dup_reads': exact_dup_reads,
            'exact_dup_groups': exact_dup_groups,
            'startend_dup_reads': startend_dup_reads,
            'startend_dup_groups': startend_dup_groups,
            'unique_inserts': unique_inserts,
        },
        'seq_ref_counts': seq_ref_counts,
        'flank_freq': flank_data,
        'detected_bobcodes': detected_bobcodes,
    }

    json_path = os.path.join(output_dir, f"{results['barcode_name']}_speciesmix_results.json")
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {json_path}")

    return results


def _persist_result_key(indiv_dir, bc, key, value):
    """Write a single key into a barcode's persisted speciesmix_results.json."""
    jp = os.path.join(indiv_dir, bc, f'{bc}_speciesmix_results.json')
    if os.path.exists(jp):
        try:
            with open(jp) as _jf:
                _jd = json.load(_jf)
            _jd[key] = value
            with open(jp, 'w') as _jf:
                json.dump(_jd, _jf, indent=2)
        except (OSError, ValueError):
            pass


def _find_bc_fastq(output_dir, bc):
    """Best-effort locate a barcode's read FASTQ for on-demand rescans. Runs
    differ in layout (fastq/<bc>/, <bc>/ at run root), so probe both; prefer a
    combined/trimmed file, else the first .fastq[.gz]."""
    cands = [os.path.join(output_dir, EXPERIMENT_CONFIG.get('fastq_dir', 'fastq'), bc),
             os.path.join(output_dir, bc)]
    for d in cands:
        if not os.path.isdir(d):
            continue
        files = os.listdir(d)
        for pref in (f'combined_{bc}.fastq', f'noadapter_combined_{bc}.fastq'):
            if pref in files:
                return os.path.join(d, pref)
        fq = [f for f in files if f.endswith(('.fastq', '.fastq.gz'))]
        if fq:
            fq.sort()
            return os.path.join(d, fq[0])
    return None


def _p5p7_in_reads(fastq_dir, bc_list, sample=3000, min_frac=0.30):
    """Evidence gate for the Illumina adapter QC: True when >= min_frac of the first
    barcode's trimmed reads carry a P5 or P7 seed. A keyword test on the free-text
    construct description would silently leave every QC row at n/a on runs whose
    description did not say P5/P7/Illumina."""
    try:
        try:
            import illumina_adapter_qc_illumina as _iqc_mod
        except ImportError:
            import illumina_adapter_qc as _iqc_mod
        p5, p7 = _iqc_mod.P5, _iqc_mod.P7
    except Exception:
        p5, p7 = 'AATGATACGGCGACCACCGAGAT', 'CAAGCAGAAGACGGCATACGAGAT'
    seeds = set()
    for m in (p5, p7):
        for s in (m, rc(m)):
            seeds.update(s[i:i + 12] for i in range(len(s) - 11))
    for bc in bc_list:
        fqd = os.path.join(fastq_dir, bc)
        fqs = [f for f in (os.listdir(fqd) if os.path.isdir(fqd) else []) if f.startswith('noadapter_') and f.endswith('.fastq')]
        if not fqs:
            continue
        n = hit = 0
        with open(os.path.join(fqd, fqs[0])) as fh:
            while n < sample:
                if not fh.readline():
                    break
                s = fh.readline().strip(); fh.readline(); fh.readline(); n += 1
                if any(s[i:i + 12] in seeds for i in range(0, max(len(s) - 11, 0), 4)):
                    hit += 1
        return n > 0 and hit / n >= min_frac
    return False



def _unit_sort_key(b):
    """Sort units numerically when they are barcodeNN, by name otherwise (per-sample units
    are sample names)."""
    t = b.replace('barcode', '')
    return (0, int(t), '') if t.isdigit() else (1, 0, b)


def _unit_short(b):
    """Display label: BC12 for a numeric barcode, the name itself otherwise."""
    t = b.replace('barcode', '')
    return ('BC' + t.zfill(2)) if t.isdigit() else b


def compute_rt_extent_all(all_results, output_dir):
    """DEFERRED RT extent (template-switch site to transcript 3' end) per barcode, from
    the STAR BAM (rt_extent.compute_rt_extent); persisted as results['rt_extent'] so PDF
    re-renders reuse it. Returns True if anything was computed."""
    bc_list = sorted(all_results.keys(), key=_unit_sort_key)
    indiv_dir = os.path.join(output_dir, 'individual-analyses')
    need = [bc for bc in bc_list if not all_results[bc].get('rt_extent')
            and os.path.exists(os.path.join(indiv_dir, bc, f'{bc}_star_combined.bam'))]
    if not need:
        return False
    ran = False
    try:
        import rt_extent_illumina as _rte
        _stx = _rte.st
        _gtf = _stx.parse_gtf(_stx.GTF); _rdna = _stx.load_rdna(_stx.RDNA_BED)
        _tx = _rte.load_tx_ends(_stx.GTF)
        for bc in need:
            r = _rte.compute_rt_extent(os.path.join(indiv_dir, bc, f'{bc}_star_combined.bam'), _gtf, _rdna, _tx)
            all_results[bc]['rt_extent'] = r
            _persist_result_key(indiv_dir, bc, 'rt_extent', r)
            ran = True
            print(f"  RT extent computed for {bc}: n={r['n']:,} median={r['median']} nt")
    except Exception as e:
        print(f"  RT extent skipped: {e}")
    return ran


def compute_illumina_qc_all(all_results, output_dir):
    """DEFERRED Illumina P5/P7 adapter QC (opposite-ends architecture, i5/i7 index
    homogeneity, P5->Read1 / P7->Read2 junction correctness). Computed per barcode
    from the noadapter FASTQ and persisted into each JSON. Returns True if anything
    was computed. Gated on the construct declaring an Illumina (P5/P7) architecture;
    a no-op (already-cached or non-Illumina run) returns False."""
    construct = (EXPERIMENT_CONFIG.get('construct_description') or '').upper()
    bc_list = sorted(all_results.keys(), key=_unit_sort_key)
    indiv_dir = os.path.join(output_dir, 'individual-analyses')
    fastq_dir = os.path.join(output_dir, 'fastq')
    if not any(k in construct for k in ('P5', 'P7', 'ILLUMINA')) and not _p5p7_in_reads(fastq_dir, bc_list):
        return False
    # Recompute if never run OR if the cache predates the current metric schema
    # (a new illumina_adapter_qc version added keys the cached dict lacks).
    _REQ_KEYS = ('idx_i5i7_opposite_ends_pct', 'truseq_r1r2_opposite_ends_pct')
    need = [bc for bc in bc_list
            if not all_results[bc].get('illumina_qc')
            or any(k not in all_results[bc]['illumina_qc'] for k in _REQ_KEYS)]
    if not need:
        return False
    ran = False
    try:
        import illumina_adapter_qc_illumina as _iqc
        for bc in need:
            fqd = os.path.join(fastq_dir, bc)
            fqs = [f for f in (os.listdir(fqd) if os.path.isdir(fqd) else [])
                   if f.startswith('noadapter_') and f.endswith('.fastq')]
            if not fqs:
                continue
            qc = _iqc.compute_illumina_qc(os.path.join(fqd, fqs[0]))
            all_results[bc]['illumina_qc'] = qc
            _persist_result_key(indiv_dir, bc, 'illumina_qc', qc)
            ran = True
            print(f"  illumina QC computed for {bc}: opp-ends={qc.get('p5p7_opposite_ends_pct')}%  "
                  f"i5={qc.get('i5_dominant')}({qc.get('i5_match_pct')}%)  "
                  f"P5->R1={qc.get('p5_read1_pct')}%  P7->R2={qc.get('p7_read2_pct')}%")
    except Exception as _e:
        print(f"  WARNING: illumina-QC computation skipped ({_e})")
    return ran


def compute_dup_funnel_all(all_results, output_dir):
    """DEFERRED confident PCR-duplicate funnel (same-gene -> +UMI -> +~start/~end within
    10bp -> +same strand) over protein-coding non-ribosomal mRNA reads with a callable
    6nt TSO UMI. The strand step isolates PCR-amplification duplicates from both-strand
    sequencing of one molecule. Persisted per barcode; gated on a UMI chemistry."""
    if 'UMI' not in (EXPERIMENT_CONFIG.get('construct_description') or '').upper():
        return False
    bc_list = sorted(all_results.keys(), key=_unit_sort_key)
    indiv_dir = os.path.join(output_dir, 'individual-analyses')
    need = [bc for bc in bc_list if 'dup_funnel' not in all_results[bc]]
    if not need:
        return False
    ran = False
    try:
        import star_taxonomy_illumina as _stax
        import umi_dedup_illumina as _ud
        _ud.configure_umi(EXPERIMENT_CONFIG.get('rt_primers_used'))
        gtf = _stax.parse_gtf(_stax.GTF)
        rdna = _stax.load_rdna(_stax.RDNA_BED)
        for bc in need:
          # Per-barcode isolation. A failure here must cost THIS barcode its
          # duplicate rows and say so, not silently truncate the loop for every barcode
          # after it. Still loud: the reason is printed per barcode, not swallowed.
          try:
            bam = os.path.join(indiv_dir, bc, f'{bc}_star_combined.bam')
            # If dedup_reads_illumina.py ran, the BAM here has ALREADY had its duplicates
            # removed — recomputing on it would report ~0% and read as a pristine library
            # rather than a filtered one. That module writes the real pre-filter metrics
            # to this sidecar before touching anything; prefer them.
            _side = os.path.join(indiv_dir, bc, f'{bc}_dup_prefilter.json')
            df = None
            if os.path.exists(_side):
                try:
                    with open(_side) as _sf:
                        df = json.load(_sf) or None
                    if df:
                        print(f"  dup metrics for {bc}: using PRE-deduplication sidecar "
                              f"({df.get('dup_pct')}% duplicates before filtering)")
                except (OSError, ValueError):
                    df = None
            if df is None:
                df = _ud.compute_dup_funnel(bam, gtf, rdna, W=10) if os.path.exists(bam) else None
            all_results[bc]['dup_funnel'] = df
            _ps = os.path.join(indiv_dir, bc, f'{bc}_prep_stats.json')
            if os.path.exists(_ps):
                try:
                    all_results[bc]['prep_stats'] = json.load(open(_ps))
                except Exception:
                    pass
            _persist_result_key(indiv_dir, bc, 'dup_funnel', df)
            ran = True
            # Do not index df['n_mrna'] unconditionally: the PRE-dedup SIDECAR has a
            # different shape and has no such key, so that would raise KeyError -- and
            # because the only except sits OUTSIDE this loop, the raise would escape and
            # every LATER barcode would silently lose its duplicate, complexity and
            # flowcell rows (invisible on a run with ONE barcode, whose rows are assigned
            # just above, before the print). Report whichever shape is actually present.
            if df:
                _fk = ('n_umi', 'n_mrna', 'dup_gene_pct', 'dup_umi_pct',
                       'dup_pos_pct', 'dup_strand_pct')
                if all(k in df for k in _fk):
                    print(f"  dup-funnel for {bc}: UMI {df['n_umi']:,}/{df['n_mrna']:,} "
                          f"| gene {df['dup_gene_pct']}% +UMI {df['dup_umi_pct']}% "
                          f"+pos {df['dup_pos_pct']}% +strand(PCR) {df['dup_strand_pct']}%")
                else:
                    print(f"  dup-funnel for {bc}: pre-dedup sidecar "
                          f"(dup_pct={df.get('dup_pct')}, "
                          f"unique={df.get('unique')}, total={df.get('total')})")
          except Exception as _bce:                      # noqa: BLE001
            print(f"  WARNING: dup metrics FAILED for {bc} ({_bce!r}); its duplicate, "
                  f"complexity and flowcell rows will be blank. Other barcodes continue.")
    except Exception as _e:
        print(f"  WARNING: dup-funnel computation skipped ({_e})")
    return ran


def _read_composition_chem():
    """Map this run to a read_end_composition chemistry key, or None to skip the
    read-composition section. Honors an explicit `read_composition_chem` guide
    field; otherwise infers 24plex from the species barcodes."""
    # ILLUMINA: the whole § Read Start & End Composition page is ONT adapter taxonomy —
    # its four panels are 'true FIRST/LAST element' and 'element AFTER/BEFORE the ONT
    # adapter'. There is no ONT adapter or native barcode in an Illumina read, and the
    # 24plex motif list it classifies against is the full TruSeq2 handle (the sequencing
    # primer, never in the read) — so it would paint the ~93% of reads that legitimately
    # start with a bobcode as "likely cDNA". Dropped entirely.
    if IS_ILLUMINA:
        return None
    c = EXPERIMENT_CONFIG.get('read_composition_chem')
    if c:
        return c
    tmap = EXPERIMENT_CONFIG.get('tso_species_map') or {}
    if _ba.is_24plex(tmap, EXPERIMENT_CONFIG.get('tso_arch')):
        return '24plex'
    return None


def compute_read_composition_all(all_results, output_dir):
    """DEFERRED read start/end composition (where each read begins/ends by
    adapter-primer element). Scans the UNTRIMMED FASTQ -- the ONT adapter must
    still be present -- so it uses combined_<bc>.fastq, NOT the noadapter file.
    Persisted per barcode; gated on a supported chemistry (read_end_composition)."""
    chem = _read_composition_chem()
    if chem is None:
        return False
    bc_list = sorted(all_results.keys(), key=_unit_sort_key)
    indiv_dir = os.path.join(output_dir, 'individual-analyses')
    fqdir = os.path.join(output_dir, EXPERIMENT_CONFIG.get('fastq_dir', 'fastq'))
    need = [bc for bc in bc_list if 'read_composition' not in all_results[bc]]
    if not need:
        return False
    ran = False
    try:
        import read_end_composition_illumina as _rec
        if chem not in _rec.CHEMISTRIES:
            print(f"  read-composition: chemistry '{chem}' not defined; section skipped")
            return False
        for bc in need:
            bcdir = os.path.join(fqdir, bc)
            fq = os.path.join(bcdir, f'combined_{bc}.fastq')
            if not os.path.exists(fq) and os.path.isdir(bcdir):
                cand = [f for f in os.listdir(bcdir) if f.endswith('.fastq')
                        and not f.startswith('noadapter_') and 'filtered' not in f]
                fq = os.path.join(bcdir, cand[0]) if cand else None
            comp = (_rec.to_json(_rec.analyse(fq, chem))
                    if fq and os.path.exists(fq) else None)
            all_results[bc]['read_composition'] = comp
            _persist_result_key(indiv_dir, bc, 'read_composition', comp)
            ran = True
            if comp:
                d = comp['diagnostics']
                print(f"  read-composition for {bc}: native-bc {d['native_bc'] or 'NONE'} "
                      f"| ONT start {d['pct_ont_start']}% end {d['pct_ont_end']}%"
                      + (f"  WARN {len(d['warnings'])}" if d['warnings'] else ""))
            else:
                print(f"  read-composition for {bc}: untrimmed FASTQ not found; skipped")
    except Exception as _e:
        print(f"  WARNING: read-composition computation skipped ({_e})")
    return ran


def generate_cross_barcode_summary(all_results, output_dir):
    """Generate landscape cross-barcode summary PDF for RT+TSO constructs."""
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import letter, landscape
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import inch
    from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, PageBreak, Image, KeepTogether
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # Sort barcodes numerically
    bc_list = sorted(all_results.keys(), key=_unit_sort_key)
    n_bc = len(bc_list)
    bc_short = {bc: _unit_short(bc) for bc in bc_list}

    # Folder tidy-up: all cross-barcode plots (+ the symm-ends cache) go into a
    # figures/ subdir instead of cluttering the run root. The PDF stays at the root
    # with a descriptive name: BarcodeSum_<YYMMDD>_BC<first>-<last>.pdf
    fig_dir = os.path.join(output_dir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)
    _run_base = os.path.basename(os.path.normpath(output_dir))
    _run_date = _run_base.split('_', 1)[0]
    if not _run_date.isdigit():
        _run_date = ''.join(ch for ch in _run_base[:6] if ch.isdigit()) or 'NA'
    # BarcodeSum_<date>_BC01-24.pdf for numeric units; named units (per-sample mode) use their names
    if all(b.replace('barcode', '').isdigit() for b in bc_list):
        _pdf_name = f"BarcodeSum_{_run_date}_BC{int(bc_list[0].replace('barcode', '')):02d}-{int(bc_list[-1].replace('barcode', '')):02d}.pdf"
    else:
        _pdf_name = f"BarcodeSum_{_run_date}_{bc_short[bc_list[0]]}-{bc_short[bc_list[-1]]}.pdf"
    # 24plex chemistry (26N-UMI RT primer: TruSeq1-26N-polydT). Mirrors the selector in
    # barcode_adapter.get_adapter — gates the reordered full-molecule structure figure.
    _is_24plex = _ba.is_24plex(EXPERIMENT_CONFIG.get('tso_species_map'),
                               EXPERIMENT_CONFIG.get('tso_arch'))

    # Barcode-accuracy filter funnel (5 cumulative steps) for the master summary.
    # Computed on-demand from each BC's STAR BAM and persisted into its JSON, so
    # it is calculated once per barcode and reused on later PDF regenerations.
    _indiv_dir = os.path.join(output_dir, 'individual-analyses')
    # Prep sidecar (raw-read accounting rows): the regenerated PDF rebuilds all_results
    # from the persisted JSON, which predates this key, so read it from disk here whenever
    # a barcode lacks it. Cheap, and independent of which pass is rendering.
    for _bc in list(all_results):
        if not (all_results[_bc] or {}).get('prep_stats'):
            _psp = os.path.join(_indiv_dir, _bc, f'{_bc}_prep_stats.json')
            if os.path.exists(_psp):
                try:
                    all_results[_bc]['prep_stats'] = json.load(open(_psp))
                except Exception:
                    pass
    _need_ff = [bc for bc in bc_list if not all_results[bc].get('filter_funnel')]
    if _need_ff:
        try:
            import star_taxonomy_illumina as _stax
            import filter_funnel_illumina as _ff
            _gtf = _stax.parse_gtf(_stax.GTF)
            _rdna = _stax.load_rdna(_stax.RDNA_BED)
            _term = _stax._get_terminal()
            for bc in _need_ff:
                bam = os.path.join(_indiv_dir, bc, f'{bc}_star_combined.bam')
                funnel = _ff.compute_filter_funnel(bam, _gtf, _rdna, _term)
                if funnel is None:
                    continue
                all_results[bc]['filter_funnel'] = funnel
                jp = os.path.join(_indiv_dir, bc, f'{bc}_speciesmix_results.json')
                if os.path.exists(jp):
                    try:
                        with open(jp) as _jf:
                            _jd = json.load(_jf)
                        _jd['filter_funnel'] = funnel
                        with open(jp, 'w') as _jf:
                            json.dump(_jd, _jf, indent=2)
                    except (OSError, ValueError):
                        pass
                print(f"  filter funnel computed for {bc}: "
                      f"final-step acc={funnel['steps'][-1].get('acc')}")
        except Exception as _e:
            print(f"  WARNING: filter-funnel computation skipped ({_e})")

    # Auto-detected bobcodes (architecture + named bobcodes from a read scan) —
    # on-demand for runs analyzed before this feature existed, so a plain PDF
    # regen backfills it. Persisted per BC; recomputed only when the key is absent.
    _need_db = [bc for bc in bc_list if 'detected_bobcodes' not in all_results[bc]]
    if _need_db:
        try:
            import scan_bobcodes_illumina as _scan
            _extra = _config_bobcode_extras()
            for bc in _need_db:
                fq = _find_bc_fastq(output_dir, bc)
                db = _scan.scan(fq, n=30000, extra=_extra) if fq else None
                all_results[bc]['detected_bobcodes'] = db
                _persist_result_key(_indiv_dir, bc, 'detected_bobcodes', db)
                if db and db.get('architecture'):
                    _top = ", ".join(f"{b['id'] or b['name']}({b['seq']})"
                                     for b in db['barcodes'][:3])
                    print(f"  bobcodes auto-detected for {bc}: {db['architecture']} | {_top}")
                elif fq is None:
                    print(f"  bobcodes auto-detect for {bc}: FASTQ not found; skipped")
        except Exception as _e:
            print(f"  WARNING: bobcode auto-detect skipped ({_e})")

    # Full ordered-molecule funnel (24plex §1 read-structure figure) — on-demand BAM
    # scan, persisted per BC like the filter funnel. insert>30 (STAR) -> not rRNA ->
    # not MT -> exact 1x 7mer -> polydT -> 26N UMI -> TruSeq1 -> TruSeq2. Replaces the
    # old TSO structure classification figure. Uses the GTF/rDNA caches for rRNA/MT.
    if _is_24plex:
        _need_fm = [bc for bc in bc_list if not all_results[bc].get('fullmol_funnel')]
        if _need_fm:
            try:
                import full_mol_funnel_illumina as _fmf
                import star_taxonomy_illumina as _stax
                _fm_gtf = _stax.parse_gtf(_stax.GTF)
                _fm_rdna = _stax.load_rdna(_stax.RDNA_BED)
                for bc in _need_fm:
                    bam = os.path.join(_indiv_dir, bc, f'{bc}_star_combined.bam')
                    fm = _fmf.compute_full_mol_funnel(bam, _fm_gtf, _fm_rdna)
                    if fm is None:
                        continue
                    all_results[bc]['fullmol_funnel'] = fm
                    _persist_result_key(_indiv_dir, bc, 'fullmol_funnel', fm)
                    print(f"  full-molecule funnel for {bc}: FULL={fm['steps'][-1]['n']:,} "
                          f"({fm['steps'][-1]['pct']}%)")
            except Exception as _e:
                print(f"  WARNING: full-molecule funnel computation skipped ({_e})")

    # Illumina P5/P7 adapter QC and the confident PCR-dup funnel are DEFERRED: they are
    # the two slow analyses, so they are NOT computed here. The driver renders the PDF
    # first (their rows read "not run yet"), then calls compute_illumina_qc_all() /
    # compute_dup_funnel_all() and regenerates this PDF with the values filled in.
    # (When already computed and cached in all_results, their rows render populated.)

    # Detect BOB-polydT chemistry — gates the BOB-polydT-specific layout tweaks
    # (RT primer rows hidden in overview, extra specificity rows on filter
    # pages, specificity-by-insert-length table, per-species top-mRNA tables).
    is_bob_polydt = any(
        all_results[bc].get('bob_polydt_classification') for bc in bc_list
    )
    # Condition labels from experiment config
    bc_condition = EXPERIMENT_CONFIG.get('condition_labels', {})
    bc_label = {bc: f"{bc_short[bc]} ({bc_condition[bc]})" if bc in bc_condition else bc_short[bc]
                for bc in bc_list}
    # Plot-axis / legend label: the sample nickname (condition label, e.g. "44C",
    # "0.25(dC)") when provided, else the BCxx short name. Used only in figures;
    # tables keep the BCxx identity.
    bc_nick = {bc: (bc_condition.get(bc) or bc_short[bc]) for bc in bc_list}

    page_w, page_h = landscape(letter)
    usable_w = page_w - 1.2 * inch  # margins

    styles = getSampleStyleSheet()
    title_style = ParagraphStyle('T', parent=styles['Title'], fontSize=14, spaceAfter=8)
    heading_style = ParagraphStyle('H', parent=styles['Heading2'], fontSize=11, spaceAfter=3, spaceBefore=6)
    body_style = ParagraphStyle('B', parent=styles['Normal'], fontSize=8, spaceAfter=3, leading=10)

    # Two-line BC column headers (BC name on top, condition wrapped below).
    # Using Paragraph + <br/> so reportlab actually wraps the text within the
    # cell width — plain strings get truncated/overlapping when conditions are
    # long (e.g. "BC09 (Total RNA 50C)") and there are many BCs.
    _hdr_font = 8 if n_bc <= 4 else (7 if n_bc <= 6 else 6)
    _header_style = ParagraphStyle(
        'BCHeader', parent=styles['Normal'],
        fontSize=_hdr_font, leading=_hdr_font + 1,
        alignment=1,  # center
        spaceAfter=0, spaceBefore=0,
        fontName='Helvetica-Bold')

    def _mk_bc_header(bc):
        # Escape ampersands etc. that may appear in condition labels.
        def esc(s):
            return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
        if bc in bc_condition:
            cond = esc(bc_condition[bc])
            return Paragraph(
                f"<b>{bc_short[bc]}</b><br/>"
                f"<font size='{max(_hdr_font - 1, 5)}'>({cond})</font>",
                _header_style)
        return Paragraph(f"<b>{bc_short[bc]}</b>", _header_style)

    bc_header = {bc: _mk_bc_header(bc) for bc in bc_list}

    def fmt(n): return f"{n:,}"
    def ipct(n, t): return f"{n/max(t,1)*100:.0f}%" if t > 0 else '-'

    def auto_col_widths(label_w, n_data, extra_label_w=0):
        data_w = (usable_w - label_w - extra_label_w) / max(n_data, 1)
        if extra_label_w > 0:
            return [extra_label_w, label_w] + [data_w] * n_data
        return [label_w] + [data_w] * n_data

    def _extract_number(val):
        """Extract a numeric value from a cell string for comparison."""
        if not isinstance(val, str) or val == '-':
            return None
        # Try to find a leading number (handles "123", "1,234", "45%", "123 (45%)", "22M 218H", etc.)
        s = val.replace(',', '')
        # Try plain number first
        m = re.match(r'^([\d.]+)', s)
        if m:
            return float(m.group(1))
        return None

    # Absolute color gradient for specificity scores (range [-1, +1]):
    #   value <= 0   -> magenta
    #   value == 0.5 -> light grey
    #   value >= 1   -> green
    # Cells in `specificity_gradient_rows` use this absolute scale (looking up
    # the numeric value from `cell_values`), instead of the default row-relative
    # ±30%-of-mean highlighting. This way the swap-rate specificity is colored
    # by its meaning (negative = anti-specific, +1 = perfect) rather than
    # by how it compares to other cells in the same row.
    def _interp(a, b, t):
        return int(a + (b - a) * t)

    def specificity_gradient_color(v):
        if v is None:
            return None
        if v <= 0:
            return colors.HexColor('#F08CC0')  # magenta
        if v >= 1.0:
            return colors.HexColor('#7BC97B')  # green
        if v <= 0.5:
            # interp magenta (F08CC0) -> light grey (E0E0E0)
            t = v / 0.5
            r = _interp(0xF0, 0xE0, t)
            g = _interp(0x8C, 0xE0, t)
            b = _interp(0xC0, 0xE0, t)
        else:
            # interp light grey (E0E0E0) -> green (7BC97B)
            t = (v - 0.5) / 0.5
            r = _interp(0xE0, 0x7B, t)
            g = _interp(0xE0, 0xC9, t)
            b = _interp(0xE0, 0x7B, t)
        return colors.HexColor(f'#{r:02x}{g:02x}{b:02x}')

    def _fit_cells(data, col_widths, font_size):
        """Wrap any plain-string cell wider than its column into a Paragraph so the row
        grows vertically instead of the text running into the neighbouring cells
        (the duplicate-rate row would otherwise do exactly that in a six-column report)."""
        if not col_widths:
            return data
        from reportlab.pdfbase.pdfmetrics import stringWidth
        from reportlab.lib.enums import TA_LEFT, TA_RIGHT
        from xml.sax.saxutils import escape as _esc
        out = []
        for ri, row in enumerate(data):
            new = []
            for ci, cell in enumerate(row):
                w = col_widths[ci] if ci < len(col_widths) else None
                if ri > 0 and isinstance(cell, str) and w and stringWidth(cell, 'Helvetica', font_size) > w - 6:
                    st = ParagraphStyle(f'fit{ri}_{ci}', fontName='Helvetica', fontSize=font_size,
                                        leading=font_size + 2, alignment=(TA_LEFT if ci == 0 else TA_RIGHT))
                    cell = Paragraph(_esc(cell), st)
                new.append(cell)
            out.append(new)
        return out

    def make_heatmap_table(data, col_widths=None, font_size=8,
                            specificity_gradient_rows=None, cell_values=None,
                            divider_rows=None, plain_rows=None, cell_colors=None,
                            grey_rows=None,
                            **kwargs):
        """Like make_table with default ±30%-of-row-mean highlighting, except
        rows listed in `specificity_gradient_rows` are colored on the absolute
        [-1, +1] specificity gradient using values pulled from `cell_values`.

        `grey_rows`: 1-based row indices to render DE-EMPHASISED (grey italic text,
        no shading) — the "not applicable on this platform" state. The row and its
        label string are KEPT so the layout is unchanged and the cross-run matrix,
        which joins on exact label text, still finds it.

        NOTE: this function takes **kwargs and ignores unknown keys, so passing e.g.
        extra_cmds here fails SILENTLY. grey_rows is therefore a real named parameter
        rather than something callers try to smuggle in.
        """
        data = _fit_cells(data, col_widths, font_size)
        t = Table(data, colWidths=col_widths)
        cmds = [
            ('FONTSIZE', (0, 0), (-1, -1), font_size),
            ('ALIGN', (1, 0), (-1, -1), 'RIGHT'), ('ALIGN', (0, 0), (0, -1), 'LEFT'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ('TOPPADDING', (0, 0), (-1, -1), 2), ('BOTTOMPADDING', (0, 0), (-1, -1), 2),
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#4472C4')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ]
        spec_rows = set(specificity_gradient_rows or [])
        spec_vals = cell_values or {}
        plain = set(plain_rows or [])
        cell_cols = cell_colors or {}
        grey = set(grey_rows or [])
        for _gr in grey:                       # de-emphasise: grey italic, no shading
            if 0 < _gr < len(data):
                cmds.append(('TEXTCOLOR', (0, _gr), (-1, _gr), colors.HexColor('#9A9A9A')))
                cmds.append(('FONTNAME', (0, _gr), (-1, _gr), 'Helvetica-Oblique'))
        for ri in range(1, len(data)):
            if ri in plain:
                # No gradient/default highlighting; apply only explicit per-cell
                # threshold colors (if any) supplied via cell_colors.
                for ci in range(1, len(data[ri])):
                    c = cell_cols.get((ri, ci))
                    if c is not None:
                        cmds.append(('BACKGROUND', (ci, ri), (ci, ri), colors.HexColor(c)))
                continue
            if ri in spec_rows:
                # Absolute gradient based on numeric value stored in cell_values
                for ci in range(1, len(data[ri])):
                    v = spec_vals.get((ri, ci))
                    if v is None:
                        continue
                    color = specificity_gradient_color(v)
                    if color is not None:
                        cmds.append(('BACKGROUND', (ci, ri), (ci, ri), color))
                continue
            # Default: row-relative ±30%-of-row-mean highlighting
            vals = []
            for ci in range(1, len(data[ri])):
                v = _extract_number(data[ri][ci])
                if v is not None:
                    vals.append((ci, v))
            if len(vals) >= 2:
                mean_v = sum(v for _, v in vals) / len(vals)
                if mean_v > 0:
                    for ci, v in vals:
                        if v > mean_v * 1.3:
                            cmds.append(('BACKGROUND', (ci, ri), (ci, ri), colors.HexColor('#C6EFCE')))
                        elif v < mean_v * 0.7:
                            cmds.append(('BACKGROUND', (ci, ri), (ci, ri), colors.HexColor('#F2CEEF')))
        for _dr in (divider_rows or ()):
            cmds.append(('LINEBELOW', (0, _dr), (-1, _dr), 1.5, colors.black))
        t.setStyle(TableStyle(cmds))
        return t

    # Accumulator for the per-filter specificity-summary page. Each entry maps
    # a filter-level label -> {bc: float}. Populated below as each section
    # computes its Overall (specificity) row, and rendered as a separate
    # summary page at the very end of the PDF.
    specificity_summary = []  # list of (label, {bc: spec_value or None})

    def make_table(data, col_widths=None, font_size=8, extra_cmds=None):
        data = _fit_cells(data, col_widths, font_size)
        t = Table(data, colWidths=col_widths)
        cmds = [
            ('FONTSIZE', (0, 0), (-1, -1), font_size),
            ('ALIGN', (1, 0), (-1, -1), 'RIGHT'), ('ALIGN', (0, 0), (0, -1), 'LEFT'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ('TOPPADDING', (0, 0), (-1, -1), 2), ('BOTTOMPADDING', (0, 0), (-1, -1), 2),
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#4472C4')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ]
        if extra_cmds:
            cmds.extend(extra_cmds)
        t.setStyle(TableStyle(cmds))
        return t

    pdf_path = os.path.join(output_dir, _pdf_name)
    doc = SimpleDocTemplate(pdf_path, pagesize=landscape(letter),
                            topMargin=0.5*inch, bottomMargin=0.5*inch,
                            leftMargin=0.6*inch, rightMargin=0.6*inch)
    elems = []

    exp_title = EXPERIMENT_CONFIG.get('experiment_title', 'Species Mix Analysis')
    construct_desc = EXPERIMENT_CONFIG.get('construct_description', '')
    loc = EXPERIMENT_CONFIG.get('barcode_location', 'BOBcode')
    primary_map = EXPERIMENT_CONFIG.get('primary_species_map', {})
    map_label = 'BOBcode barcodes' if loc == 'BOBcode' else 'TSO barcodes'
    map_desc = ', '.join(f'<b>{code}</b> ({sp})' for code, sp in primary_map.items())
    cond_desc = ' → '.join(f'{bc_short.get(bc, bc)} ({lbl})' for bc, lbl in bc_condition.items()
                            if bc in bc_list) if bc_condition else ''
    elems.append(Paragraph(exp_title, title_style))
    elems.append(Paragraph(
        "<i>How to read this report: the bold rows at the top of the Master Summary are the headline numbers; "
        "each later section gives the detail behind one of them. Where a table shows a metric at several filter "
        "levels, the row marked 'reported' is the one carried into the summary. Definitions of every number are in "
        "Methods &amp; Settings at the end.</i>", body_style))
    elems.append(Spacer(1, 6))
    # The report OPENS with the Master Summary table (inserted at the front
    # below). The construct/provenance line, the auto-detected-bobcodes table, and
    # the Master Summary methods text are collected here and appended as an
    # APPENDIX at the very back of the report.
    _construct_para = Paragraph(
        f"Construct: {construct_desc}&nbsp;&nbsp;|&nbsp;&nbsp;"
        f"{map_label}: {map_desc}"
        + (f"&nbsp;&nbsp;|&nbsp;&nbsp;{cond_desc}" if cond_desc else ''),
        body_style)

    # ===== Auto-detected bobcodes (config-free, from the reads) =====
    # Which chemistry (24plex / tso-bob / BOB-J) and which NAMED bobcodes each
    # sample used, read straight from the reads (scan_bobcodes) and matched to the
    # reference sheets. Collected here, rendered in the appendix at the back.
    _bobcode_top_elems = []
    _db_by_bc = {bc: all_results[bc].get('detected_bobcodes') for bc in bc_list}
    # Barcode-free chemistries (e.g. the NEB template-switch oligo) carry NO
    # sample barcode -- species is called by alignment. Some of those TSOs embed a
    # fixed motif (AGTACAT) that the tso-bob naming later reused as a "barcode", so
    # the scanner would otherwise mislabel that TSO constant as a sample bobcode.
    # When the config declares the run barcode-free, say so instead.
    _barcode_free = any(k in construct_desc.lower()
                        for k in ('barcode-free', 'barcode free', 'no tso barcode'))
    if any((_db_by_bc.get(bc) or {}).get('architecture') for bc in bc_list):
        _bobcode_top_elems.append(Paragraph(
            "Auto-detected bobcodes "
            "<i>&mdash; chemistry and named bobcodes identified from the reads "
            "(anchor motif &rarr; barcode matched to the reference set, &le;1 mismatch); "
            "% = share of scanned reads bearing that bobcode."
            + (" This library is <b>barcode-free</b> (species by alignment); the "
               "motif shown is the fixed TSO sequence, not a sample barcode."
               if _barcode_free else "") + "</i>", heading_style))
        _bl = EXPERIMENT_CONFIG.get('bobcode_labels') or {}
        # One row per bobcode, one '% reads' column per ONT barcode (no ONT sample or
        # chemistry columns; the chemistry is named in the heading).
        # The pooled samples sit behind every ONT barcode, so the codes are merged across them.
        _archs = sorted({(_db_by_bc.get(bc) or {}).get('architecture') for bc in bc_list if (_db_by_bc.get(bc) or {}).get('architecture')})
        _bobcode_top_elems.append(Paragraph(f"<i>Chemistry detected from the reads: <b>{', '.join(_archs) or 'not detected'}</b>. "
            "% reads = share of scanned reads carrying the bobcode, per ONT barcode.</i>", body_style))
        _codes = {}                     # seq -> {'name','id','species','frac': {bc: frac}}
        _cut_by_bc = {}
        for bc in bc_list:
            db = _db_by_bc.get(bc) or {}
            bars = db.get('barcodes') or []
            if not bars:
                continue
            _cut = float(EXPERIMENT_CONFIG.get('bobcode_show_frac', 0.10)) * bars[0]['frac']
            _cap = int(EXPERIMENT_CONFIG.get('bobcode_show_max', 4))
            _show = [b for b in bars if b['frac'] >= _cut][:_cap] or bars[:1]
            for b in _show:
                e = _codes.setdefault(b['seq'], {'name': b['name'], 'id': b['id'] or '', 'species': b['species'] or primary_map.get(b['seq']) or '—', 'frac': {}})
                e['frac'][bc] = b['frac']
        if _barcode_free and _codes:
            _codes = dict(list(_codes.items())[:1])
        _dbh = ['Sample' if _bl else 'Bobcode', 'ID', 'Sequence', 'Species'] + [f"% reads {bc_short.get(bc, bc)}" for bc in bc_list]
        _dbrows = [_dbh]
        for seq, e in sorted(_codes.items(), key=lambda kv: -max(kv[1]['frac'].values())):
            _who = 'barcode-free (no sample bobcode)' if _barcode_free else (_bl.get((seq or '').upper()) or e['name'])
            _dbrows.append([_who, e['id'], seq, 'by alignment' if _barcode_free else e['species']]
                           + [f"{e['frac'][bc]*100:.1f}%" if bc in e['frac'] else '—' for bc in bc_list])
        if len(_dbrows) == 1:
            _dbrows.append(['none matched', '', '', ''] + [''] * len(bc_list))
        _dbw = [1.9 * inch, 0.75 * inch, 1.3 * inch, 0.9 * inch] + [0.9 * inch] * len(bc_list)
        _bobcode_top_elems.append(make_table(
            _dbrows, col_widths=_dbw, font_size=7,
            extra_cmds=[('ALIGN', (0, 1), (3, -1), 'LEFT'),
                        ('FONTNAME', (2, 1), (2, -1), 'Courier')]))
        _bobcode_top_elems.append(Spacer(1, 10))

    # 1. Read Structure Classification
    # Dynamic section numbering: sections are numbered via _secn() so a
    # chemistry-conditional section (e.g. the TSO-only "TSO Structure" page) keeps
    # the numbering contiguous in both reports.
    _sec = [0]
    def _secn():
        _sec[0] += 1
        return _sec[0]

    # ===== Read Start & End Composition (deferred; directly after Master Summary) =====
    # Where each read begins/ends by adapter-primer element -- an ONT-library QC that
    # reads out ligation orientation and how completely reads span the molecule. Gated
    # on a supported chemistry; the data are computed by compute_read_composition_all
    # (deferred), so on the first render this shows "pending" and populates on regen.
    _rc_chem = _read_composition_chem()
    if _rc_chem is not None:
        elems.append(Paragraph(f"{_secn()}. Read Start &amp; End Composition", heading_style))
        elems.append(Paragraph(
            "<i>Where each read starts and ends, read directly from the untrimmed reads. Top row: the very first and "
                "very last element of the read. Bottom row: the first library element after the ONT adapter block and the "
                "last one before it, which tells which end of the molecule was read and in which direction. Elements are "
                "recognised by short exact seeds so that ONT errors do not hide them; the ONT native barcode and the i5/i7 "
                "indices are learned from the reads themselves (table below). 'likely cDNA' = no adapter found near that end.</i>", body_style))
        _rc_data = {bc: all_results[bc].get('read_composition') for bc in bc_list}
        if any(_rc_data.values()):
            try:
                from read_composition_section_illumina import render_grid
                _rcp = render_grid(all_results, bc_list, bc_short, fig_dir)
                if _rcp and os.path.exists(_rcp):
                    elems.append(Image(_rcp, width=usable_w, height=usable_w * 0.545))
            except Exception as _e:
                print(f"  (read-composition figure skipped: {_e})")
            elems.append(Spacer(1, 6))
            _c0 = next((r['diagnostics']['constants'] for r in _rc_data.values() if r), {})
            elems.append(Paragraph(f"<b>Learned from the reads</b> <i>— per barcode, not taken from the run config. "
                f"ONT adapter {_c0.get('ont_adapter','?')}; native-barcode flanks "
                f"{_c0.get('nb_flank_5','?')}/{_c0.get('nb_flank_3','?')}.</i>", body_style))
            # Indices in the sample-sheet convention: the 10-nt core upstream of
            # the Read1/Read2 primer from the adapter QC, i5 forward, i7 reverse-complemented
            # (= what the sequencing provider lists) with the in-primer form beside it. The
            # 12-nt window after the P5/P7 handle mixes linker and index, so it is not printed.
            _pv = [['BC', 'ONT native barcode', 'i5 index (forward)', 'i7 index (sample sheet / in primer)', 'ONT@start']]
            for bc in bc_list:
                _d = ((_rc_data[bc] or {}).get('diagnostics') or {})
                _q = all_results[bc].get('illumina_qc') or {}
                _i5 = (f"{_q['i5_dominant']} ({_q.get('i5_match_pct', 0):.0f}%)" if _q.get('i5_dominant') else '—')
                _i7 = (f"{_q['i7_dominant']} / {_q.get('i7_dominant_oligo', '?')} ({_q.get('i7_match_pct', 0):.0f}%)" if _q.get('i7_dominant') else '—')
                _pv.append([bc_short[bc], _d.get('native_bc') or '—', _i5, _i7, f"{_d.get('pct_ont_start','—')}%"])
            elems.append(make_table(_pv, col_widths=[0.6 * inch, 2.4 * inch, 1.8 * inch, 2.6 * inch, 0.9 * inch]))
            _warns = [f"{bc_short[bc]}: {w}" for bc in bc_list
                      for w in (((_rc_data[bc] or {}).get('diagnostics') or {}).get('warnings') or [])]
            if _warns:
                elems.append(Spacer(1, 4))
                elems.append(Paragraph("<b>Warnings:</b> " + " | ".join(_warns[:4]), body_style))
        else:
            elems.append(Paragraph("<i>Pending — computed on PDF regeneration (deferred: it "
                "scans the untrimmed FASTQ). Re-run to populate.</i>", body_style))
        elems.append(PageBreak())

    elems.append(Paragraph(f"{_secn()}. Read Structure Classification", heading_style))
    _is_tso = (EXPERIMENT_CONFIG.get('barcode_location') == 'TSO')
    # >=2-TSO chimeras can only be excluded where the TSO multiplicity scan ran (full-length
    # reads); on Illumina r2_positional it is disabled (star_taxonomy), so no label may
    # promise an exclusion that was not applied.
    _has_mult = _is_tso and any((all_results[bc].get('tso_multiplicity') or {}).get('n_total')
                                for bc in bc_list)
    # RT-primer-end label for this group ('C28', or 'R1' for R1_polydT_VN
    # groups). Substituted into the §1 caption + structure figure.
    _is_24 = _is_tso and _ba.is_24plex(EXPERIMENT_CONFIG.get('tso_species_map'), EXPERIMENT_CONFIG.get('tso_arch'))
    # 24plex has no C28 primer: its RT end is the R1 poly-dT.
    _rt_lbl = EXPERIMENT_CONFIG.get('rt_end_label') or ('R1 poly-dT' if _is_24 else 'C28')
    _tm24 = EXPERIMENT_CONFIG.get('tso_species_map') or {}
    _n_h24 = sum(1 for v in _tm24.values() if v == 'human'); _n_m24 = sum(1 for v in _tm24.values() if v == 'mouse')
    # Chemistry-aware descriptors reused in the Master-Summary caption + Methods.
    _full_struct_desc = ("a declared species code + STAR insert ≥50 bp + R1 poly-dT RT end (GGG not required" + ("; ≥2-TSO chimeras excluded" if _has_mult else "") + ")" if _is_24 else
                         f"single 7-mer barcode (≥5 bp Nextera-R1) + "
                         f"insert ≥50 bp + polyT/A ≥10 + {_rt_lbl} RT-end" + (" (≥2-TSO chimeras excluded)" if _has_mult else "") if _is_tso else
                         "code + polyT ≥5 + insert ≥20 bp + C28 ≥5 bp (TSO not required)")
    _codes_desc = ((", ".join(f"{k}→{v}" for k, v in sorted(_tm24.items())) if len(_tm24) <= 4 else
                    f"the run's species map: {_n_h24} human and {_n_m24} mouse codes, listed in the 'Barcodes used in this run' table") if _is_24 else
                   "ATCGAAA→human, CAGTTGA→mouse" if _is_tso else
                   "ATCGCT→human, TCGATA→mouse")
    _bc_loc_desc = ("a 7-mer (or longer) code in the TSO immediately 3′ of the TruSeq-R1 anchor CTCTTCCGATCT "
                    "(the end of the TruSeq2 handle), before the 6N UMI + spacer and the GGG template-switch" if _is_24 else
                    "a 7-mer in the TSO (after the Nextera-R1 backbone "
                    "GTGACTGGAGTTCAGACGTG, before the GGG template-switch)" if _is_tso else
                    "a 6-mer between the bob-c primer (…ATGTC) and the polyT")
    _fuzzy_desc = ("within ≤1 mismatch of a single declared code, longest code first, ties rejected (the 24plex "
                   "set has minimum pairwise Hamming distance 4, so ≤1 is unambiguous)" if _is_24 else
                   "fuzzy within ≤2 mismatches of a single 7-mer (the two codes are "
                   "Hamming-6 apart, so ≤2 of one means ≥4 from the other — unambiguous)" if _is_tso else
                   "fuzzy within ≤2 mismatches of a single code (the two codes are 6 "
                   "mismatches apart, so ≤2 is unambiguous)")
    # Barcode-correctness exclusion list (TSO additionally drops >=2-TSO chimeras,
    # which carry two conflicting barcodes). Reused in the §-notes and Methods.
    _excl_desc = "rRNA + MT + intergenic + low-complexity" + (
        " + ≥2-TSO chimeras" if _has_mult else "")

    # Redesigned for BOB-polydT: structural-completeness + "Retained" keep-set.
    # ONT clips the outer ends, so each terminal element is graded full/partial/
    # absent. Full structure = code + polyT>=5 + insert>=20bp + TSO>=8nt + C28>=5bp.
    _ret_chem = any(all_results[bc].get('bob_polydt_retention') for bc in bc_list)

    # ---- Unified chemistry-adaptive structure funnel -------------------------
    # One ordered step list per barcode: 24plex from fullmol_funnel (insert → not
    # rRNA → not MT → exact-1×-barcode → polydT → 26N UMI → TruSeq1 → TruSeq2);
    # tso-bob from tso7_retention (counts) + star_funnel_spec (accuracy) (has-7mer
    # → GGG → insert≥50 → C28 end). Absent elements are simply not in the source,
    # so they drop out. Drives BOTH the dedicated funnel tables below and the
    # Master-Summary per-step rows, so the same chemistry-appropriate steps show
    # everywhere.
    # Fixed CANONICAL funnel order — every run shows all rows; a step a chemistry
    # doesn't use renders blank. 24plex fills insert..ts2; tso-bob fills insert/
    # mrna/has7/exact/polyt/c28. Labels are display-side (override producer labels).
    _CANON_FUNNEL = [
        ('has7',   'has one barcode'),
        ('mrna',   '+ mRNA only'),
        ('dist',   '+ G-run present (distance enforced at prep)'),
        ('grun3',  '+ G-run ≥3'),
        ('insert', '+ insert >30 bp (STAR) = full structure'),
    ]
    def _uf(bc):
        """Canonical funnel for one barcode: {n_total, steps:[{key,label,present,
        n,acc,acc_n}]} in the fixed canonical order; a step the chemistry lacks has
        present=False (renders blank). 24plex from fullmol_funnel; tso-bob from
        star_funnel_spec (both keyed by canonical keys)."""
        r = all_results[bc]
        fm = r.get('fullmol_funnel')
        if fm and fm.get('steps'):
            src = {st.get('key'): st for st in fm['steps']}
            n_total = fm.get('n_total') or 0
            def _get(k):
                s = src.get(k)
                return None if s is None else {'n': s.get('n', 0), 'acc': s.get('acc'), 'acc_n': s.get('acc_n', 0)}
        else:
            sfs = r.get('star_funnel_spec') or {}
            n_total = (r.get('tso7_retention') or {}).get('n_total') or 0
            def _get(k):
                s = sfs.get(k)
                return None if s is None else {'n': s.get('n_reach', 0), 'acc': s.get('specificity'), 'acc_n': s.get('n', 0)}
        steps = []
        for k, lbl in _CANON_FUNNEL:
            d = _get(k)
            steps.append({'key': k, 'label': lbl, 'present': d is not None,
                          'n': (d or {}).get('n'), 'acc': (d or {}).get('acc'),
                          'acc_n': (d or {}).get('acc_n', 0)})
        return {'n_total': n_total, 'steps': steps}

    # Magenta→yellow→green ramp for single-step retention % (0% magenta, 50%
    # yellow, 100% green; continuous linear interpolation).
    def _fm_retain_color(p):
        try:
            _t = max(0.0, min(1.0, float(p) / 100.0))
        except (TypeError, ValueError):
            return colors.HexColor('#FFFFFF')
        _MAG = (0xD6, 0x19, 0x9C); _YEL = (0xFF, 0xE8, 0x4D); _GRN = (0x1A, 0x98, 0x50)
        if _t < 0.5:
            _s2 = _t / 0.5; _a, _b = _MAG, _YEL
        else:
            _s2 = (_t - 0.5) / 0.5; _a, _b = _YEL, _GRN
        _rgb = tuple(int(round(_a[_i] + (_b[_i] - _a[_i]) * _s2)) for _i in range(3))
        return colors.HexColor('#%02X%02X%02X' % _rgb)

    # Accuracy heatmap for the "barcode % correct by funnel level" table. Piecewise-
    # linear over anchor stops: magenta ≤90%, yellow @95%, light green @99%, deep
    # green ≥99.5%. Flat magenta below 90 and flat deep green above 99.5.
    def _fm_acc_color(p):
        try:
            _v = float(p)
        except (TypeError, ValueError):
            return colors.HexColor('#FFFFFF')
        _stops = [(90.0, (0xD6, 0x19, 0x9C)),   # magenta
                  (95.0, (0xFF, 0xE8, 0x4D)),   # yellow
                  (99.0, (0xB2, 0xE2, 0x9A)),   # light green
                  (99.5, (0x1A, 0x98, 0x50))]   # deep green
        if _v <= _stops[0][0]:
            _rgb = _stops[0][1]
        elif _v >= _stops[-1][0]:
            _rgb = _stops[-1][1]
        else:
            _rgb = _stops[0][1]
            for _i in range(len(_stops) - 1):
                _lo, _ca = _stops[_i]; _hi, _cb = _stops[_i + 1]
                if _lo <= _v <= _hi:
                    _s2 = (_v - _lo) / (_hi - _lo) if _hi > _lo else 0.0
                    _rgb = tuple(int(round(_ca[_j] + (_cb[_j] - _ca[_j]) * _s2)) for _j in range(3))
                    break
        return colors.HexColor('#%02X%02X%02X' % _rgb)

    def _emit_unified_funnel():
        """Canonical funnel — the SAME fixed step rows for every chemistry, blank where a
        step doesn't apply. ONE table: per barcode, reads at the level (single-step
        retention %, magenta→yellow→green) and barcode accuracy at the level (n scored,
        magenta ≤90 → yellow 95 → light 99 → deep green ≥99.5)."""
        _ufb = {bc: _uf(bc) for bc in bc_list}
        elems.append(Paragraph("<b>Full-structure funnel with barcode accuracy</b> <i>— per barcode, "
            "the reads passing each successive criterion (the % is single-step retention: reads "
            "surviving this filter ÷ reads at the previous applicable row, not cumulative) and, of the "
            "reads at that level that aligned to one species (" + _excl_desc + " excluded — see Methods), the "
            "fraction whose barcode species matches the aligned species (n = reads scored). Steps this "
            "construct does not have are not shown; the last row is the full-structure keep-set. "
            "Retention shades magenta (low) → yellow (~50%) → green (high); accuracy shades magenta ≤90% "
            "→ yellow 95% → light green 99% → deep green ≥99.5%.</i>", body_style))
        _nk = lambda bc: bc_nick.get(bc, bc_short.get(bc, bc))
        _fd = [['Step'] + sum([[f"{_nk(bc)} reads", f"{_nk(bc)} % correct"] for bc in bc_list], [])]
        _fd.append(['Total reads'] + sum([[f"{_ufb[bc]['n_total']:,}", ''] for bc in bc_list], []))
        _bg = []
        _prev = {bc: _ufb[bc]['n_total'] for bc in bc_list}   # last applicable count, per bc
        _ri = 1                                                # table row index of the last appended row
        _present = {_k for _si, (_k, _l) in enumerate(_CANON_FUNNEL)
                    if any(_ufb[bc]['steps'][_si]['present'] for bc in bc_list)}
        _rep_key = next((k for k in ('dist', 'mrna', 'has7') if k in _present), None)   # the Master Summary accuracy level
        for _si, (_k, _lbl) in enumerate(_CANON_FUNNEL):
            # Steps this construct does not have are not shown: an always-empty
            # row (e.g. the C28 end on 24plex) only raises questions.
            if not any(_ufb[bc]['steps'][_si]['present'] for bc in bc_list):
                continue
            _row = [_lbl + (' (reported)' if _k == _rep_key else '')]; _ri += 1
            for _bj, bc in enumerate(bc_list):
                s = _ufb[bc]['steps'][_si]
                if s['present'] and s['n'] is not None:
                    _n = s['n']; _p = _prev[bc]
                    _sp = (100.0 * _n / _p) if _p else 0.0
                    _row.append(f"{_n:,} ({_sp:.0f}%)")
                    _bg.append(('BACKGROUND', (1 + 2 * _bj, _ri), (1 + 2 * _bj, _ri), _fm_retain_color(_sp)))
                    _prev[bc] = _n
                else:
                    _row.append('-')
                _av = s.get('acc') if s['present'] else None
                if s['present'] and _av is not None:
                    _row.append(f"{_av:.2f}% (n={s.get('acc_n', 0):,})")
                    _bg.append(('BACKGROUND', (2 + 2 * _bj, _ri), (2 + 2 * _bj, _ri), _fm_acc_color(_av)))
                else:
                    _row.append('-')
            _fd.append(_row)
        elems.append(make_table(_fd, col_widths=auto_col_widths(150, 2 * n_bc), extra_cmds=_bg))
        elems.append(Spacer(1, 6))
        elems.append(Paragraph("<i>Cross-run and cross-chemistry comparisons use the chemistry-general funnel "
            "(one implementation for every chemistry) written to the run's funnel_comparable files, not this "
            "table.</i>", body_style))
        elems.append(Spacer(1, 10))

    # Chemistry-adaptive description of the Master-Summary funnel steps (first 4),
    # reused in the Methods + appendix captions so their prose matches the actual rows.
    _ms_step_desc = " → ".join(f"<b>{_s['label']}</b>"
                               for _s in (_uf(bc_list[0])['steps'][:4] if bc_list else [])) or "the funnel steps"

    if _is_tso:
        # TSO-bob §1 (pure-STAR): structure figure + funnel-count + funnel-correctness,
        # all from star_taxonomy's tso7_retention + star_funnel_spec (no per_read.json).
        _tr = {bc: (all_results[bc].get('tso7_retention') or {}) for bc in bc_list}
        _have_fm = _is_24plex and any(all_results[bc].get('fullmol_funnel') for bc in bc_list)
        if any(r.get('n_total') for r in _tr.values()) or _have_fm:
            if _have_fm:
                # 24plex: full ordered-molecule funnel (insert →…→ TruSeq2).
                # No explanatory preamble — the funnel row labels carry the step names,
                # and the ONT-specific prose (partial TruSeq matching because ONT clips
                # the outer ends) does not apply to these reads.
                elems.append(Spacer(1, 6))
            else:
                elems.append(Paragraph(
                    "<i>Each read is graded by how completely the canonical TSO molecule "
                    "survived: 7-mer barcode (≥5 bp Nextera-R1 backbone) "
                    f"+ cDNA insert ≥50 bp + polyT/A ≥10 + {_rt_lbl} RT-end (GGG template-switch "
                    "no longer required — see Methods). A read carrying ≥2 TSO "
                    "units (two 7-mer barcodes) is a chimera and is split out (brown), never "
                    "counted as Full structure. Pure-STAR: structure from the FASTQ, insert "
                    "length from the STAR alignment.</i>", body_style))
                elems.append(Spacer(1, 6))
                try:
                    from tso_structure_section_illumina import render_tso_structure_figure
                    _ps = render_tso_structure_figure(all_results, bc_list, bc_short, fig_dir,
                                                      rt_label=_rt_lbl)
                    if _ps and os.path.exists(_ps):
                        elems.append(Image(_ps, width=usable_w * 0.82, height=usable_w * 0.82 * 0.73))
                except Exception as _e:
                    print(f'  (TSO structure figure skipped: {_e})')
                elems.append(PageBreak())
            # UNIFIED chemistry-adaptive funnel: counts (single-step retention +
            # magenta→yellow→green gradient) + barcode-accuracy-by-level, both from
            # _uf(). Only the ordered step list differs by chemistry (24plex vs tso-bob).
            _emit_unified_funnel()
            elems.append(Paragraph("<i>Not applicable to this chemistry: the per-component completeness "
                "table and the TSO / C28 cumulative detected-length figure exist for the tso-bob and "
                "BOB-polydT constructs only (the 24plex construct has no C28 primer, and its TSO is the "
                "TruSeq2 handle); they are omitted here rather than rendered at zero.</i>", body_style))
            elems.append(Spacer(1, 6))
        else:
            elems.append(Paragraph("<i>TSO structural QC unavailable.</i>", body_style))
    elif not _ret_chem:
        elems.append(Paragraph(
            "<i>Structural QC metrics unavailable — run the structure step, "
            "then regenerate.</i>", body_style))
    else:
        elems.append(Paragraph(
            "<i>Each read is graded by how completely the canonical molecule "
            "survived. ONT clips the outer ends (5′ C28 primer, 3′ TSO), so each "
            "terminal element is scored <b>full / partial / absent</b> rather "
            "than present/absent. A read has the <b>full structure</b> (kept for "
            "downstream analysis) only if it has all of: a <b>code</b> (bobcode), "
            "<b>polyT ≥5</b>, <b>insert ≥20 bp</b>, and <b>C28 ≥5 bp</b>. The TSO is "
            "measured but <b>not required</b> (it is usually clipped at the 3′ "
            "end).</i>", body_style))
        elems.append(Spacer(1, 6))

        # Figure: construct schematic + per-BC stacked outcome + polyT histogram
        try:
            import sys as _sys
            _sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
            from polydt_structure_section_illumina import render_polydt_structure_figure
            _ps_png = render_polydt_structure_figure(
                all_results, bc_list, bc_short, fig_dir)
            if _ps_png and os.path.exists(_ps_png):
                elems.append(Image(_ps_png, width=usable_w * 0.82,
                                   height=usable_w * 0.82 * 0.73))
        except Exception as _e:
            print(f'  (polydt structure figure skipped: {_e})')

        # TSO + C28 cumulative detected-length (per BC), from the inside out:
        # TSO from the GGG, C28 from the barcode -> % of reads reaching >= each length
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as _plt
            _fig, _axs = _plt.subplots(1, 2, figsize=(11, 2.8))
            _cmap = _plt.get_cmap('tab10' if n_bc <= 10 else 'tab20')
            for _key, _ax, _ttl, _xmax, _full, _inner in (
                    ('tso', _axs[0], 'TSO', 25, 25, 'GGG'),
                    ('c28', _axs[1], 'C28', 20, 18, 'barcode')):
                for _i, bc in enumerate(bc_list):
                    h = ((all_results[bc].get('star_len_hist', {}) or {}).get(_key, {}) or {})
                    if not h:
                        continue
                    counts = {int(k): v for k, v in h.items()}
                    tot = sum(counts.values()) or 1
                    xs = list(range(0, _xmax + 1))
                    # cumulative from the inside out: % of reads reaching >= x nt
                    ys = [sum(c for L, c in counts.items() if L >= x) / tot * 100
                          for x in xs]
                    _ax.plot(xs, ys, '-', lw=1.1, color=_cmap(_i % 20), label=bc_nick[bc])
                _ax.set_xlim(0, _xmax); _ax.set_ylim(0, 100)
                _ax.set_xlabel(f'{_ttl} length from {_inner} outward (nt); full = {_full}')
                _ax.set_ylabel('% of reads reaching ≥ length')
                _ax.set_title(f'{_ttl} cumulative detected length per BC', fontsize=9,
                              fontweight='bold')
                for _s in ('top', 'right'):
                    _ax.spines[_s].set_visible(False)
            _axs[0].legend(fontsize=6, ncol=2, frameon=False)
            _plt.tight_layout()
            _tcp = os.path.join(fig_dir, 'star_tso_c28_length.png')
            _fig.savefig(_tcp, dpi=150, bbox_inches='tight'); _plt.close(_fig)
            elems.append(Paragraph("<b>TSO &amp; C28 cumulative detected length</b> "
                "<i>— ONT clips read ends, so the TSO (3′) and C28 (5′) adapters are "
                "often partial. Cumulative from the inside out (TSO from the GGG; C28 "
                "from the barcode): each line = % of code-bearing reads whose adapter "
                "extends to at least that length.</i>", body_style))
            elems.append(Image(_tcp, width=7.4 * inch, height=1.9 * inch))
            elems.append(Spacer(1, 8))
        except Exception as _e:
            print(f'  (TSO/C28 length plot skipped: {_e})')

        def _ret(bc):
            return all_results[bc].get('bob_polydt_retention', {}) or {}

        def _binpct(r, binkey, sub):
            b = r.get(binkey, {}) or {}
            denom = sum(b.values())
            return ipct(b.get(sub, 0), denom) if denom else '-'

        def _polyt_median(r):
            hist = r.get('polyt_hist', {}) or {}
            tot = sum(hist.values())
            if not tot:
                return '-'
            target, cum = tot / 2.0, 0
            for L in sorted(hist, key=lambda k: int(k)):
                cum += hist[L]
                if cum >= target:
                    return str(int(L))
            return '-'

        # Retention funnel
        elems.append(PageBreak())
        elems.append(Paragraph("<b>Full-structure funnel</b> "
            "<i>— reads passing each successive criterion (% of total reads); "
            "the final row is the full-structure keep-set.</i>", body_style))
        funnel_specs = [
            ('Total reads',             'n_total',           'raw'),
            ('Has code (bobcode)',      'code',              'funnel'),
            ('  + polyT ≥5',            'code_polyt',        'funnel'),
            ('  + insert ≥20 bp',    'code_polyt_insert', 'funnel'),
        ]
        fn_data = [['Metric'] + [bc_header[bc] for bc in bc_list]]
        for label, key, src in funnel_specs:
            row = [label]
            for bc in bc_list:
                r = _ret(bc); tot = r.get('n_total', 0)
                if src == 'raw':
                    row.append(f"{r.get(key, 0):,}")
                else:
                    v = (r.get('funnel', {}) or {}).get(key, 0)
                    row.append(f"{v:,} ({ipct(v, tot)})")
            fn_data.append(row)
        elems.append(make_table(fn_data, col_widths=auto_col_widths(180, n_bc)))
        elems.append(Spacer(1, 10))

        # Parallel table: barcode specificity (% correct) at each funnel level —
        # does tightening the retention filter improve barcode specificity?
        elems.append(Paragraph("<b>Barcode % correct by funnel level</b> "
            "<i>— % of species-typed reads (STAR; " + _excl_desc + " "
            "excluded — see Methods) correctly barcoded by species at each funnel level "
            "(n = reads scored). STAR-based, so it matches the Master Summary.</i>",
            body_style))
        lvl_specs = [
            ('Has code (bobcode)',       'code'),
            ('  + polyT ≥5',             'code_polyt'),
            ('  + insert ≥20 bp',     'code_polyt_insert'),
        ]
        bs_data = [['Read subset'] + [bc_header[bc] for bc in bc_list]]
        for label, lv in lvl_specs:
            row = [label]
            for bc in bc_list:
                sl = ((all_results[bc].get('star_funnel_spec', {}) or {}).get(lv, {}) or {})
                n = sl.get('n', 0); spec = sl.get('specificity')
                row.append(f'{spec:.2f}% (n={n:,})' if (n and spec is not None) else '-')
            bs_data.append(row)
        elems.append(make_table(bs_data, col_widths=auto_col_widths(180, n_bc)))
        elems.append(Spacer(1, 10))

        # Per-component completeness
        elems.append(Paragraph("<b>Per-component completeness</b> "
            "<i>— among code-bearing reads. C28: full ≥18/20 · partial 5–17 · "
            "absent &lt;5.  TSO: full ≥22/25 · partial 8–21 · absent &lt;8.</i>",
            body_style))
        comp_specs = [
            ('bobcode detected',            lambda r: ipct((r.get('funnel', {}) or {}).get('code', 0), r.get('n_total', 0))),
            ('polyT present (≥5)',          lambda r: ipct(r.get('polyt_present', 0), (r.get('funnel', {}) or {}).get('code', 0))),
            ('polyT median (nt)',           _polyt_median),
            ('TSO — full (≥22)',            lambda r: _binpct(r, 'tso_bins', 'full')),
            ('TSO — partial (8–21)',        lambda r: _binpct(r, 'tso_bins', 'partial')),
            ('TSO — absent (<8)',           lambda r: _binpct(r, 'tso_bins', 'absent')),
        ]
        comp_data = [['Component'] + [bc_header[bc] for bc in bc_list]]
        for label, fn in comp_specs:
            comp_data.append([label] + [fn(_ret(bc)) for bc in bc_list])
        elems.append(make_table(comp_data, col_widths=auto_col_widths(180, n_bc)))
        elems.append(Spacer(1, 8))

    # Cell shading: white (0%) -> yellow (10%) -> green (>=20%, clamped).
    def _wyg_0_20(pct):
        t = max(0.0, min(1.0, pct / 20.0))
        if t <= 0.5:
            u = t / 0.5
            r, g, b = 255.0, 255.0 - 14.0 * u, 255.0 - 137.0 * u
        else:
            u = (t - 0.5) / 0.5
            r, g, b = 255.0 - 132.0 * u, 241.0 - 40.0 * u, 118.0 + 5.0 * u
        return colors.Color(r / 255.0, g / 255.0, b / 255.0)

    def _wyg_0_100(pct):
        # Same white(0)->yellow(mid)->green(max) ramp as _wyg_0_20 but on a 0-100 scale.
        t = max(0.0, min(1.0, pct / 100.0))
        if t <= 0.5:
            u = t / 0.5
            r, g, b = 255.0, 255.0 - 14.0 * u, 255.0 - 137.0 * u
        else:
            u = (t - 0.5) / 0.5
            r, g, b = 255.0 - 132.0 * u, 241.0 - 40.0 * u, 118.0 + 5.0 * u
        return colors.Color(r / 255.0, g / 255.0, b / 255.0)

    _tcat_style = ParagraphStyle('TaxCat', parent=styles['Normal'],
                                 fontSize=7, leading=8, alignment=0)

    def _wrap_first_col(data):
        for _ri in range(1, len(data)):
            _lab = data[_ri][0]
            if isinstance(_lab, str):
                data[_ri][0] = Paragraph(
                    _lab.replace('&', '&amp;').replace('<', '&lt;')
                        .replace('>', '&gt;'), _tcat_style)


    # ============ TSO Structure & Multiplicity (TSO-bob only) ===============
    # ILLUMINA: removed entirely. Multiplicity/structure needs to see TSO units INSIDE
    # one read; a 150 bp mate only ever sees one end of the molecule, so every number on
    # this page is either structurally zero or an artifact of chance 7-mer matches.
    # What a 'TSO unit' is depends on the chemistry: 24plex = find_24plex_units (exact TruSeq-R1
    # anchor + declared code within 1 mismatch); classic tso-bob = the Nextera backbone scan.
    _tso_unit_desc = (
        "the TruSeq-R1 anchor CTCTTCCGATCT (exact) followed by a declared species code "
        "(≤1 mismatch, longest code first, ties rejected; hits within 27 nt on one strand collapse to one unit)"
        if _ba.is_24plex(EXPERIMENT_CONFIG.get('tso_species_map'), EXPERIMENT_CONFIG.get('tso_arch')) else
        "a Nextera-R1 backbone (≥5 bp) carrying a canonical 7-mer barcode (≤2 mismatches)")
    if (not IS_ILLUMINA) and _is_tso and any((all_results[bc].get('tso_multiplicity') or {}).get('n_total')
                       for bc in bc_list):
        _MM = {bc: (all_results[bc].get('tso_multiplicity') or {}) for bc in bc_list}
        elems.append(PageBreak())
        elems.append(Paragraph(f"{_secn()}. TSO Structure &amp; Multiplicity", heading_style))
        elems.append(Paragraph(
            "<i>Each read is scanned for <b>TSO units</b> — " + _tso_unit_desc + ". The canonical "
            "molecule has exactly <b>one</b> TSO (5′ end); reads with <b>≥2</b> TSO are "
            "TSO–TSO chimeras. For 2-TSO reads we ask whether the two barcodes are the <b>same</b> "
            "barcode twice or <b>two different</b> barcodes (any two samples). The black tick in the "
            "right panel is the % different expected if molecules paired at random (1 − Σ p<sub>i</sub>"
            "<super>2</super> over the barcode frequencies): at the tick, chimeras join random molecules; "
            "below it, molecules of one sample join preferentially.</i>", body_style))
        elems.append(Spacer(1, 6))
        try:
            from tso_multiplicity_section_illumina import render_tso_multiplicity_figure
            _mp = render_tso_multiplicity_figure(all_results, bc_list, bc_short, fig_dir)
            if _mp and os.path.exists(_mp):
                from reportlab.lib.utils import ImageReader as _IR
                _iw, _ih = _IR(_mp).getSize()
                elems.append(Image(_mp, width=usable_w, height=usable_w * _ih / _iw))
        except Exception as _e:
            print(f'  (TSO multiplicity figure skipped: {_e})')
        elems.append(Spacer(1, 8))

        def _cd(bc, k):
            m = _MM[bc]; tot = max(m.get('n_total', 0), 1)
            return 100.0 * (m.get('count_dist', {}).get(k, 0)) / tot

        def _pp(v):
            return '-' if v is None else f'{v:.0f}%'
        _td = [['Metric'] + [bc_header[bc] for bc in bc_list]]
        _td.append(['Total reads'] + [f"{_MM[bc].get('n_total', 0):,}" for bc in bc_list])
        _td.append(['1 TSO — canonical'] + [f"{_cd(bc, '1'):.0f}%" for bc in bc_list])
        _td.append(['2 TSO — chimera'] + [f"{_cd(bc, '2'):.0f}%" for bc in bc_list])
        _td.append(['≥3 TSO'] + [f"{_cd(bc, '3+'):.0f}%" for bc in bc_list])
        _td.append(['  of 2-TSO: two different barcodes'] + [_pp(_MM[bc].get('diff_code2_pct')) for bc in bc_list])
        _td.append(['  of 2-TSO: the same barcode twice'] + [_pp(_MM[bc].get('same_code2_pct')) for bc in bc_list])
        _td.append(['  expected different if molecules pair at random'] + [_pp(_MM[bc].get('exp_diff_code2_pct')) for bc in bc_list])
        if any(_MM[bc].get('species_split_informative') for bc in bc_list):   # both species >= 10% of codes
            _td.append(['  of 2-TSO: different species'] + [_pp(_MM[bc].get('diff2_pct')) for bc in bc_list])
            _td.append(['  expected different species (random)'] + [_pp(_MM[bc].get('exp_diff2_pct')) for bc in bc_list])
        elems.append(make_table(_td, col_widths=auto_col_widths(200, n_bc)))
        elems.append(Spacer(1, 6))
        _m0 = _MM[bc_list[0]]; _nk = bc_nick.get(bc_list[0], bc_short.get(bc_list[0], bc_list[0]))
        _r2 = (_m0.get('prod_ge2') or {}).get('retained_pct'); _r1 = (_m0.get('prod_one') or {}).get('retained_pct')
        _d2 = _m0.get('diff_code2_pct'); _e2 = _m0.get('exp_diff_code2_pct')
        # Every claim below is read from this run's table; nothing is hard-coded from
        # a single run.
        _pair_note = ('' if _d2 is None or _e2 is None else
                      'far above the random expectation, i.e. molecules pair preferentially across samples' if _d2 >= _e2 + 10 else
                      'at the random expectation, i.e. any two molecules can join' if abs(_d2 - _e2) < 10 else
                      'below the random expectation, i.e. molecules of the same sample join preferentially')
        elems.append(Paragraph(
            "<i>≥2-TSO reads carry two barcodes and are head-to-head TSO–TSO chimeras (one 5′ TSO + one "
            "reverse-complement 3′ TSO). On " + _nk + f": {_pp(_d2)} of 2-TSO reads carry <b>two different barcodes</b> "
            f"against {_pp(_e2)} expected if molecules paired at random ({_pair_note})"
            + ". Because their two barcodes conflict, ≥2-TSO reads are <b>excluded from all barcode-correctness "
            "metrics</b> (the Master Summary accuracy row and the accuracy tables of the funnel and RNA-class sections). TSO unit = " + _tso_unit_desc + ".</i>", body_style))

    # ---- Pure-STAR RNA-class composition & rRNA-excluded swap --------------
    # Independent of minimap2 (results['star_taxonomy'] / ['star_swap'] from
    # star_taxonomy_illumina.py): species = STAR combined-genome chromosome prefix
    # (HUMAN_/MOUSE_), reads multimapping across both species = ambiguous;
    # categories from the combined GTF; rRNA flagged via data-driven rdna_loci.bed
    # and EXCLUDED from the swap metric (rDNA is asymmetric: GRCh38 has it, GRCm39
    # mostly does not, so STAR otherwise force-types rRNA reads to human).
    _has_star = any(all_results[bc].get('star_taxonomy') for bc in bc_list)
    _star_order = next((all_results[bc].get('star_taxonomy_order')
                        for bc in bc_list
                        if all_results[bc].get('star_taxonomy_order')), [])
    if _has_star:
        elems.append(PageBreak())
        elems.append(Paragraph(
            f"{_secn()}. RNA-class Composition &amp; Species-barcode Swap (STAR)",
            heading_style))
        elems.append(Paragraph(
            f"<i>Reads aligned with {METHODS_ALIGNER_BIN} to the combined human+mouse genome. Species = the genome a read's alignments are on; a "
            "read that aligns to both genomes is <b>ambiguous</b>. Each mapped read gets one RNA class, decided in this order: "
            "any alignment on a known rDNA locus → rRNA; mitochondrial contig → Mitochondrial; no gene → Intergenic; a gene but "
            "no exon → pre-mRNA / intronic; an rRNA gene → rRNA; a protein-coding gene → mRNA, or Ribosomal-protein mRNA when "
            "the gene is RP*/MRP*; anything else → Other ncRNA. The mRNA class therefore excludes ribosomal-protein genes, "
            "while the barcode-accuracy rows count every protein-coding gene. Unmapped reads are not classified.</i>",
            body_style))

        def _star_cat_total(bc, cat):
            c = (all_results[bc].get('star_taxonomy', {}) or {}).get(cat, {}) or {}
            return sum(c.values())

        def _star_mapped(bc):
            t = all_results[bc].get('star_taxonomy', {}) or {}
            return sum(sum(v.values()) for v in t.values())

        # ---- STAR Table 1: molecular composition per barcode (mapped reads),
        #      with STAR Table 3 (barcode specificity by category) right below it
        #      so the two share a page. -----------------------------------------
        _comp_cap = Paragraph(
            "<b>RNA classes per barcode (mapped reads)</b> <i>— every mapped read sits in exactly one class, so a column "
            "adds up to ~100%. Cell = "
            "reads (% of mapped). Shading white 0% → yellow 10% → green ≥20%. "
            "Unmapped reads are not categorised and have no row; mapped reads in this table: "
            + '; '.join(f"{bc_nick.get(bc, bc_short.get(bc, bc))} {_star_mapped(bc):,} of {(all_results[bc].get('dup_funnel') or {}).get('dedup_kept') or all_results[bc].get('total_reads', 0) or 0:,} reads after deduplication "
                        f"({100 * _star_mapped(bc) / max((all_results[bc].get('dup_funnel') or {}).get('dedup_kept') or all_results[bc].get('total_reads', 0) or 0, 1):.1f}%)" for bc in bc_list) + ".</i>",
            body_style)
        st_data = [['Molecular category (STAR)'] + [bc_header[bc] for bc in bc_list]]
        st_pcts = {}
        for cat in [c for c in _star_order if c != 'Unmapped / contaminant']:
            row = [cat]
            for ci, bc in enumerate(bc_list):
                n = _star_cat_total(bc, cat)
                m = _star_mapped(bc)
                st_pcts[(len(st_data), ci + 1)] = n / max(m, 1) * 100
                row.append(f"{n:,} ({ipct(n, m)})")
            st_data.append(row)
        _wrap_first_col(st_data)
        st_table = make_table(st_data, col_widths=auto_col_widths(160, n_bc))
        st_table.setStyle(TableStyle(
            [('BACKGROUND', (ci, ri), (ci, ri), _wyg_0_20(pct))
             for (ri, ci), pct in st_pcts.items()]))

        # ---- mRNA species breakdown (human / mouse / ambiguous) per barcode ----
        # Splits the 'mRNA (protein-coding)' row of the composition table above by
        # STAR-assigned species (columns sum to ~100% of mRNA reads) — the clean
        # species-mixing ratio, restricted to protein-coding mRNA.
        # The 'mRNA species breakdown' table is not rendered: its numbers are the Master
        # Summary '% mRNA human / mouse' rows and the mRNA row of the assignability table.
        _sp_cap = Paragraph(
            "<b>Barcode accuracy per RNA class</b> <i>— of the barcoded reads in each class that aligned to one species, "
            "the share whose barcode names that species. The <b>Overall</b> row is the headline accuracy: it leaves out rRNA "
            "(cannot be assigned to a species on the combined genome), mitochondrial reads, intergenic reads (poly-dT priming "
            "on genomic A/T tracts maps to either species), low-complexity alignments"
            + (" and reads with two or more TSO units (two conflicting barcodes)" if _has_mult else "")
            + ". The class rows above it are unfiltered, so the effect of each exclusion is visible. "
            "Magenta 0 → grey 50 → green 100.</i>",
            body_style)
        _spec_cats = [c for c in _star_order
                      if c not in ('rRNA', 'Unmapped / contaminant')]
        sp_data = [['Molecular category (STAR)'] + [bc_header[bc] for bc in bc_list]]
        sp_cells = {}; sp_grad = set()
        for _ri, cat in enumerate(_spec_cats, start=1):
            row = [cat]
            for _ci, bc in enumerate(bc_list):
                pc = ((all_results[bc].get('star_swap', {}) or {}).get(
                    'per_category', {}) or {}).get(cat, {})
                s = pc.get('specificity')
                row.append(f"{s:.2f}%" if isinstance(s, (int, float)) else '-')
                if isinstance(s, (int, float)):
                    sp_cells[(_ri, _ci + 1)] = s / 100.0
            sp_data.append(row); sp_grad.add(_ri)
        _ov_ri = len(sp_data)
        ov_row = ['Overall (rRNA/MT/intergenic/low-cx' + ('/2-TSO' if _has_mult else '') + ' excl.)']
        for _ci, bc in enumerate(bc_list):
            ov = ((all_results[bc].get('star_swap', {}) or {}).get(
                'overall_rrna_excluded', {}) or {})
            s = ov.get('specificity')
            ov_row.append(f"{s:.2f}%" if isinstance(s, (int, float)) else '-')
            if isinstance(s, (int, float)):
                sp_cells[(_ov_ri, _ci + 1)] = s / 100.0
        sp_data.append(ov_row); sp_grad.add(_ov_ri)
        _wrap_first_col(sp_data)
        sp_data[_ov_ri][0] = Paragraph('<b>Overall (rRNA/MT/intergenic/low-cx'
                                       + ('/2-TSO' if _has_mult else '') + ' excl.) (reported)</b>', _tcat_style)
        sp_table = make_heatmap_table(
            sp_data, col_widths=auto_col_widths(165, n_bc), font_size=7,
            specificity_gradient_rows=sp_grad, cell_values=sp_cells,
            divider_rows=(len(_spec_cats),))   # line above the Overall row
        elems.append(KeepTogether([_comp_cap, st_table, Spacer(1, 6),
                                   _sp_cap, sp_table]))
        elems.append(Spacer(1, 8))

        # ---- RNA-class composition stacked bar (per BC) ---------------------
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as _plt
            _CATCOL = {
                'mRNA (protein-coding)': '#2CA02C',
                'Ribosomal-protein mRNA': '#8C9F4F',
                'Other ncRNA (lncRNA/NMD/ret-intron/pseudo)': '#A4D69E',
                'rRNA': '#FF7F0E', 'Mitochondrial': '#9467BD',
                'pre-mRNA / intronic': '#1F77B4', 'Intergenic / genomic': '#8C564B',
                'Unmapped / contaminant': '#1A1A1A',
            }
            _fig, _ax = _plt.subplots(figsize=(8, 3.0))
            _xs = range(len(bc_list)); _bottom = [0] * len(bc_list)
            for cat in reversed(_star_order):
                if cat == 'Unmapped / contaminant':
                    continue                       # drop unmapped from the histogram
                vals = []
                for bc in bc_list:
                    t = all_results[bc].get('star_taxonomy', {}) or {}
                    m = sum(sum(v.values()) for v in t.values()) or 1
                    vals.append(sum((t.get(cat, {}) or {}).values()) / m * 100)
                _ax.bar(_xs, vals, bottom=_bottom, color=_CATCOL.get(cat, '#cccccc'),
                        label=cat, width=0.7)
                _bottom = [b + v for b, v in zip(_bottom, vals)]
            _ax.set_xticks(list(_xs))
            _ax.set_xticklabels([bc_nick[bc] for bc in bc_list], fontsize=7)
            _ax.set_ylabel('% of mapped reads'); _ax.set_ylim(0, 100)
            _h, _l = _ax.get_legend_handles_labels()
            _ax.legend(_h[::-1], _l[::-1], fontsize=6, loc='center left',
                       bbox_to_anchor=(1.0, 0.5), frameon=False)
            _ax.set_title('RNA-class composition by barcode (STAR, mapped reads)',
                          fontsize=9, fontweight='bold')
            for _s in ('top', 'right'):
                _ax.spines[_s].set_visible(False)
            _plt.tight_layout()
            _ccp = os.path.join(fig_dir, 'star_composition_bar.png')
            _fig.savefig(_ccp, dpi=150, bbox_inches='tight'); _plt.close(_fig)
            elems.append(Image(_ccp, width=7 * inch, height=2.6 * inch))
            elems.append(Spacer(1, 8))
        except Exception as _e:
            print(f'  (composition bar skipped: {_e})')

        # ---- STAR Table 2: species assignability (pooled across BCs) ---------
        elems.append(PageBreak())
        elems.append(Paragraph(
            "<b>Which genome the reads aligned to, by RNA class (all barcodes pooled)</b> <i>— "
            "<b>human</b> / <b>mouse</b>: every alignment of the read is on that genome; "
            "<b>ambiguous</b>: the read aligns to both genomes, so its species cannot be called. "
            "Cells give reads and their share of the class (row); the total column gives the class "
            "size as a share of all mapped reads; the bottom row is the species split of the whole "
            "library. Two things to read off: the ambiguous column is almost entirely rRNA (rRNA is "
            "near-identical between the two species and the mouse assembly lacks most of its rDNA), "
            "which is why rRNA is excluded from every barcode-accuracy number; and the mRNA row is "
            "the cleanest estimate of the species mix of real transcripts.</i>", body_style))
        _SPK2 = ('human', 'mouse', 'ambiguous')      # 'n.a.' column dropped: always 0 (mapped reads only)
        _spool = {cat: {k: 0 for k in _SPK2} for cat in _star_order}
        for bc in bc_list:
            tx = all_results[bc].get('star_taxonomy', {}) or {}
            for cat in [c for c in _star_order if c != 'Unmapped / contaminant']:
                for k, v in (tx.get(cat, {}) or {}).items():
                    if k in _spool[cat]:
                        _spool[cat][k] += v
        _sgrand = sum(sum(c.values()) for c in _spool.values()) or 1
        ss_data = [['RNA class', 'human (% of class)', 'mouse (% of class)', 'ambiguous (% of class)', 'total (% of mapped)']]
        ss_shade = {}
        _scoltot = {k: 0 for k in _SPK2}
        for cat in [c for c in _star_order if c != 'Unmapped / contaminant']:
            c = _spool[cat]
            rowtot = sum(c.values())
            ri = len(ss_data)
            cells = [cat]
            for ci, k in enumerate(_SPK2, start=1):
                n = c[k]
                _scoltot[k] += n
                pct = n / _sgrand * 100
                cells.append(f"{n:,} ({100 * n / rowtot:.0f}%)" if n else '–')   # share of the class
                # no abundance shading: it would encode read counts, not assignability
            cells.append(f"{rowtot:,} ({rowtot / _sgrand * 100:.1f}% of mapped)")
            ss_data.append(cells)
        _sct_row = ['Column total']
        for k in _SPK2:
            _sct_row.append(f"{_scoltot[k]:,} ({_scoltot[k] / _sgrand * 100:.1f}%)")
        _sct_row.append(f"{_sgrand:,} (100%)")
        ss_data.append(_sct_row)
        _sct_ri = len(ss_data) - 1
        _wrap_first_col(ss_data)
        ss_data[_sct_ri][0] = Paragraph('<b>Column total</b>', _tcat_style)
        ss_table = make_table(ss_data, col_widths=auto_col_widths(165, 4))
        _ss_style = [('BACKGROUND', (ci, ri), (ci, ri), _wyg_0_20(pct))
                     for (ri, ci), pct in ss_shade.items()]
        _ss_style += [
            ('LINEABOVE', (0, _sct_ri), (-1, _sct_ri), 1.0, colors.black),
            ('FONTNAME', (1, _sct_ri), (-1, _sct_ri), 'Helvetica-Bold'),
        ]
        ss_table.setStyle(TableStyle(_ss_style))
        elems.append(ss_table)
        elems.append(Spacer(1, 8))

        # (Barcode specificity by molecular category table moved up to sit
        #  directly under the composition table — see STAR Table 1 above.)

        # (No 3C barcode × species crosstab — its rRNA-included raw split would
        #  duplicate the rRNA-excluded per-category specificity table above.
        #  No barcode % correct across filters table either — redundant with the
        #  per-category specificity table above and the §1 funnel-level table.)

        # ---- Aligned insert length (per BC, STAR): violin + overlaid density profile ----
        # ILLUMINA: removed. "Aligned insert length" is a READ length capped at 150 bp,
        # not the molecule length these plots were designed to show (axes run to 1600 bp).
        try:
            if IS_ILLUMINA:
                raise _IllRemoved('removed for Illumina: aligned length is a 150bp READ length, not molecule length')
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as _plt
            import numpy as _np
            try:
                from scipy.stats import gaussian_kde as _gkde
            except Exception:
                _gkde = None
            _lens = {}
            for bc in bc_list:
                h = (all_results[bc].get('star_summary', {}) or {}).get('aln_hist', {}) or {}
                arr = []
                for k, v in h.items():
                    k = int(k); v = int(v)
                    if v <= 0:
                        continue
                    # Spread the v reads uniformly across the 50-bp bin [k, k+50)
                    # (deterministic) rather than stacking them all at the bin centre.
                    # Stacking at centres made the KDE place a bump on every 50-bp
                    # centre → spurious ~50-bp "waviness" for peaky BCs. Mean/median of
                    # each bin are unchanged, so the violin markers are unaffected.
                    arr.extend(k + (j + 0.5) * (50.0 / v) for j in range(v))
                _lens[bc] = _np.array(arr)
            _cmap = _plt.get_cmap('tab10' if n_bc <= 10 else 'tab20')
            _order = [bc for bc in bc_list if len(_lens[bc])]

            _figv, _axv = _plt.subplots(figsize=(8, 3.2))
            _parts = _axv.violinplot([_lens[bc] for bc in _order], showmeans=True,
                                     showmedians=True, widths=0.85)
            for _b in _parts['bodies']:
                _b.set_facecolor('#4a90d9'); _b.set_alpha(0.55); _b.set_edgecolor('#274a6b')
            if _parts.get('cmeans') is not None:
                _parts['cmeans'].set_color('crimson'); _parts['cmeans'].set_linewidth(1.6)
            if _parts.get('cmedians') is not None:
                _parts['cmedians'].set_color('black')
            _axv.set_xticks(range(1, len(_order) + 1))
            _axv.set_xticklabels([bc_nick[bc] for bc in _order], fontsize=7)
            _axv.set_ylabel('STAR aligned insert length (bp)'); _axv.set_ylim(0, 1600)
            _axv.set_title('Aligned insert length per BC — violin (red=mean, black=median)',
                           fontsize=9, fontweight='bold')
            for _s in ('top', 'right'):
                _axv.spines[_s].set_visible(False)
            _plt.tight_layout()
            _ilv = os.path.join(fig_dir, 'star_insert_length_violin.png')
            _figv.savefig(_ilv, dpi=150, bbox_inches='tight'); _plt.close(_figv)
            elems.append(Image(_ilv, width=7 * inch, height=2.8 * inch))
            elems.append(Spacer(1, 6))

            _figp, _axp = _plt.subplots(figsize=(8, 3.0))
            _x = _np.linspace(0, 1500, 400)
            for _i, bc in enumerate(_order):
                L = _lens[bc]
                if _gkde is not None and len(L) > 5 and L.std() > 0:
                    _y = _gkde(L, bw_method=0.15)(_x) * 100 * 50
                else:
                    _c, _e = _np.histogram(L, bins=_np.arange(0, 1550, 50))
                    _y = _np.interp(_x, _e[:-1] + 25, _c / max(1, len(L)) * 100)
                _axp.plot(_x, _y, '-', lw=1.6, color=_cmap(_i % 20), label=bc_nick[bc])
                _axp.fill_between(_x, _y, color=_cmap(_i % 20), alpha=0.06)
            _axp.set_xlim(0, 1500); _axp.set_ylim(bottom=0)
            _axp.set_xlabel('STAR aligned insert length (bp)')
            _axp.set_ylabel('% of species-assigned reads per 50 bp')
            _axp.set_title('Aligned insert-length profile — all BCs overlaid',
                           fontsize=9, fontweight='bold')
            _axp.legend(fontsize=6.5, ncol=2, frameon=False)
            for _s in ('top', 'right'):
                _axp.spines[_s].set_visible(False)
            _plt.tight_layout()
            _ilo = os.path.join(fig_dir, 'star_insert_length_overlay.png')
            _figp.savefig(_ilo, dpi=150, bbox_inches='tight'); _plt.close(_figp)
            elems.append(Image(_ilo, width=7 * inch, height=2.6 * inch))
            elems.append(Spacer(1, 6))
        except Exception as _e:
            print(f'  (insert-length plots skipped: {_e})')

        # ---- RT extent: template-switch site to transcript 3' end ----
        # Reads only results['rt_extent'] (computed by compute_rt_extent_all from the BAM).
        try:
            _rte = {bc: all_results[bc].get('rt_extent') for bc in bc_list}
            _rte = {bc: v for bc, v in _rte.items() if v and v.get('n')}
            if _rte:
                import numpy as _np
                import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as _plt
                _cmap = _plt.get_cmap('tab10' if n_bc <= 10 else 'tab20')
                _figr, _axr = _plt.subplots(1, 2, figsize=(9.6, 3.4))
                for _i, bc in enumerate(_rte):
                    v = _rte[bc]
                    _h = _np.array(v['hist'], float) / v['n']; _x = _np.arange(len(_h)) * v['binw'] + v['binw'] / 2
                    _axr[0].plot(_x, _np.cumsum(_h) * 100, '-', lw=1.6, color=_cmap(_i % 20), label=bc_nick[bc])
                    _fh = _np.array(v['fhist'], float) / v['n']; _fx = _np.arange(len(_fh)) * v['fbinw'] + v['fbinw'] / 2
                    _axr[1].plot(_fx, _fh * 100, '-', lw=1.5, color=_cmap(_i % 20), label=bc_nick[bc])
                _axr[0].set_xscale('log'); _axr[0].set_xlim(10, 3000)
                _axr[0].set_xticks([10, 30, 100, 300, 1000, 3000]); _axr[0].set_xticklabels(['10', '30', '100', '300', '1000', '3000'])
                _axr[0].set_ylim(0, 100); _axr[0].set_ylabel('% of reads (cumulative)')
                _axr[0].set_xlabel("template-switch site to transcript 3' end (nt, log scale)")
                _axr[1].set_xlim(0, 400); _axr[1].set_ylim(bottom=0); _axr[1].set_ylabel('% of reads per bin')
                _axr[1].set_xlabel('same distance, 0-400 nt, 10-nt bins')
                _axr[0].legend(fontsize=6.5, frameon=False, loc='upper left')
                for _a in _axr:
                    for _s in ('top', 'right'):
                        _a.spines[_s].set_visible(False)
                    _a.grid(alpha=0.15, lw=0.5)
                _figr.suptitle("How far reverse transcription got: template-switch site to transcript 3' end",
                               fontsize=9, fontweight='bold')
                _plt.tight_layout()
                _rtp = os.path.join(fig_dir, 'star_rt_extent.png')
                _figr.savefig(_rtp, dpi=150, bbox_inches='tight'); _plt.close(_figr)
                elems.append(Image(_rtp, width=7 * inch, height=2.6 * inch))
                elems.append(Paragraph(
                    "<i>Reverse transcription starts at the poly(A) tail, so the distance from the template-switch site to the "
                    "transcript 3' end is the length of the RT product (the cDNA insert). Left: share of reads whose RT product is "
                    "at most that long; right: the 0-400 nt range in 10-nt bins. Reads counted: uniquely mapped reads inside a "
                    "named protein-coding gene (rRNA and mitochondrial reads excluded); the site is the aligned end nearest the "
                    "transcript 5' end, and the distance runs to the nearest annotated transcript end of that gene downstream of "
                    "it. If labelled sites stop RT the curves shift left with labelling density; if RT reads through they reflect "
                    "processivity and transcript length. Per barcode:</i>", body_style))
                _rtt = [['Barcode', 'reads scored', 'median (nt)', 'IQR (nt)', 'under 150 nt', 'over 2 kb']]
                for bc, v in _rte.items():
                    _rtt.append([bc_nick[bc], f"{v['n']:,}", f"{v['median']:,}", f"{v['q25']:,}-{v['q75']:,}",
                                 f"{v['pct_lt150']:.0f}%", f"{v['pct_gt2000']:.0f}%"])
                elems.append(make_table(_rtt, col_widths=[1.3 * inch] + [1.0 * inch] * 5, font_size=7))
                elems.append(Spacer(1, 6))
        except Exception as _e:
            print(f'  (RT-extent plot skipped: {_e})')

        # ---- Incorrect symmetric-ends (P5-P5 vs P7-P7) % vs insert length -----
        # Illumina-only, GATED: rendered only when the Illumina QC has run AND P5 and
        # P7 are each detected in >=10% of reads. Recomputed from the BAMs (per-read
        # P5/P7 adapter hits x STAR aligned insert length, pooled across BCs), cached
        # by a BAM size+mtime signature so PDF re-renders don't rescan the BAMs.
        try:
            _iqs = [all_results[bc].get('illumina_qc') for bc in bc_list]
            _iqs = [q for q in _iqs if q]
            _p5m = [q['partial_p5_pct'] for q in _iqs if q.get('partial_p5_pct') is not None]
            _p7m = [q['partial_p7_pct'] for q in _iqs if q.get('partial_p7_pct') is not None]
            _gate = (_p5m and _p7m and (sum(_p5m) / len(_p5m) >= 10.0)
                     and (sum(_p7m) / len(_p7m) >= 10.0))
            if _gate:
                import subprocess as _sp, re as _re, json as _json
                import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as _plt
                import illumina_adapter_qc_illumina as _IQ
                _CIGp = _re.compile(r'(\d+)([MIDN=X])')
                def _qaln(c):
                    return sum(int(x) for x, o in _CIGp.findall(c) if o in 'MI=X')
                _BINW = 25; _MAXX = 1000; _nb = _MAXX // _BINW
                _bams = [(bc, os.path.join(output_dir, 'individual-analyses', bc,
                                           f'{bc}_star_combined.bam')) for bc in bc_list]
                _bams = [(bc, p) for bc, p in _bams if os.path.exists(p)]
                _sig = ';'.join(f'{os.path.getsize(p)}:{int(os.path.getmtime(p))}'
                                for _, p in _bams)
                _cache = os.path.join(fig_dir, 'star_symm_ends_vs_insert_p5_p7.json')
                _bins = None
                if os.path.exists(_cache):
                    try:
                        _cj = _json.load(open(_cache))
                        if _cj.get('sig') == _sig:
                            _bins = _cj
                    except Exception:
                        _bins = None
                if _bins is None:
                    _cn = [0] * _nb; _c5 = [0] * _nb; _c7 = [0] * _nb; _cb = [0] * _nb
                    for _bc, _p in _bams:
                        _pr = _sp.Popen(['samtools', 'view', '-F', '0x904', _p],
                                        stdout=_sp.PIPE, text=True)
                        for _line in _pr.stdout:
                            _f = _line.split('\t'); _ch = _f[2]
                            if not (_ch.startswith('HUMAN_') or _ch.startswith('MOUSE_')):
                                continue
                            _b = _qaln(_f[5]) // _BINW
                            if _b >= _nb:
                                continue
                            _cn[_b] += 1
                            _s = _f[9]; _L = len(_s)
                            _h5 = _IQ._hits(_s, _IQ._P5spec); _h7 = _IQ._hits(_s, _IQ._P7spec)
                            _e5 = _IQ._both_ends_same(_h5, _L)
                            _e7 = _IQ._both_ends_same(_h7, _L)
                            if _e5:
                                _c5[_b] += 1
                            if _e7:
                                _c7[_b] += 1
                            if _e5 or _e7:
                                _cb[_b] += 1
                        _pr.wait()
                    _bins = {'sig': _sig, 'binw': _BINW, 'n': _cn, 'n5': _c5, 'n7': _c7, 'nboth': _cb}
                    try:
                        _json.dump(_bins, open(_cache, 'w'))
                    except Exception:
                        pass
                _cn = _bins['n']; _c5 = _bins['n5']; _c7 = _bins['n7']; _cb = _bins['nboth']
                _MINN = 40; _xs = []; _y5 = []; _y7 = []; _yb = []
                for _b in range(_nb):
                    if _cn[_b] >= _MINN:
                        _xs.append(_b * _BINW + _BINW / 2)
                        _y5.append(100 * _c5[_b] / _cn[_b])
                        _y7.append(100 * _c7[_b] / _cn[_b])
                        _yb.append(100 * _cb[_b] / _cn[_b])
                if _xs:
                    _figs, _axs = _plt.subplots(figsize=(8, 3.0))
                    _axs.plot(_xs, _y5, '-o', ms=3, color='#d62728', lw=1.5, label='P5-P5 (both ends P5)')
                    _axs.plot(_xs, _y7, '-o', ms=3, color='#1f77b4', lw=1.5, label='P7-P7 (both ends P7)')
                    _axs.plot(_xs, _yb, '--', color='#7f7f7f', lw=1.1, label='either (combined)')
                    _axs.set_xlim(0, _MAXX); _axs.set_ylim(0, max(max(_yb), 1) * 1.15)
                    _axs.set_xlabel('STAR aligned insert length (bp, 25-bp bins)')
                    _axs.set_ylabel('% reads, incorrect symmetric ends')
                    _axs.set_title('Incorrect symmetric ends (Illumina) vs insert length — '
                                   'P5-P5 vs P7-P7 (pooled)', fontsize=9, fontweight='bold')
                    _axs.legend(fontsize=6.5, frameon=False)
                    for _s in ('top', 'right'):
                        _axs.spines[_s].set_visible(False)
                    _plt.tight_layout()
                    _syp = os.path.join(fig_dir, 'star_symm_ends_vs_insert_p5_p7.png')
                    _figs.savefig(_syp, dpi=150, bbox_inches='tight'); _plt.close(_figs)
                    # Illumina adapter/index QC table — the rows moved out of the
                    # Master Summary. Built straight from illumina_qc (guaranteed
                    # present here, since the plot gate requires it) and placed on the
                    # SAME new page as the plot. _IQ imported above for the NEB match.
                    _il_moved = [
                        ('P5+P7 partial, opposite ends (Illumina %)',        'p5p7_opposite_ends_pct'),
                        ('i5+i7 index seqs, opposite ends (Illumina %)',     'idx_i5i7_opposite_ends_pct'),
                        ('TruSeq R1+R2 handles, opposite ends (Illumina %)', 'truseq_r1r2_opposite_ends_pct'),
                        ('Index detected (NEBNext Set-1)',                   '__neb__'),
                        ('P5->Read1 contiguous (correct %)',                 'p5_read1_pct'),
                        ('P5->Read2 contiguous (WRONG %)',                   'p5_read2_pct'),
                        ('P7->Read2 contiguous (correct %)',                 'p7_read2_pct'),
                        ('P7->Read1 contiguous (WRONG %)',                   'p7_read1_pct'),
                    ]
                    def _ilm_cell(_bc, _key):
                        _q = all_results[_bc].get('illumina_qc') or {}
                        if _key == '__neb__':
                            _i5 = _IQ.match_neb_index(_q.get('i5_dominant'), 'i5') or 'unknown'
                            _i7 = _IQ.match_neb_index(_q.get('i7_dominant'), 'i7') or 'unknown'
                            return f"i5={_i5} / i7={_i7}"
                        _v = _q.get(_key)
                        return '-' if _v is None else f"{_v:.1f}%"
                    _ilm_data = [['Illumina adapter / index QC'] + [bc_header[bc] for bc in bc_list]]
                    for _lbl, _key in _il_moved:
                        _ilm_data.append([_lbl] + [_ilm_cell(bc, _key) for bc in bc_list])
                    elems.append(PageBreak())
                    elems.append(KeepTogether([
                        Paragraph(
                            "<i>Molecules carrying the same Illumina adapter at both ends (P5-P5 or P7-P7), a library-construction "
                            "artefact, as a share of reads in each 25-bp bin of aligned insert length (bins with ≥40 reads, all "
                            "barcodes pooled; dashed grey = either). Shown when P5 and P7 are each found in ≥10% of reads.</i>", body_style),
                        Image(_syp, width=7 * inch, height=2.6 * inch),
                    ]))
                    elems.append(Spacer(1, 6))
        except Exception as _e:
            print(f'  (P5/P7 symm-vs-insert plot skipped: {_e})')

        # ---- Barcode % correct vs aligned read length (per BC, STAR) --------
        # ILLUMINA: removed (both the raw and the moving-average panel). Read length is
        # capped at 150 bp and binned by 50, so only ~3 bins can ever exist — the trend
        # line, and the Spearman that accompanies it, are not meaningful.
        try:
            if IS_ILLUMINA:
                raise _IllRemoved('removed for Illumina: read length capped at 150bp gives only ~3 bins')
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as _plt
            _MINN = 20   # min species-confident barcoded reads per 50-bp bin
            _figc, _axc = _plt.subplots(figsize=(8, 3.0))
            _cmapc = _plt.get_cmap('tab10' if n_bc <= 10 else 'tab20')
            _any_cvl = False; _ymin_all = 100.0
            # Prefer the distance-filtered length/accuracy (24plex full-molecule funnel:
            # reads that reached the barcode-insert distance step); fall back to the
            # pre-funnel star_summary population for chemistries without that funnel.
            def _cvl_src(bc):
                return ((all_results[bc].get('fullmol_funnel') or {}).get('correct_vs_len')
                        or (all_results[bc].get('star_summary', {}) or {}).get('correct_vs_len', {}) or {})
            _use_ff_cvl = any((all_results[bc].get('fullmol_funnel') or {}).get('correct_vs_len')
                              for bc in bc_list)
            _cvl_note = ('reads through the barcode-insert distance filter'
                         if _use_ff_cvl else _excl_desc + ' excluded')
            for _i, bc in enumerate(bc_list):
                cvl = _cvl_src(bc)
                pts = []
                for _k, _v in cvl.items():
                    _c, _t = _v
                    if _t >= _MINN:
                        pts.append((int(_k), _c / _t * 100))
                pts.sort()
                if len(pts) < 2:
                    continue
                _any_cvl = True
                _xs = [k + 25 for k, _ in pts]
                _ys = [y for _, y in pts]
                _axc.plot(_xs, _ys, '-', lw=0.9, alpha=0.45, color=_cmapc(_i % 20))
                _W = 5   # moving-average window (# of 50-bp bins ~ 250 bp), drawn over the bins
                _ysm = [sum(_ys[max(0, j - _W // 2): j + _W // 2 + 1]) / len(_ys[max(0, j - _W // 2): j + _W // 2 + 1]) for j in range(len(_ys))]
                _axc.plot(_xs, _ysm, '-', lw=2.0, color=_cmapc(_i % 20), label=bc_nick[bc])
                _ymin_all = min(_ymin_all, min(_ys))
            if _any_cvl:
                _axc.set_xlim(0, 1000)
                _axc.set_ylim(max(0.0, min(_ymin_all, 98.0) - 1.0), 100.2)
                _axc.set_xlabel('STAR aligned read length (bp, 50-bp bins)')
                _axc.set_ylabel('% barcoded reads correct')
                _axc.set_title('Barcode % correct vs aligned read length per BC (thin: 50-bp bins; thick: 5-bin moving average)',
                               fontsize=9, fontweight='bold')
                _axc.legend(fontsize=6.5, ncol=2, frameon=False)
                for _s in ('top', 'right'):
                    _axc.spines[_s].set_visible(False)
                _plt.tight_layout()
                _clp = os.path.join(fig_dir, 'star_correct_vs_length.png')
                _figc.savefig(_clp, dpi=150, bbox_inches='tight'); _plt.close(_figc)
                elems.append(KeepTogether([
                    Paragraph(
                        "<i>Per-BC barcode correctness vs STAR aligned read length (50-bp bins; "
                        "bins with ≥20 species-confident barcoded reads; " + _cvl_note + " — "
                        "see Methods). Thin line = per-bin value, thick line = 5-bin (250-bp) moving average. A downward "
                        "slope means longer reads barcode less accurately; Spearman rank correlation of accuracy with length: "
                        + ('; '.join(f"{bc_nick[bc]} ρ = {(all_results[bc].get('star_summary') or {}).get('length_rho'):.2f}" for bc in bc_list if (all_results[bc].get('star_summary') or {}).get('length_rho') is not None) or 'n/a')
                        + ".</i>", body_style),
                    Image(_clp, width=7 * inch, height=2.6 * inch),
                ]))
                elems.append(Spacer(1, 6))

                # ---- Same plot as a MOVING AVERAGE, y-axis 98-100% ----
            else:
                _plt.close(_figc)
        except Exception as _e:
            print(f'  (correct-vs-length plot skipped: {_e})')

    # ============ Library Saturation (exact, from the duplicity histogram) ============
    # The Master Summary already carries the scalar
    # complexity rows; this is the same quantity as a CURVE, which is what shows whether
    # a library is sampled out or still climbing. Computed exactly from the duplicity
    # histogram (no subsampling, no RNG) by umi_dedup.saturation_curve.
    try:
        import umi_dedup_illumina as _ud_sc
        _ud_sc.configure_umi(EXPERIMENT_CONFIG.get('rt_primers_used'))
        _sc_series = {}
        for _bc in bc_list:
            _d = (all_results[_bc].get('dup_funnel') or {})
            _h = _d.get('duplicity_hist') or _d.get('duplicity_hist_aligned')   # definition D first
            if _h:
                _pts = _ud_sc.saturation_curve(_h)
                if _pts:
                    _sc_series[_bc] = _pts
        if _sc_series:
            import matplotlib.pyplot as _plt
            _figs, _axs = _plt.subplots(figsize=(9.5, 3.6))
            for _bc, _pts in _sc_series.items():
                _x = [_p[0] / 1e6 for _p in _pts]
                _y = [_p[1] / 1e6 for _p in _pts]
                _lbl = f"{bc_nick.get(_bc, bc_short.get(_bc, _bc))} ({_y[-1]:.2f}M)"   # library name, as in every table
                _axs.plot(_x, _y, lw=1.8, marker='o', ms=2.6, label=_lbl)
            _axs.set_xlabel('reads sampled (millions)', fontsize=8)
            _axs.set_ylabel('distinct molecules (millions)', fontsize=8)
            _axs.set_title('Library saturation: molecules recovered vs reads sampled '
                           '(exact, from the duplicity histogram)',
                           fontsize=9, fontweight='bold')
            # The libraries were sequenced to DIFFERENT depths, so their right-hand end
            # points are not comparable to each other. Mark the shallowest end depth:
            # left of this line every library is being read at the same depth, which is
            # the only place a fair molecules-per-read comparison can be made.
            _xmin_end = min(_p[-1][0] for _p in _sc_series.values()) / 1e6
            if len(_sc_series) > 1:
                _axs.axvline(_xmin_end, ls='--', lw=1.0, color='#5E6B68', zorder=1)
                _axs.text(_xmin_end, _axs.get_ylim()[1] * 0.02, ' common depth',
                          fontsize=6.5, color='#5E6B68', va='bottom', ha='left')
            _axs.legend(fontsize=7, frameon=False)
            _axs.tick_params(labelsize=7)
            _axs.grid(True, ls=':', lw=0.5, color='#C8CFCD')
            _axs.set_axisbelow(True)
            for _sp in ('top', 'right'):
                _axs.spines[_sp].set_visible(False)
            _plt.tight_layout()
            _satp = os.path.join(fig_dir, 'library_saturation.png')
            _figs.savefig(_satp, dpi=150, bbox_inches='tight'); _plt.close(_figs)
            elems.append(PageBreak())
            elems.append(Paragraph(f"{_secn()}. Library Saturation", heading_style))
            elems.append(KeepTogether([
                Paragraph(
                    "<i>Distinct molecules recovered as a function of reads sampled, computed "
                    "EXACTLY from the observed duplicity histogram rather than by subsampling: "
                    "for a molecule seen k times in N reads, the chance that d sampled reads all "
                    "miss it is C(N−k,d)/C(N,d), evaluated in log-gamma space. A curve still "
                    "rising at its right-hand end means the library was NOT sequenced out and "
                    "more depth would still return new molecules; a curve that flattens means "
                    "the pool is exhausted, which for a negative control is expected and is a "
                    "statement about pool SIZE, not quality. End-point molecule counts are in "
                    "the legend. Same read scope as the duplicate rate.</i>",
                    body_style),
                Image(_satp, width=7.2 * inch, height=2.7 * inch),
            ]))
            elems.append(Spacer(1, 6))
    except Exception as _e:                                   # noqa: BLE001
        print(f'  (library saturation plot skipped: {_e})')

    # ============ 3. Transcript Diversity (STAR gene assignment) ============
    elems.append(PageBreak())
    elems.append(Paragraph(f"{_secn()}. Transcript Diversity", heading_style))
    _divs = {bc: (all_results[bc].get('star_transcript_diversity') or {}) for bc in bc_list}
    if any(d.get('total_mapped', 0) > 0 for d in _divs.values()):
        _ug = [d.get('unique_genes', 0) for d in _divs.values() if d.get('total_mapped')]
        _ug_rng = f"{min(_ug):,}" + (f"–{max(_ug):,}" if max(_ug) != min(_ug) else "")
        _si = [d.get('simpson', 0) for d in _divs.values() if d.get('total_mapped')]
        elems.append(Paragraph(
            "<i>Diversity over <b>protein-coding mRNAs only</b> (ribosomal-protein and "
            "mitochondrial genes excluded): each mapped read is assigned to its max-overlap "
            "gene from the genome annotation; intronic reads inside a protein-coding gene and reads whose "
            "alignments hit both genomes count here although they do not in the composition table's mRNA "
            "row, so this universe is larger than that row. "
                f"Unique mRNAs {_ug_rng}. A high share of reads in the top genes flags a low-complexity library.</i>",
            body_style))
        def _dv(r):
            return r.get('star_transcript_diversity') or {}
        def _mrna_sp_pct(r, sp):
            # per-species protein-coding gene-assigned reads (rRNA-masked), as % of all
            # mapped reads. Uses the SAME gene-diversity basis as the 'mRNA reads' total
            # and the unique-gene rows (top_genes[sp]['total']), so the two sub-rows sum
            # to the total row. (Every combined-genome gene is HUMAN_/MOUSE_, so
            # human+mouse == total_mapped exactly.) rRNA on rDNA-in-gene loci is already
            # excluded upstream in _star_gene_stats (see its rRNA MASK note).
            tx = r.get('star_taxonomy') or {}
            mapped = sum(sum(x for x in v.values() if isinstance(x, (int, float)))
                         for k, v in tx.items()
                         if k != 'Unmapped / contaminant' and isinstance(v, dict))
            c = ((r.get('star_top_genes') or {}).get(sp) or {}).get('total', 0) or 0
            return f"{c:,} ({100 * c / mapped:.1f}%)" if mapped else f"{c:,}"
        _dmetrics = [
            ('mRNA reads (PC, non-ribo/MT)', lambda r: f"{_dv(r).get('total_mapped', 0):,}"),
            ('Unique mRNAs', lambda r: f"{_dv(r).get('unique_genes', 0):,}"),
            # Split by the GENE's species (chrom prefix), not the read's barcode. These
            # can sum to more than 'Unique mRNAs' above, which pools by gene name and so
            # merges orthologues that share a name (ACTB/Actb).
            ('→ unique human genes', lambda r: f"{_dv(r).get('unique_genes_human', 0):,}"),
            ('→ unique mouse genes', lambda r: f"{_dv(r).get('unique_genes_mouse', 0):,}"),
            ('Top 1 mRNA %', lambda r: f"{_dv(r).get('top1_pct', 0)}%"),
            ('Top 5 mRNAs %', lambda r: f"{_dv(r).get('top5_pct', 0)}%"),
            ('Top 10 mRNAs %', lambda r: f"{_dv(r).get('top10_pct', 0)}%"),
        ]
        _ddata = [['Metric'] + [bc_header[bc] for bc in bc_list]]
        for _lab, _fn in _dmetrics:
            _ddata.append([_lab] + [_fn(all_results[bc]) for bc in bc_list])
        elems.append(make_heatmap_table(_ddata, col_widths=auto_col_widths(120, n_bc)))
        elems.append(Spacer(1, 8))
        if any(all_results[bc].get('star_top_genes') for bc in bc_list):
            _tg_small = ParagraphStyle('TGsmall', parent=body_style, fontSize=7, leading=8)

            def _topgene_table(sp_key, sp_label):
                _data = [['#'] + [bc_header[bc] for bc in bc_list]]
                for _rank in range(5):
                    _row = [str(_rank + 1)]
                    for bc in bc_list:
                        _sp = (all_results[bc].get('star_top_genes') or {}).get(sp_key) or {}
                        _top = _sp.get('top') or []
                        _tot = _sp.get('total') or 1
                        if _rank < len(_top):
                            _gn, _gc = _top[_rank][:2]
                            _mm = _top[_rank][2] if len(_top[_rank]) > 2 else 0
                            _row.append(Paragraph(f"{_gn}{'†' if _mm >= 0.9 else ''} ({_gc / _tot * 100:.1f}%)", _tg_small))
                        else:
                            _row.append('-')
                    _data.append(_row)
                _cap = Paragraph(
                    f"<b>Top 5 {sp_label} protein-coding mRNAs by barcode</b> "
                    f"<i>(non-ribosomal, non-mitochondrial; % of {sp_label} mRNA reads). "
                    "† ≥90% of this gene's reads are STAR multimappers assigned by their primary alignment: "
                    "a multi-copy or repeat-like locus, not uniquely attributable expression.</i>",
                    body_style)
                return KeepTogether([_cap,
                                     make_table(_data, col_widths=auto_col_widths(20, n_bc)),
                                     Spacer(1, 6)])

            elems.append(_topgene_table('human', 'human'))
            elems.append(_topgene_table('mouse', 'mouse'))
    else:
        elems.append(Paragraph("<i>No STAR gene assignments available.</i>", body_style))

    # ============ 4. Sample PCA (correct vs chemistry-swap, STAR) ============
    # The 'Sample PCA — correct vs chemistry-swapped reads' section is not rendered:
    # the swap pseudo-samples are the wrong-barcode reads, tens to a few hundred per
    # barcode on any good run, so the decomposition is noise exactly when the run is good. The
    # per-gene counts (star_pca_gene_counts) stay in the JSON.

    # ILLUMINA: section removed. It tallies reads containing each ONT/BOB reference
    # sequence (adapters, TSO backbone, polydT primer, ONT motifs) searched in BOTH
    # orientations. On Illumina R2 those handles are sequencing primers that are not
    # in the read, the orientation is fixed, and short motifs hit by chance — so the
    # table reports either structural zeros or noise.
    if not IS_ILLUMINA:
        elems.append(PageBreak())
        elems.append(Paragraph(f"{_secn()}. Sequence Reference Detection", heading_style))
        elems.append(Paragraph(
            "<i>Reads containing each element of this run's construct (exact match, either orientation), as a share of all reads. "
            "A read can contain several elements, so the rows do not add up. Short elements pick up chance matches. "
            "Elements found in no read are omitted; the full oligo list lives in sequence_list.txt.</i>",
            body_style))
        # On 24plex runs only the elements of this construct are listed; the
        # legacy oligo list stays in sequence_list.txt for the legacy chemistry paths.
        _is_24_seq = _ba.is_24plex(EXPERIMENT_CONFIG.get('tso_species_map'), EXPERIMENT_CONFIG.get('tso_arch'))
        _seq_desc = {name: desc for sec in SEQUENCE_SECTIONS.values() for (name, _q, desc, _u) in sec}
        _seq_allowed = ({name for sec in ('24plex Construct', 'ONT Adapters') for (name, _q, _d, _u) in SEQUENCE_SECTIONS.get(sec, [])}
                        if _is_24_seq else None)
        seq_data = [['Sequence', 'What it is'] + [bc_header[bc] for bc in bc_list]]
        seq_pcts = {}  # (row, col) -> pct for gradient coloring
        _out_ri = 0    # displayed-row index (0-zero rows are dropped, so it != enumerate)
        _n_omitted = 0
        for (label, _seq) in SEQUENCE_REFERENCE:
            _nm = label.rsplit(' (', 1)[0]
            if (_seq_allowed is not None and _nm not in _seq_allowed) or _nm == 'ONT_native_adapter':   # trimmed reads never hold the full adapter
                continue
            row = [label, Paragraph(_seq_desc.get(_nm, ''), ParagraphStyle('sqd', parent=body_style, fontSize=6.5, leading=8))]
            cell_pcts = []
            row_total = 0
            for bc in bc_list:
                r = all_results[bc]
                count = r.get('seq_ref_counts', {}).get(label, 0)
                total_r = r.get('total_reads', 1)
                pct = count / max(total_r, 1) * 100
                cell_pcts.append(pct)
                row_total += count
                row.append(f"{count} ({pct:.0f}%)")
            if row_total == 0:
                _n_omitted += 1          # not found in any read of any barcode — omit row
                continue
            _out_ri += 1
            for ci, pct in enumerate(cell_pcts):
                seq_pcts[(_out_ri, ci + 2)] = pct  # +1 header row; two label columns
            seq_data.append(row)
        if len(seq_data) == 1:           # everything omitted (no reference sequence detected)
            seq_data.append(['(no reference sequence detected in any barcode)', '']
                            + [''] * n_bc)

        cw = [150, 190] + [(usable_w - 340) / max(n_bc, 1)] * n_bc
        seq_t = Table(seq_data, colWidths=cw)
        seq_cmds = [
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('FONTSIZE', (0, 0), (-1, -1), 7),
            ('TOPPADDING', (0, 0), (-1, -1), 2),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 2),
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#4472C4')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('ALIGN', (2, 0), (-1, -1), 'CENTER'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ]
        # Gradient: 0% = white, 50% = yellow, 100% = green
        def pct_to_gradient(pct):
            pct = min(max(pct, 0), 100)
            if pct <= 50:
                # White (255,255,255) → Light Yellow (255,255,180)
                t = pct / 50.0
                return colors.Color(1.0, 1.0, 1.0 - t * (75/255))
            else:
                # Light Yellow (255,255,180) → Light Green (180,230,180)
                t = (pct - 50) / 50.0
                return colors.Color(1.0 - t * (75/255), 1.0 - t * (25/255), 180/255)

        for (ri, ci), pct in seq_pcts.items():
            bg = pct_to_gradient(pct)
            seq_cmds.append(('BACKGROUND', (ci, ri), (ci, ri), bg))
        seq_t.setStyle(TableStyle(seq_cmds))
        elems.append(seq_t)
        elems.append(Spacer(1, 8))

        # 6. Methods & Settings
    # ILLUMINA: the ONT Methods section is not rendered (reasons below), but the report
    # still needs a Methods section for its "(see Methods)" cross-references. This is a
    # compact replacement carrying only what applies to this path, and every value is
    # DERIVED from the run rather than written down, so it cannot drift out of step
    # with the analysis (see METHODS PROVENANCE above).
    if IS_ILLUMINA:
        elems.append(PageBreak())
        elems.append(Paragraph(f"{_secn()}. Methods &amp; Settings", heading_style))
        _mst = ParagraphStyle('MethodsI', parent=body_style, fontSize=7.5,
                              leading=10, spaceAfter=4)
        _msb = ParagraphStyle('MethodsIBold', parent=_mst, fontName='Helvetica-Bold',
                              fontSize=8, spaceBefore=6)

        def _mp(txt, style=None):
            elems.append(Paragraph(txt, style or _mst))

        _smap = EXPERIMENT_CONFIG.get('primary_species_map') or {}
        _smap_s = ', '.join(f'{k}={v}' for k, v in sorted(_smap.items())) or 'not declared'
        _offs = EXPERIMENT_CONFIG.get('grun_offsets') or {}
        _offs_s = ', '.join(f'{k}+{v}' for k, v in sorted(_offs.items())) or 'not declared'
        _unmeas = [x.strip() for x in
                   str(EXPERIMENT_CONFIG.get('unmeasurable', '') or '').split(',') if x.strip()]

        _mp('Sequencing &amp; alignment', _msb)
        _mp(f'Sequencer: {METHODS_SEQUENCER}. Aligner: {METHODS_ALIGNER} with '
            f'{METHODS_ALIGNER_TUNING}, against a combined human+mouse genome index '
            f'(GRCh38.113 + GRCm39.113, chromosomes prefixed HUMAN_/MOUSE_). Species is '
            f'assigned by chromosome prefix, so both genomes must share one index; rDNA '
            f'loci are excluded because rDNA is asymmetric between the assemblies and '
            f'would otherwise force rRNA reads to human.')

        _mp('Species barcode', _msb)
        _mp(f'Read 2 begins at the bobcode because the TruSeq2 handle IS the read-2 '
            f'sequencing primer and is therefore absent from the read. The 7-mer is '
            f'called POSITIONALLY at R2[0:7] (<= 1 mismatch), not by searching for an '
            f'upstream anchor. Species map for this run: {_smap_s}. Non-templated G-run '
            f'offsets from the 7-mer end: {_offs_s}. The map is an expectation, not a '
            f'fact: it is re-derived from the alignment and reported in '
            f'species_map_check.txt, which is the authority.')

        _mp('Duplicates', _msb)
        _mp(f'Definition D: two reads are one molecule when they share bobcode, contig, strand, 5′ anchor '
            f'position{DEDUP_KEY_EXTRA} and their UMIs cluster ({METHODS_DUPLICATE}). The duplicate-rate, complexity and '
            f'saturation rows are over ALIGNED reads of the PRE-deduplication population (unmapped reads are '
            f'deduplicated too but not scored), so they describe the library as sequenced rather than the '
            f'filtered set the rest of the report scores.')

        _mp('Barcode correctness', _msb)
        _mp('Barcode accuracy is the fraction of reads whose bobcode species matches the '
            'species of the genome the read aligned to. It is reported cumulatively down '
            'a funnel; the mRNA-only level is the meaningful one, because ribosomal RNA '
            'is deeply conserved between the two genomes and cross-maps, which '
            'manufactures apparent error unrelated to barcoding.')

        if _unmeas:
            _mp('Not measurable at this read length', _msb)
            _UNMEAS_LABEL = {'polyt': 'the poly(dT) run', 'polya': 'the poly(A) tail', 'ge18': 'the 18-nt minimum insert',
                             'truseq1': 'the TruSeq1 handle', 'truseq2': 'the TruSeq2 handle', 'insert': 'the insert length',
                             'dist': 'the barcode-to-insert distance', 'grun': 'the G-run', 'c28': 'the C28 primer'}
            _mp('This run\'s reads cannot observe: <b>' + ', '.join(_UNMEAS_LABEL.get(x.lower(), x) for x in _unmeas) + '</b>. '
                'Those steps are dropped from the structure funnel rather than retained '
                'at zero, which would read as a total structural failure instead of an '
                'inapplicable filter. Declared in the run config and cross-checked '
                'against the measured read lengths before the run starts.')

    # ILLUMINA: '. Methods & Settings' section removed. Its prose documents ONT
    # constructs, adapter/backbone detection and long-read filters that no longer
    # apply, and every table it renders was already dropped on this path.
    if not IS_ILLUMINA:
        elems.append(PageBreak())
        elems.append(Paragraph(f"{_secn()}. Methods & Settings", heading_style))

        methods_style = ParagraphStyle('Methods', parent=body_style, fontSize=7.5,
                                        leading=10, spaceAfter=4)
        methods_bold = ParagraphStyle('MethodsBold', parent=methods_style,
                                       fontName='Helvetica-Bold', fontSize=8, spaceBefore=6)

        elems.append(Paragraph("Analysis Pipeline", methods_bold))
        elems.append(Paragraph(
            ("Each barcode is processed in three stages: (1) barcode and structure — every read is scanned on both "
             "strands for the TruSeq-R1 anchor + a declared code, TSO units are counted, and the full-molecule funnel "
             "grades structure from the FASTQ (poly-dT, TruSeq1/TruSeq2 handles) and from the STAR alignment (insert, "
             "barcode-insert distance); (2) species + RNA-class assignment — align full reads with "
             f"{METHODS_ALIGNER_BIN} to the combined human+mouse genome and assign species by chromosome prefix (cross-species "
             "multimappers = ambiguous), molecular category from the combined GTF, and rRNA from the rDNA loci learned from the alignments "
             "loci; (3) report generation — this cross-barcode summary (plus per-sample tables on pooled runs).") if _is_24 else
            ("Each barcode is processed in three stages: (1) read structure analysis — detect primer, polyT/A, "
             "and TSO components by exact sequence match in both orientations, classify read structures, extract "
             "cDNA inserts, detect PCR duplicates; (2) species + RNA-class assignment — align full reads with "
             f"{METHODS_ALIGNER_BIN} to the combined human+mouse genome and assign species by chromosome prefix (cross-species "
             "multimappers = ambiguous), molecular category from the combined GTF, and rRNA from the rDNA loci learned from the alignments "
             "loci; (3) report generation — individual barcode PDF and cross-barcode summary."),
            methods_style))

        elems.append(Paragraph("Alignment, Species &amp; RNA-class Assignment (STAR)", methods_bold))
        elems.append(Paragraph(
            f"Aligner: {METHODS_ALIGNER_BIN} (official 2.7.11b build) against a combined human+mouse <i>genome</i> index "
            f"(GRCh38.113 + GRCm39.113, chromosomes prefixed HUMAN_/MOUSE_) with {METHODS_ALIGNER_TUNING}. "
            "Species = chromosome prefix of a read's alignments: a read on one species only is assigned to it; "
            "a read multimapping across <b>both</b> species is <b>ambiguous</b>. Molecular category is assigned "
            "from the combined GTF (gene_biotype + exon overlap: protein-coding / ribosomal-protein / other "
            "ncRNA / pre-mRNA-intronic / intergenic / mitochondrial). rRNA reads are identified via the "
            "rDNA loci (rdna_loci.bed) and <b>excluded</b> from the species-barcode swap metric: the rDNA repeat "
            "is represented in GRCh38 but largely absent from GRCm39, so on the combined genome STAR otherwise "
            "force-types rRNA reads to human. Composition counts are over mapped reads.",
            methods_style))

        elems.append(Paragraph("Barcode Correctness (species-barcode swap)", methods_bold))
        elems.append(Paragraph(
            "<b>Barcode % correct</b> = the fraction of <i>barcoded</i> reads whose STAR-assigned "
            "species matches the species encoded by their barcode (" + _codes_desc + "). "
            "A read is scored only if it passes ALL of: (1) a species barcode is detected — exact, "
            "or " + _fuzzy_desc + "; (2) it aligns to exactly one species (reads multimapping across "
            "both species are <b>ambiguous</b> and excluded); (3) its molecular category is not an "
            "excluded class; (4) its STAR-aligned bases are not low-complexity (see below)"
            + ("; and (5) it carries exactly <b>one TSO</b> — reads with ≥2 TSO units are TSO–TSO "
               "chimeras bearing two conflicting barcodes (usually two different samples) and are "
               "excluded (see the TSO Structure &amp; Multiplicity section)" if _is_tso else "")
            + ". "
            "<b>Excluded categories:</b> (i) <b>rRNA</b> — rDNA is in GRCh38 but "
            "~absent from GRCm39, so rRNA reads are force-typed human and cannot be species-resolved; "
            "(ii) <b>Mitochondrial</b> — MT transcripts are excluded from the correctness denominator; "
            "(iii) <b>Intergenic / genomic</b> — dominated by insert-less / internally poly-dT-primed "
            "molecules whose polyT tail maps to genomic poly-A/T tracts present in BOTH species, "
            "making the species call arbitrary (verified by per-read inspection: such reads carry "
            "EXACT barcodes yet align with zero mismatches to genomic poly-A stretches); "
            "(iv) <b>unmapped</b> reads have no species (excluded implicitly). "
            "<b>Low-complexity filter:</b> a read is additionally dropped if the bases STAR "
            "actually aligned (the sequence that decides the species call) have a normalized "
            "dinucleotide Shannon entropy below 0.65 — i.e. the call rests on a homopolymer / "
            "simple-repeat (e.g. a polyT tail mapping to a genomic poly-A/T tract), regardless of "
            "read length. This is a sequence-complexity test, not a length cutoff; the 0.65 "
            "threshold was set on the development runs, where it removed about 70% of suspect "
            "short-swap reads while dropping about 0.1% of confidently mapped long reads. The identical "
            "exclusions apply to Ms/Hu % correct, the funnel-level table, and the "
            "barcode-correctness-vs-length plot, so every correctness number in this report is "
            "computed the same way. "
            ""
            "<b>Master-Summary barcode accuracy:</b> one row reports barcode % correct (barcode species vs "
            "STAR-aligned species, reads scored in parentheses) at the deepest barcode-core level of this chemistry; "
            "the full funnel (" + _ms_step_desc + " …) with accuracy at every level is in the Read Structure section. "
            "<b>Mouse/Human % correct</b> split the last barcode level "
            "by aligned species (of reads aligning to that species, the fraction correctly barcoded). "
            "<b>% reads ending in ONT adapter</b> = reads whose last 110 nt contain the ONT adapter block (the ONT@end rule of the read-composition section).",
            methods_style))

        elems.append(Paragraph("Read Structure Detection", methods_bold))
        elems.append(Paragraph(
            ("The component detector of the earlier chemistries (RT primer, polyT/A ≥10, 18-bp TSO primer, Read1-f; gaps ≥30 bp = insert; "
             "first + last 20 bp duplicate fingerprint) still runs, but on this chemistry its outputs are diagnostic files only: "
             "every structure row in this report comes from the full-molecule funnel and the STAR alignment "
             "(Read Structure Classification section), and duplicates from the gene → UMI → position → strand funnel.") if _is_24 else
            ("Components detected: RT primer (RT_HR), polyT/A (≥10 consecutive bases), TSO primer (18bp), "
             "and Read1-f (Illumina adapter). Searched in both forward and reverse complement. "
             "Components sorted by position, adjacent same-type deduplicated. Gaps ≥30bp classified as cDNA insert. "
             "Structure strings normalized: canonical form = alphabetically first of forward/reverse."),
            methods_style))

        elems.append(Paragraph("Species Barcode Extraction", methods_bold))
        if _is_24:
            elems.append(Paragraph(
                "The species code sits in the TSO immediately 3′ of the TruSeq-R1 anchor CTCTTCCGATCT (the end "
                "of the TruSeq2 handle), before the 6N UMI + spacer and the GGG template-switch; detected on both "
                "strands. It is assigned to the nearest declared code within ≤1 mismatch, longest code first, ties "
                "rejected (codes: " + _codes_desc + "). The declared map is re-derived from the alignment "
                "(species_map_check.txt), which is the authority.",
                methods_style))
        elif _is_tso:
            elems.append(Paragraph(
                "The species 7-mer sits in the TSO, immediately after the 20bp Nextera-R1 backbone "
                "(GTGACTGGAGTTCAGACGTG) and before the GGG template-switch; detected in both "
                "orientations. The 7-mer is assigned to the nearest of ATCGAAA (human) / CAGTTGA "
                "(mouse) within ≤2 mismatches — the two codes are Hamming-6 apart, so the call is "
                "unambiguous. Species ground truth for the swap metric comes from this 7-mer.",
                methods_style))
        else:
            elems.append(Paragraph(
                "The species 6-mer sits between the bob-c primer (ending in the ATGTC anchor) and "
                "the polyT; detected in both orientations. The 6-mer is assigned to the nearest of "
                "ATCGCT (human) / TCGATA (mouse) within ≤2 mismatches — the two codes are 6 "
                "mismatches apart, so the call is unambiguous.",
                methods_style))

        elems.append(Paragraph("RT Primer Subtypes", methods_bold))
        elems.append(Paragraph(
            "RT_HR detected = ATACTCGTGAC (11bp). "
            "RT_w_c28_6T = RT_HR + linker (ATCAGCTTCC) + T + 5×T (polyT priming). "
            "RT_w_c28_6N = RT_HR + linker + T + 6 bases (not all T; random N9 priming).",
            methods_style))

        elems.append(Paragraph("PCR Duplicate Detection (confident-duplicate funnel)", methods_bold))
        elems.append(Paragraph(
            "PCR duplicates are called with a stepwise confident-duplicate funnel over protein-coding, "
            "non-ribosomal mRNA reads carrying a callable 6 nt TSO UMI (the assessable universe; reads "
            "without a UMI cannot be scored). Reads are progressively collapsed into original molecules "
            "by adding one identity criterion at a time — (1) same gene; (2) + identical 6 nt UMI; "
            "(3) + 5′-start and 3′-end reference coordinates agreeing within ±10 bp (a tolerance that "
            "absorbs ONT indel/soft-clip jitter, so true copies are not split); (4) + same strand — and "
            "the duplicate fraction after the last step is the reported PCR duplicate rate (the intermediate steps are kept in the run's JSON). Strand is the discriminator: the library "
            "is double-stranded and each amplicon can be read from either strand, so genuine PCR-"
            "amplification duplicates (independent copies) split ~50/50 same/opposite strand, whereas the "
            "two complementary strands of one molecule are obligately opposite-strand. Step 3 therefore "
            "measures total molecular redundancy (both-strand sequencing + PCR); the final same-strand "
            "step isolates the PCR-amplification duplicate rate. Duplicates are counted, not removed.",
            methods_style))

        elems.append(Paragraph("Sequence Matching", methods_bold))
        elems.append(Paragraph(
            "Sequence reference detection (Sequence Reference Detection section) uses exact matching (forward and reverse complement). "
            "The species barcode is matched to the nearest declared code within 1 mismatch (see Species "
            "Barcode Extraction).",
            methods_style))

        elems.append(Paragraph("Color Coding", methods_bold))
        elems.append(Paragraph(
            "Composition table: fixed scale white 0% → yellow 10% → green ≥20% (the assignability table is not shaded). Barcode-specificity "
            "tables: magenta 0 → grey 50 → green 100%. Funnel tables: magenta (low retention) → yellow (~50%) → green. "
            "Older cross-barcode tables: white (0%) → light yellow (midpoint) → light green (maximum value). "
            "Sequence reference detection: fixed white (0%) → light yellow (50%) → light green (100%).",
            methods_style))

        # ILLUMINA: appendix 'Sequence Reference (from analysis guide ground truth)'
        # removed. It lists every ONT/BOB reference sequence (TSO backbone, RT/polydT
        # primer, C28 linker, ONT adapters, splint oligos) as the search vocabulary for
        # detection tables that no longer exist on this path — none of those handles are
        # present in an Illumina read.
        if not IS_ILLUMINA:
            # Comprehensive sequence reference table grouped by category
            elems.append(Paragraph("Sequence Reference (from analysis guide ground truth)", methods_bold))
            elems.append(Paragraph(
                "The elements of this run's construct, from sequence_list.txt (the full oligo list stays there for the "
                "earlier chemistries). Each is searched in both orientations for the Sequence Reference Detection section.",
                methods_style))
            seq_ref_data = [['Category', 'Name', 'Sequence', 'Length', 'What it is']]
            # Categorize ground truth entries by name prefix
            def categorize_gt(name):
                n = name.lower()
                if n.startswith('tso'): return 'TSO'
                if n in ('bob_hr-rt_hr',) or n.startswith('splint'): return 'Ligation'
                if n.startswith('bob') or n.startswith('bob-'): return 'BOB'
                if any(n.startswith(p) for p in ['rt_', 'polydt', 'c28']): return 'RT'
                if any(n.startswith(p) for p in ['read', 'adaptor']): return 'Illumina'
                if any(n.startswith(p) for p in ['truseq', 'r1_anchor', 'r2_full', 'p5', 'p7']): return 'Illumina'
                if n.startswith('ont'): return 'ONT'
                return 'Other'
            seq_ref_entries = []
            # Every entry of sequence_list.txt (the full list is kept in the report on purpose),
            # this run's construct elements first.
            for gt_name, gt_seq in sorted(GROUND_TRUTH.items(), key=lambda kv: 0 if (_seq_allowed and kv[0] in _seq_allowed) else 1):
                cat = categorize_gt(gt_name)
                seq_ref_entries.append((cat, gt_name, gt_seq, f'{len(gt_seq)}bp'))
            for cat, name, seq, length in seq_ref_entries:
                seq_ref_data.append([cat, name, seq, length, Paragraph(_seq_desc.get(name, ''), ParagraphStyle('sqd2', parent=methods_style, fontSize=6.5, leading=8))])
            seq_t = Table(seq_ref_data, colWidths=[50, 100, 215, 35, 170])
            seq_cmds = [
                ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
                ('FONTSIZE', (0, 0), (-1, -1), 6.5),
                ('FONTNAME', (2, 1), (2, -1), 'Courier'),
                ('TOPPADDING', (0, 0), (-1, -1), 2),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 2),
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#4472C4')),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('ALIGN', (3, 0), (3, -1), 'CENTER'),
            ]
            # Group separators: thick lines between categories
            prev_cat = None
            for ri, (cat, _, _, _) in enumerate(seq_ref_entries):
                if prev_cat and cat != prev_cat:
                    seq_cmds.append(('LINEABOVE', (0, ri + 1), (-1, ri + 1), 1.5, colors.black))
                prev_cat = cat
            seq_t.setStyle(TableStyle(seq_cmds))
            elems.append(seq_t)
            elems.append(Spacer(1, 8))

        elems.append(Paragraph("Key Parameters", methods_bold))
        params_data = [
            ['Parameter', 'Value'],
            ['Alignment tool', f'{METHODS_ALIGNER} (combined genome)'],
            ['Reference genome', 'combined GRCh38.113 + GRCm39.113 (HUMAN_/MOUSE_)'],
            ['Annotation', 'combined_genome.gtf (gene_biotype + exons)'],
            ['rRNA loci', 'rdna_loci.bed (learned from the alignments; excluded from swap)'],
            ['Species rule', 'chrom prefix; cross-species multimap = ambiguous'],
            ] + ([
                ['Species code call', 'TruSeq-R1 anchor CTCTTCCGATCT + declared code, ≤1 mismatch, longest code first, ties rejected; both strands'],
                ['TSO unit', 'anchor + code; hits within 27 nt on one strand collapse; ≥2 units = chimera (excluded from correctness)'],
                ['Full-molecule funnel', 'exact single unit → protein-coding (incl. RP) → barcode-insert distance ≤25 nt → STAR insert >30 bp → poly-dT ≥10 → ≥18 nt 5′ of it → TruSeq1 partial (anchor + ≥5/10 of CTACACGACG, 18-40 nt before the poly-dT) → TruSeq2 partial (anchor + ≥5/9 of AGACGTGTG)'],
                ['Full molecule', 'every funnel step above present; on 24plex through the TruSeq2 partial handle'],
                ['G-run', 'G count at the code end + set offset (A +8, B +11, C +14), ±1 nt; mapped barcoded reads'],
                ['Low-complexity filter', 'normalised dinucleotide entropy of the aligned bases < 0.65 (correctness rows only)'],
                ['Duplicate rule (counted, not removed)', 'same gene → exact 6-nt UMI → start/end within ±10 bp (UMI ed ≤1) → same strand'],
            ] if _is_24 else [
                ['PolyT/A threshold', '≥10 consecutive bases'],
                ['Insert gap threshold', '≥30bp between components'],
                ['Duplicate fingerprint', 'First + last 20bp of insert'],
            ]) + [
            ['Sequencer', METHODS_SEQUENCER],
        ]
        # long values (funnel rules) must wrap inside the 280-pt cell
        params_data = [params_data[0]] + [[k, Paragraph(str(v), methods_style)] for k, v in params_data[1:]]
        params_t = Table(params_data, colWidths=[160, 280])
        params_cmds = [
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('FONTSIZE', (0, 0), (-1, -1), 7),
            ('TOPPADDING', (0, 0), (-1, -1), 2),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 2),
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#4472C4')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ]
        params_t.setStyle(TableStyle(params_cmds))
        elems.append(params_t)

    # ===================================================================
    # Master summary — 21 per-BC headline metrics, inserted at the top of
    # the PDF (page 1). Species/swap rows are pure-STAR (results['star_summary']
    # / ['star_swap']); structure/adapter rows come from the read-structure
    # analysis already in each JSON.
    # ===================================================================
    def _g_pct_at_dn(bc, pos_idx):
        ff = (all_results[bc].get('flank_freq', {}) or {}).get('TSO-seq1-variant') or {}
        dn = ff.get('downstream') or []
        if pos_idx >= len(dn):
            return None
        d = dn[pos_idx]; tot = sum(d.values())
        return (d.get('G', 0) / tot * 100) if tot else None

    def _ms_vals(bc):
        r = all_results[bc]
        ss = r.get('star_summary', {}) or {}
        total = r.get('total_reads', 0) or 0
        tx = r.get('star_taxonomy', {}) or {}
        # RNA-class composition is over MAPPED reads (identical denominator to the
        # §3 composition table's _star_mapped), so the two tables agree.
        _mapped = sum(sum((v or {}).values()) for v in tx.values() if isinstance(v, dict))
        def _catpct(cat):
            n = sum((tx.get(cat, {}) or {}).values())   # all species in the category
            return (n / _mapped * 100) if _mapped else None
        _sd = (r.get('star_duplicates') or r.get('duplicates') or {})
        n_dup = _sd.get('exact_dup_reads', 0) or 0
        _sd_mr = _sd.get('mrna') or {}
        n_dup_mrna = _sd_mr.get('exact_dup_reads', 0) or 0
        n_mrna_fp = _sd_mr.get('n', 0) or 0
        canon = (r.get('structure_insert_stats', {}) or {}).get('TSO > insert > bob-c') or {}
        def _rate(k):
            return (r.get(k, 0) / total * 100) if total else None
        return {
            'star_taxonomy': tx,          # for the '% reads mapped' row
            'pct_correct': ss.get('pct_correct'),
            'pct_correct_retained': ss.get('pct_correct_retained'),
            'terminal_correct': ss.get('terminal_correct'),
            'terminal_pct': ss.get('terminal_pct'),
            'ms_spec': ss.get('ms_spec'), 'hu_spec': ss.get('hu_spec'),
            'total_reads': total,   # needed by _fm fallback (tso-bob) for the Reads(total) row
            'unified_funnel': _uf(bc),   # chemistry-adaptive per-step funnel (accuracy rows)
            'filter_funnel': r.get('filter_funnel'),
            'star_funnel_spec': r.get('star_funnel_spec'),   # _fm fallback (tso-bob) Mouse/Human % @ mRNA level
            'illumina_qc': r.get('illumina_qc'),
            'dup_funnel': r.get('dup_funnel'),
            'fullmol_funnel': r.get('fullmol_funnel'),
            **({'read_composition': r['read_composition']} if 'read_composition' in r else {}),
            'length_p': ss.get('length_p'), 'length_rho': ss.get('length_rho'),
            'g4': _g_pct_at_dn(bc, 3), 'g5': _g_pct_at_dn(bc, 4),
            'n_dup': n_dup, 'dup_pct': (n_dup / total * 100) if total else 0,
            'n_dup_mrna': n_dup_mrna, 'dup_pct_mrna': (n_dup_mrna / n_mrna_fp * 100) if n_mrna_fp else 0,
            'mrna_pct': _catpct('mRNA (protein-coding)'),
            'mrna_hu_n': (tx.get('mRNA (protein-coding)', {}) or {}).get('human', 0) or 0,
            'mrna_ms_n': (tx.get('mRNA (protein-coding)', {}) or {}).get('mouse', 0) or 0,
            'mrna_sp_tot': sum(v for v in (tx.get('mRNA (protein-coding)', {}) or {}).values()
                               if isinstance(v, (int, float))),
            'rrna_pct': _catpct('rRNA'),
            'mt_pct': _catpct('Mitochondrial'),
            'canon_med': canon.get('median'),
            'aln_median': ss.get('aln_median'), 'aln_mean': ss.get('aln_mean'),
            'tso_ge2_pct': ss.get('tso_ge2_pct'),
            'struct_pct': ss.get('struct_pct') or {},
            'hu_aln': ss.get('hu_aln'), 'ms_aln': ss.get('ms_aln'),
            'hu_mean': ss.get('hu_mean'), 'ms_mean': ss.get('ms_mean'),
            'err_pct': ss.get('err_pct'),
            'pct_polyt': _rate('reads_with_polyt'), 'pct_bobc': _rate('reads_with_bob_c'),
            'pct_tso': _rate('reads_with_tso'), 'pct_2xtso': _rate('reads_multi_tso'),
            'empty_pct_bc': (r.get('empty_products') or {}).get('pct_bc'),
            'empty_pct_all': (r.get('empty_products') or {}).get('pct_all'),
            'polydtless_pct_all': (r.get('polydt_less') or {}).get('pct_all'),
        }
    bc_summaries = {bc: _ms_vals(bc) for bc in bc_list}
    # Raw-read accounting rows read the prep sidecar through this dict.
    for _bc in bc_list:
        _p = (all_results.get(_bc) or {}).get('prep_stats')
        if _p and isinstance(bc_summaries.get(_bc), dict):
            bc_summaries[_bc]['prep_stats'] = _p

    def _f1(v, sfx=''):
        return '-' if v is None else f'{v:.2f}{sfx}'
    def _fi(v):
        return '-' if v is None else f'{int(round(v))}'
    def _fi_mm(mean, median):
        """'mean / median' insert length; falls back to whichever is present."""
        if mean is None and median is None:
            return '-'
        if mean is None or median is None:
            return _fi(mean if mean is not None else median)
        return f'{int(round(mean))} / {int(round(median))}'
    def _pct_int(v):
        return '-' if v is None else f'{v * 100:.1f}%'
    def _fmt_length_p(s):
        p = s.get('length_p'); rho = s.get('length_rho')
        if p is None:
            return '-'
        return f'p={p:.2g} (ρ={rho:+.2f})' if rho is not None else f'p={p:.2g}'

    # Chemistry-specific structural rows (the STAR/correctness rows are shared).
    _is_tso = (EXPERIMENT_CONFIG.get('barcode_location') == 'TSO')
    if _is_tso:
        def _sp(s, k):
            return (s.get('struct_pct') or {}).get(k)
        _g_rows = []   # the GGG template-switch is shown in the element rows below
        _canon_row = ('Insert mean / median (barcoded scorable reads, bp)',
                      lambda s: _fi_mm(s.get('aln_mean'), s.get('aln_median')),
                      'minmax', lambda s: s.get('aln_median'))
        _elem_rows = [
            ('% reads with the full R1 end (22 nt, barcode or not)', lambda s: _f1(_sp(s, 'backbone'), '%'), 'minmax', lambda s: _sp(s, 'backbone')),
            ('% reads with TSO barcode', lambda s: _f1(_sp(s, 'barcode'), '%'), 'minmax', lambda s: _sp(s, 'barcode')),
            ('% reads with ≥2 TSO (chimera)', lambda s: _f1(s.get('tso_ge2_pct'), '%'), 'minmax_inv', lambda s: s.get('tso_ge2_pct')),
        ]
    else:
        _g_rows = [
            ('G4 % (non-templated G at TSO +4)', lambda s: _f1(s['g4'], '%'), 'pct33_100', lambda s: s['g4']),
            ('G5 % (non-templated G at TSO +5)', lambda s: _f1(s['g5'], '%'), 'pct33_100', lambda s: s['g5']),
        ]
        _canon_row = ('Avg insert (TSO-insert-BobC, bp)', lambda s: _fi(s['canon_med']),
                      'minmax', lambda s: s['canon_med'])
        _elem_rows = [
            ('% reads with BobC', lambda s: _f1(s['pct_bobc'], '%'), 'minmax', lambda s: s['pct_bobc']),
            ('% reads with TSO', lambda s: _f1(s['pct_tso'], '%'), 'minmax', lambda s: s['pct_tso']),
            ('% reads with 2x TSO', lambda s: _f1(s['pct_2xtso'], '%'), 'minmax', lambda s: s['pct_2xtso']),
        ]
    # --- Barcode-accuracy filter funnel accessors (5 cumulative steps) --------
    # Accuracy after each swap-inclusion filter step; the Ms/Hu/terminal rows are
    # the filtered analogues (computed on the reads surviving all 4 filters).
    def _ff(s):
        return s.get('filter_funnel') or {}

    def _grun_fmt(s):
        ff = _ff(s); m = ff.get('grun_mean')
        if m is None:
            return '-'
        md = ff.get('grun_median')
        return f"{m:.2f} / {md:.2f}" if md is not None else f"{m:.2f}"

    def _grun0_fmt(s):
        # % of barcode-callable reads with a 0-length G-run = NO template switch.
        # Denominator = all reads in grun_hist (incl. the 0 bin). Keys may be str/int.
        hist = _ff(s).get('grun_hist') or {}
        tot = sum(hist.values())
        if not tot:
            return '-'
        z = hist.get('0', hist.get(0, 0))
        return f"{100 * z / tot:.1f}%"

    def _ff_step_acc(s, i):
        steps = _ff(s).get('steps') or []
        return steps[i].get('acc') if i < len(steps) else None

    def _ff_step_n(s, i):
        steps = _ff(s).get('steps') or []
        return steps[i].get('n') if i < len(steps) else None

    def _ff_fmt_step(s, i):
        a = _ff_step_acc(s, i)
        if a is None:
            return '-'
        n = _ff_step_n(s, i)
        return f"{a:.2f}% ({n:,})" if n is not None else f"{a:.2f}%"

    def _ff_val_step(s, i):
        a = _ff_step_acc(s, i)
        return (a / 100.0) if a is not None else None

    def _ff_fmt(s, key):
        v = _ff(s).get(key)
        return '-' if v is None else f"{v:.2f}%"

    def _ff_val(s, key):
        v = _ff(s).get(key)
        return (v / 100.0) if v is not None else None

    # Fullmol accuracy-funnel accessors for the Master Summary barcode-accuracy rows
    # (insert>30 -> +not rRNA -> +not MT -> +exact 1x barcode, then a mouse/human
    # split). Each cell is "acc% (n = reads scored at that level)".
    def _fm(s):
        fm = s.get('fullmol_funnel')
        if fm:
            return fm
        # Non-24plex fallback: fullmol_funnel is 24plex-only, so for tso-bob /
        # old-TSO runs expose the headline accuracy from the (universal) filter
        # funnel — exact-barcode row <- terminal_correct, Mouse/Human % correct
        # <- ms/hu_correct, Reads(total) <- total_reads. The intermediate funnel
        # steps (insert / +not rRNA / +not MT) are 24plex-only and stay blank.
        sfs = s.get('star_funnel_spec') or {}
        ff = s.get('filter_funnel') or {}
        if not sfs and not ff:
            return {}
        nf = ff.get('n_final', 0) or 0
        tc = ff.get('terminal_correct')
        return {
            'n_total': s.get('total_reads'),
            # Mouse/Human % correct at the mRNA level (canonical funnel), consistent with 24plex.
            'ms_correct': sfs.get('_ms_correct'), 'ms_n': sfs.get('_ms_n', 0),
            'hu_correct': sfs.get('_hu_correct'), 'hu_n': sfs.get('_hu_n', 0),
            'steps': ([{'key': 'barcode', 'acc': tc, 'acc_n': nf}] if tc is not None else []),
        }

    def _fm_total(s):
        # Reads EXAMINED. Under dedup-first, _fm(s)['n_total'] counts only the reads that
        # SURVIVED deduplication, which would silently understate the sequencing depth by
        # ~20x. The pre-filter count is captured before any filtering; prefer it.
        t = _df(s).get('fastq_reads') or _fm(s).get('n_total')
        return f"{t:,}" if t else '-'

    def _fm_step(s, key):
        return next((st for st in (_fm(s).get('steps') or []) if st.get('key') == key), None)

    def _fm_fmt(s, key):
        st = _fm_step(s, key)
        if not st or st.get('acc') is None:
            return '-'
        return f"{st['acc']:.2f}% (n={st.get('acc_n', 0):,})"

    def _fm_val(s, key):
        st = _fm_step(s, key)
        a = st.get('acc') if st else None
        return (a / 100.0) if a is not None else None

    # 0-100 accuracy (for the fixed-threshold color functions, which take percentages).
    def _fm_pct(s, key):
        st = _fm_step(s, key)
        return st.get('acc') if st else None

    def _fm_sp_fmt(s, pfx):
        a = _fm(s).get(f'{pfx}_correct')
        return '-' if a is None else f"{a:.2f}% (n={_fm(s).get(f'{pfx}_n', 0):,})"

    def _fm_sp_val(s, pfx):
        a = _fm(s).get(f'{pfx}_correct')
        return (a / 100.0) if a is not None else None

    def _fm_sp_pct(s, pfx):
        return _fm(s).get(f'{pfx}_correct')

    # % of reads whose literal 3' terminus IS the ONT adapter (read_composition,
    # deferred). Reads "pending" until the deferred pass has run.
    def _rc_ont_end_pct(s):
        r = s.get('read_composition') or {}
        n = r.get('n') or 0
        # Same rule as ONT@end in the read-composition provenance table (end_boundary: the
        # native-barcode flank or an ONT-adapter seed within the last 110 nt), so the row and
        # the table share one definition rather than classifying the literal terminus.
        v = (r.get('n_ont_end') or 0)
        return (100.0 * v / n) if n else None

    def _rc_ont_end_fmt(s):
        if 'read_composition' not in s:
            return _PENDING
        p = _rc_ont_end_pct(s)
        return '-' if p is None else f"{p:.1f}%"

    # Illumina QC / dup funnel are DEFERRED (computed after the first PDF, then the
    # PDF is regenerated). Until then their rows read "not run yet" (if the chemistry
    # is applicable) or "n/a" (if it isn't).
    _iq_app = (any(k in (EXPERIMENT_CONFIG.get('construct_description') or '').upper()
                   for k in ('P5', 'P7', 'ILLUMINA'))
               or any(_r.get('illumina_qc') is not None for _r in all_results.values()))
    _ud_app = 'UMI' in (EXPERIMENT_CONFIG.get('construct_description') or '').upper()
    _PENDING = 'not run yet'

    # --- Threshold color palettes for Master Summary rows -------------------
    # These bypass the min->max gradient: a raw value maps directly to one of a
    # few fixed colors by fixed cutoffs. LY=light-yellow, LG=light-green,
    # DG=deep-green, MG=magenta, RD=red.
    _C_LY, _C_LG, _C_DG, _C_MG, _C_RD = '#FFF2A8', '#C6EFCE', '#63BE7B', '#F08CC0', '#E06666'

    def _col_barcode(v):   # accuracy %: <99 yellow, 99-99.9 light green, >99.9 deep green
        if v is None: return None
        if v < 99.0:  return _C_LY
        if v <= 99.9: return _C_LG
        return _C_DG

    def _col_ont_end(v):   # >75 green, 50-75 yellow, <50 magenta
        if v is None: return None
        if v > 75:    return _C_DG
        if v >= 50:   return _C_LY
        return _C_MG

    def _col_mrna(v):      # >75 green, 50-75 light green, 25-50 yellow, <25 magenta
        if v is None: return None
        if v > 75:    return _C_DG
        if v > 50:    return _C_LG
        if v >= 25:   return _C_LY
        return _C_MG

    def _col_rrna(v):      # <5 green, 5-15 light green, 15-30 yellow, >30 red
        if v is None: return None
        if v < 5:     return _C_DG
        if v < 15:    return _C_LG
        if v <= 30:   return _C_LY
        return _C_RD

    def _col_mt(v):        # <2 green, 2-5 light green, 5-10 yellow, >10 magenta
        if v is None: return None
        if v < 2:     return _C_DG
        if v < 5:     return _C_LG
        if v <= 10:   return _C_LY
        return _C_MG

    def _col_artifact(v):  # white (0%) -> magenta (>=20%); high = more PCR artifact
        if v is None: return None
        x = max(0.0, min(1.0, v / 20.0))
        r = int(round(0xFF + (0xF0 - 0xFF) * x))
        g = int(round(0xFF + (0x8C - 0xFF) * x))
        b = int(round(0xFF + (0xC0 - 0xFF) * x))
        return '#%02X%02X%02X' % (r, g, b)

    def _col_end_incorrect(v):
        # Cross-orientation "incorrect %" rows (lower = better). v is a FRACTION
        # (_iq_val returns pct/100). Continuous gradient: green ≤1% → light green
        # 2.5% → yellow 10% → magenta ≥15% (so anything >10% trends magenta),
        # piecewise-linear between the anchor stops. Palette matches the table's
        # discrete _C_* colors so it blends with the rest of the summary.
        if v is None:
            return None
        _p = v * 100.0
        _stops = [(1.0,  (0x63, 0xBE, 0x7B)),   # green       (_C_DG)
                  (2.5,  (0xC6, 0xEF, 0xCE)),   # light green (_C_LG)
                  (10.0, (0xFF, 0xF2, 0xA8)),   # yellow      (_C_LY)
                  (15.0, (0xF0, 0x8C, 0xC0))]   # magenta     (_C_MG)
        if _p <= _stops[0][0]:
            _rgb = _stops[0][1]
        elif _p >= _stops[-1][0]:
            _rgb = _stops[-1][1]
        else:
            _rgb = _stops[0][1]
            for _i in range(len(_stops) - 1):
                _lo, _ca = _stops[_i]; _hi, _cb = _stops[_i + 1]
                if _lo <= _p <= _hi:
                    _t = (_p - _lo) / (_hi - _lo) if _hi > _lo else 0.0
                    _rgb = tuple(int(round(_ca[_j] + (_cb[_j] - _ca[_j]) * _t)) for _j in range(3))
                    break
        return '#%02X%02X%02X' % _rgb

    def _iq(s):
        return s.get('illumina_qc') or {}
    def _iq_pct(s, k):
        if s.get('illumina_qc') is None:
            return _PENDING if _iq_app else 'n/a'
        v = _iq(s).get(k)
        return '-' if v is None else f"{v:.1f}%"
    def _iq_val(s, k):
        v = _iq(s).get(k)
        return (v / 100.0) if v is not None else None
    def _iq_idx(s, which):
        if s.get('illumina_qc') is None:
            return _PENDING if _iq_app else 'n/a'
        d = _iq(s).get(f'{which}_dominant')
        m = _iq(s).get(f'{which}_match_pct')
        return '-' if m is None else f"{m:.1f}% [{d}]"

    def _iq_neb(s):
        """Identify the i5/i7 primer by matching the detected index window to the
        saved NEBNext index list — 'i501/i703' etc., or 'unknown' (non-NEB/TruSeq)."""
        if s.get('illumina_qc') is None:
            return _PENDING if _iq_app else 'n/a'
        import illumina_adapter_qc_illumina as _il
        i5 = _il.match_neb_index(_iq(s).get('i5_dominant'), 'i5') or 'unknown'
        i7 = _il.match_neb_index(_iq(s).get('i7_dominant'), 'i7') or 'unknown'
        return f"i5={i5} / i7={i7}"

    # Confident PCR-duplicate FUNNEL (gene -> +UMI -> +~position -> +strand) over
    # protein-coding non-ribo mRNA reads with a callable 6nt TSO UMI. DEFERRED like
    # the UMI/Illumina blocks: reads "not run yet" until the PDF is regenerated.
    def _df(s):
        return s.get('dup_funnel') or {}

    _sat_cache = {}

    def _sat(s):
        """Good-Turing saturation for this barcode, from the duplicity histogram in the
        dup sidecar. Cached: the summary is recomputed per row otherwise.
        Returns {} for runs whose sidecar predates the histogram, so those rows read '-'
        rather than inventing a number."""
        # SCOPE: prefer the aligned-scope histogram so these rows share a denominator
        # with the duplicate rate (which is aligned-only, the Picard / CellRanger
        # convention). Runs made before that key existed fall back to the FASTQ-scope
        # histogram, which counts unmapped molecules too and so reads a little lower.
        _d = _df(s) or {}
        h = _d.get('duplicity_hist') or _d.get('duplicity_hist_aligned')   # definition D first
        if not h:
            return {}
        key = id(s)
        if key not in _sat_cache:
            try:
                import umi_dedup_illumina as _ud_sat
                _ud_sat.configure_umi(EXPERIMENT_CONFIG.get('rt_primers_used'))
                _sat_cache[key] = _ud_sat.saturation_summary(h) or {}
            except Exception as _e:                       # noqa: BLE001
                print(f"  WARNING: saturation summary failed ({_e!r})")
                _sat_cache[key] = {}
        return _sat_cache[key]
    def _df_umi(s):
        if s.get('dup_funnel') is None:
            return _PENDING if _ud_app else 'n/a'
        d = _df(s)
        if not d or d.get('n_umi') is None:
            return '-'
        return f"{d['n_umi']:,} / {d['n_mrna']:,} ({d['umi_pct']:.1f}%)"
    def _df_fmt(s, key):
        if s.get('dup_funnel') is None:
            return _PENDING if _ud_app else 'n/a'
        d = _df(s)
        if not d or d.get(key) is None:
            return '-'
        return f"{d[key]:,} ({d[key + '_pct']:.2f}%)"
    def _df_val(s, key):
        v = _df(s).get(key)
        return (v / 100.0) if v is not None else None

    def _df_count_fmt(s, key):
        """Plain read COUNT from the dup funnel (no %). Same deferred-value handling as
        _df_fmt: the dup funnel is computed after the first PDF render."""
        if s.get('dup_funnel') is None:
            return _PENDING if _ud_app else 'n/a'
        v = _df(s).get(key)
        return f"{v:,}" if v is not None else '-'

    # Illumina P5/P7 adapter QC — the headline metrics kept in the Master Summary
    # (below a thick divider) via `+ _illumina_rows`: the 3 color-coded incorrect-%
    # rows + the i5/i7 index-match rows. The verbose diagnostics (opposite-ends %,
    # NEBNext index identity, P5/P7→Read contiguity) were MOVED to their own table
    # on the "incorrect symmetric ends" plot page (see _il_moved above). DEFERRED:
    # reads "not run yet" in the first PDF, populated when the PDF is regenerated.
    # ILLUMINA: these are REMOVED, not greyed. Counter-intuitively
    # they get LESS meaningful on real Illumina data, not more: every metric here
    # (_opp_ends / _both_ends_same) requires BOTH molecule termini inside ONE read,
    # which only an ONT read spanning the whole molecule ever has — a 150 bp mate sees
    # one end by construction. Measured on a 2x150 library: partial_p5 0.5%, partial_p7
    # 0.1%, p5p7_opposite_ends 0.0%, i7_n = 1. Index hopping is also unmeasurable
    # because the FASTQs arrive already demultiplexed, so the index pair is constant.
    # Reporting 0.0% here would read as "no adapter artifacts" rather than "not measurable".
    _illumina_rows = [] if IS_ILLUMINA else [
        ('P5-P5 or P7-P7 ends (incorrect %)',          lambda s: _iq_pct(s, 'p5p5_or_p7p7_pct'), ('thresh', _col_end_incorrect), lambda s: _iq_val(s, 'p5p5_or_p7p7_pct')),
        ('i5-i5 or i7-i7 index seqs (incorrect %)',      lambda s: _iq_pct(s, 'idx_i5i5_or_i7i7_pct'), ('thresh', _col_end_incorrect), lambda s: _iq_val(s, 'idx_i5i5_or_i7i7_pct')),
        ('TruSeq R1-R1 or R2-R2 ends (incorrect %)',     lambda s: _iq_pct(s, 'truseq_r1r1_or_r2r2_pct'), ('thresh', _col_end_incorrect), lambda s: _iq_val(s, 'truseq_r1r1_or_r2r2_pct')),
    ]

    # mRNA-species split (of protein-coding mRNA reads; columns sum to ~100%), shaded
    # white 0% -> yellow 50% -> green 100% like the standalone mRNA-species table.
    def _col_mrna_split(v):
        if v is None:
            return None
        c = _wyg_0_100(v)   # colors.Color -> hex string (make_heatmap_table wraps in HexColor)
        return '#%02X%02X%02X' % (round(c.red * 255), round(c.green * 255), round(c.blue * 255))
    def _mrna_split_fmt(s, key):
        _n = s.get(key, 0) or 0; _t = s.get('mrna_sp_tot', 0) or 0
        return f"{_n:,} ({ipct(_n, _t)})" if _t else f"{_n:,}"
    def _mrna_split_val(s, key):
        _n = s.get(key, 0) or 0; _t = s.get('mrna_sp_tot', 0) or 0
        return (100.0 * _n / _t) if _t else None

    # Chemistry-adaptive per-step barcode-accuracy rows: the first up-to-4 unified
    # funnel steps (24plex: insert/notrRNA/notMT/exact-1×; tso-bob: has-7mer/GGG/
    # insert/C28). Labels + values come from each barcode's unified_funnel, so the
    # Master-Summary "% correct by step" auto-switches with the chemistry.
    def _ufacc_fmt(s, key):
        _st = next((x for x in ((s.get('unified_funnel') or {}).get('steps') or []) if x['key'] == key), None)
        if not _st or _st.get('acc') is None:
            return '-'
        return f"{_st['acc']:.2f}% (n={_st.get('acc_n', 0):,})"
    def _ufacc_pct(s, key):
        _st = next((x for x in ((s.get('unified_funnel') or {}).get('steps') or []) if x['key'] == key), None)
        return _st.get('acc') if _st else None
    # First 4 PRESENT canonical steps (has-7mer/exact/mRNA/barcode-insert-distance — the
    # barcode-accuracy core; tso-bob lacks the distance step so it shows its own 4th step).
    # ONE headline accuracy row: the deepest barcode-core step (has-barcode / mRNA / distance)
    # this chemistry has. The per-level view lives in section 1; repeating the per-level
    # rows here would duplicate that table verbatim.
    _core = ([s for s in _uf(bc_list[0])['steps'] if s.get('present') and s['key'] in ('has7', 'mrna', 'dist')]
             if bc_list else [])
    _acc_step = _core[-1] if _core else None
    _funnel_rows = ([(("Barcode % correct - mRNA + barcode-insert distance level (all levels in section 1)" if _acc_step['key'] == 'dist'
                       else "Barcode % correct - mRNA level (all levels in section 1)"),
                      (lambda s, _k=_acc_step['key']: _ufacc_fmt(s, _k)),
                      ('thresh', _col_barcode),
                      (lambda s, _k=_acc_step['key']: _ufacc_pct(s, _k)))] if _acc_step else [])
    _last_funnel_label = _funnel_rows[0][0] if _funnel_rows else None

    # Mouse/Human % correct are distance-filtered only for 24plex (fullmol funnel); for
    # tso-bob the values fall back to the mRNA level, so the label must NOT claim distance.
    _hums_lbl = ('protein-coding incl. RP + distance-filtered'
                 if _ba.is_24plex(EXPERIMENT_CONFIG.get('tso_species_map'), EXPERIMENT_CONFIG.get('tso_arch'))
                 else 'protein-coding incl. RP')
    def _ps(s):
        return s.get('prep_stats') or {}
    def _ps_n(s):
        n = _ps(s).get('n'); return f"{n:,}" if n else '-'
    def _ps_frac(s, key):
        p = _ps(s); n = p.get('n') or 0; v = p.get(key)
        return '-' if not n or v is None else f"{v:,} ({100.0 * v / n:.1f}%)"
    def _ps_pct(s, key):
        p = _ps(s); n = p.get('n') or 0; v = p.get(key)
        return (100.0 * v / n) if n and v is not None else None
    def _mapped_n(s):
        _tx = s.get('star_taxonomy') or {}
        return sum(sum(v.values()) for v in _tx.values() if isinstance(v, dict))
    def _mapped_pct(s):
        _d = ((s.get('dup_funnel') or {}).get('dedup_kept') or s.get('total_reads') or 0)
        return (100.0 * _mapped_n(s) / _d) if _d else None
    def _mapped_fmt(s):
        _d = ((s.get('dup_funnel') or {}).get('dedup_kept') or s.get('total_reads') or 0); _p = _mapped_pct(s)
        return '-' if _p is None else f"{_p:.1f}% ({_mapped_n(s):,} of {_d:,})"
    summary_rows = [
        # Raw-read accounting FIRST: 'Reads (total)' below is what survived
        # prep (bobcode called, poly(A) read-through trimmed, insert >= 30 nt), and a
        # library can lose most of its pairs there without any of the rows below saying
        # so. Some libraries lose 80-89% of pairs to reads that run from the bobcode
        # straight into poly(A) (no insert). Present only when the prep sidecar is staged.
        ('Read pairs sequenced (raw)',                    _ps_n, 'none', lambda s: _ps(s).get('n')),
        ('Usable read pairs after prep (bobcode + insert >= 30 nt)',
         lambda s: _ps_frac(s, 'written'), 'none', lambda s: _ps_pct(s, 'written')),
        ('Dropped: insert < 30 nt after poly(A) trim (bobcode into poly(A))',
         lambda s: _ps_frac(s, 'too_short_after_trim'), 'none', lambda s: _ps_pct(s, 'too_short_after_trim')),
        ('Dropped: no bobcode called',
         lambda s: _ps_frac(s, 'no_bobcode'), 'none', lambda s: _ps_pct(s, 'no_bobcode')),
        # Duplicate rate: it frames every count below it. Its denominator is the
        # PRE-deduplication aligned read set, so it keeps describing the real library
        # even though every other row is now computed on deduplicated reads.
        # ONE duplicate definition in the report: definition D (same bobcode, contig,
        # strand, 5' position, UMI within the chemistry's edit distance), over ALIGNED
        # reads with a UMI call. A position-independent UMI cluster (`dup_pct`/`unique`)
        # is a different definition from the keep-set every other row is computed on,
        # and the two cannot be reconciled from the page; that value stays in the
        # sidecar as `dup_pct` (UMI-only) but is not shown.
        ('Duplicate rate (definition D: same position + fragment length + UMI, aligned reads)' if DEDUP_KEY_EXTRA
         else 'Duplicate rate (definition D: same position + UMI, aligned reads)',
         lambda s: (lambda d: '-' if not d or d.get('dup_pct_defD_aligned') is None else
                    f"{d['dup_pct_defD_aligned']:.1f}%  ({d['dedup_kept_aligned']:,} molecules / "
                    f"{d['mapped_with_umi']:,} aligned reads)")(_df(s)),
         'minmax_inv', lambda s: (_df(s) or {}).get('dup_pct_defD_aligned')),
        # No 'Reads (total, after prep)' row: it is the same count as 'Usable read pairs after
        # prep' (that row also carries the %). The funnel table still prints Total reads.
        # Count only, no accuracy — deduplication is not a structural filter, so it sits
        # outside the barcode-accuracy funnel rows that follow. This is the number of
        # reads EVERY row below is computed on (the dedup filter's output), which is
        # slightly larger than the 'unique' figure in the duplicate-rate row above:
        # that one counts only reads aligned to HUMAN_/MOUSE_, this one is built from the
        # FASTQ so unmapped reads are represented too.
        ('Reads kept after deduplication (all prepped reads, unmapped included)',
         lambda s: _df_count_fmt(s, 'dedup_kept'),
         'none', lambda s: _df(s).get('dedup_kept')),
        ('% reads mapped (STAR, of reads after deduplication)', _mapped_fmt, 'minmax', _mapped_pct),
    ] + _funnel_rows + [
        ('Mouse % correct (aligned, ' + _hums_lbl + ')',  lambda s: _fm_sp_fmt(s, 'ms'), ('thresh', _col_barcode), lambda s: _fm_sp_pct(s, 'ms')),
        ('Human % correct (aligned, ' + _hums_lbl + ')',  lambda s: _fm_sp_fmt(s, 'hu'), ('thresh', _col_barcode), lambda s: _fm_sp_pct(s, 'hu')),
        ('% mRNA human (of mRNA, STAR)',               lambda s: _mrna_split_fmt(s, 'mrna_hu_n'), ('thresh', _col_mrna_split), lambda s: _mrna_split_val(s, 'mrna_hu_n')),
        ('% mRNA mouse (of mRNA, STAR)',               lambda s: _mrna_split_fmt(s, 'mrna_ms_n'), ('thresh', _col_mrna_split), lambda s: _mrna_split_val(s, 'mrna_ms_n')),
        ('mRNA % (protein-coding, RP excluded, of mapped)', lambda s: _f1(s['mrna_pct'], '%'), ('thresh', _col_mrna), lambda s: s['mrna_pct']),
        ('rRNA % (of mapped)', lambda s: _f1(s['rrna_pct'], '%'), ('thresh', _col_rrna), lambda s: s['rrna_pct']),
        ('Mitochondrial % (of mapped)', lambda s: _f1(s['mt_pct'], '%'), ('thresh', _col_mt), lambda s: s['mt_pct']),
        ("% reads ending in ONT adapter (block in the last 110 nt)", _rc_ont_end_fmt, ('thresh', _col_ont_end), _rc_ont_end_pct),
        ('G-run mean / median (non-templated G after TSO; mapped barcoded reads)', _grun_fmt, 'none', lambda s: _ff(s).get('grun_mean')),
        ('No template switch (G-run = 0; mapped barcoded reads)', _grun0_fmt, 'none', lambda s: None),
    ] + _g_rows + [
        ('Empty products % (of all reads)', lambda s: _f1(s.get('empty_pct_all'), '%'),
         ('thresh', _col_artifact), lambda s: s.get('empty_pct_all')),
        ('Poly(dT)-less artifacts % (of all reads)', lambda s: _f1(s.get('polydtless_pct_all'), '%'),
         ('thresh', _col_artifact), lambda s: s.get('polydtless_pct_all')),
        # No Lander-Waterman "Estimated library complexity" row: LW assumes uniform
        # sampling, which RNA-seq violates. Measured on one library it called the pool
        # exhausted (1.17M estimated vs 1.17M observed) while a THIRD of its molecules
        # had been seen exactly once. Across assays -- which is the point of a benchmark
        # -- comparing two LW numbers compares their expression skew as much as their
        # libraries. The sidecar still records lib_size_M; the exact saturation rows
        # below and the saturation curve carry the information without the model.
        # GOOD-TURING saturation, computed from the duplicity histogram the dedup step
        # records. These are EXACT and assume nothing about abundance:
        #   coverage       = 1 - f1/N   fraction of the molecule POOL already sampled
        #   marginal yield = f1/N       new molecules per additional read
        # NOTE coverage is a SATURATION measure, not a quality one: a poor library with
        # a tiny pool samples out quickly and scores HIGH (a negative control can read
        # 98.5% on 164k molecules). Read it beside the molecule count, never alone.
        # Coverage and marginal yield are not shown as rows: both are f1/N, the same fact
        # as the singleton row below; the saturation curve carries the trend.
        ('Molecules seen exactly once (singletons)',
         lambda s: (lambda d: '-' if not d or d.get('singletons') is None else
                    f"{d['singletons']:,} ({d['singleton_pct']:.1f}%)")(_sat(s)),
         'minmax', lambda s: (_sat(s) or {}).get('singleton_pct')),
        ('Reads to reach half the molecules observed',
         lambda s: (lambda d: '-' if not d or d.get('reads_to_half_of_observed') is None
                    else f"{d['reads_to_half_of_observed']:,}")(_sat(s)),
         'minmax', lambda s: (_sat(s) or {}).get('reads_to_half_of_observed')),
        # Duplicates whose copies are physical NEIGHBOURS on the flowcell (same tile,
        # within OPTICAL_PIXEL_DIST px) = optical + NovaSeq ExAmp pad-hopping. A sequencer
        # artifact, not library PCR duplication — subtract it from the row below.
        # Only measurable on CONTIGUOUS reads; a strided subsample reports ~0 regardless.
        ('Flowcell duplicates (optical/ExAmp) %',
         lambda s: (lambda d: '-' if not d or d.get('dup_flowcell_pct') is None else
                    f"{d['dup_flowcell_pct']:.2f}%  ({d['dup_flowcell']:,} reads)")(_df(s)),
         'minmax_inv', lambda s: (_df(s) or {}).get('dup_flowcell_pct')),
        _canon_row,
        ('ONT error rate (aligned, %)', lambda s: _f1(s.get('err_pct'), '%'), 'minmax_inv', lambda s: s.get('err_pct')),
    ] + _elem_rows + _illumina_rows
    # Headline band (answers 'which one is the final value'): the rows the
    # lab reads first come first, in bold, above a thick line; everything else is detail below.
    _HEAD = ('Reads (total)', 'Read pairs sequenced', 'Usable read pairs', 'Duplicate rate', '% reads mapped',
             'Barcode % correct', '% mRNA human', '% mRNA mouse', 'mRNA % (', 'rRNA %', 'Mitochondrial %',
             '% reads with TSO barcode', '% reads with ≥2 TSO')
    def _head_rank(r):
        return next((i for i, pfx in enumerate(_HEAD) if isinstance(r[0], str) and r[0].startswith(pfx)), None)
    _head = sorted((r for r in summary_rows if _head_rank(r) is not None), key=_head_rank)
    summary_rows = _head + [r for r in summary_rows if _head_rank(r) is None]
    _n_head = len(_head)

    # ---------------- ILLUMINA row policy (applied as a post-pass) -------------
    # Done here, not inside the literals above, so the ONT row definitions stay
    # untouched and every Illumina decision is in ONE greppable place.
    #
    #  DROP  (#11) - insert-length. On single-end R2 the aligned length is a READ
    #        length capped at 150, not the molecule/fragment length the ONT rows
    #        reported. Keeping it would silently redefine the metric.
    #  GREY  (#9)  - still listed (label strings preserved for the cross-run join)
    #        but rendered 'n/a (Illumina)' in grey with no threshold shading, because
    #        the underlying property does not exist on an Illumina read.
    if IS_ILLUMINA:
        _DROP_ROWS = {
            # Insert length: on single-end R2 this is a READ length capped at 150, not
            # the molecule/fragment length the ONT rows reported.
            'Insert mean / median (barcoded scorable reads, bp)',
            'Insert mean / median (human-aligned scorable reads, bp)',
            'Insert mean / median (mouse-aligned scorable reads, bp)',
            # --- ONT-only properties that do not exist on an Illumina read ---------
            # ONT basecaller error; the Illumina mismatch rate is a different quantity
            # (0.35% here vs ~4-5% ONT) against ONT-calibrated thresholds.
            'ONT error rate (aligned, %)',
            "% reads ending in ONT adapter (block in the last 110 nt)",
            # TruSeq2 IS the R2 sequencing primer, so the backbone can never appear in
            # the read; 0% would read as "adapter degraded" rather than "unobservable".
            '% reads with the full R1 end (22 nt, barcode or not)',
            '% partial TSO backbone (barcode found, R1 end truncated)',
            # polyT is on R1; on R2 an A/T run means 3' polyA READ-THROUGH (a short
            # insert), so the row inverts its own meaning.
            '% reads with polyT',
            # Circular on Illumina: prep only emits reads that HAD a callable bobcode,
            # so any %-of-reads computed downstream is 100% by construction.
            '% reads with TSO barcode',
            # Chimera/multiplicity depends on finding two TSO units inside one read;
            # a 150 bp mate sees one end of the molecule, so this is not measurable.
            '% reads with ≥2 TSO (chimera)',
            # Both are ONT insert-less/artifact classes keyed on the full-molecule read.
            'Empty products % (of all reads)',
            'Poly(dT)-less artifacts % (of all reads)',
        }
        _GREY_ROWS = {
        }
        _kept = []
        for _r in summary_rows:
            _label = _r[0]
            if _label in _DROP_ROWS:
                continue
            if _label in _GREY_ROWS:
                _kept.append((_label, (lambda s: 'n/a (Illumina)'), 'grey', (lambda s: None)))
            else:
                _kept.append(_r)
        summary_rows = _kept

    summary_data = [['Metric'] + [bc_header[bc] for bc in bc_list]]
    summary_cell_values = {}
    summary_cell_colors = {}
    summary_grad_rows = set()
    summary_plain_rows = set()
    summary_grey_rows = set()   # ILLUMINA: rows kept but rendered 'n/a' in grey
    # The Dup-funnel block and every row below it get NO cell coloring (plain).
    _dup_start = next((i for i, r in enumerate(summary_rows, start=1)
                       if isinstance(r[0], str) and r[0].startswith(('Dup funnel', 'Duplicate rate'))), None)
    for ri_off, (label, formatter, mode, value_fn) in enumerate(summary_rows, start=1):
        raw_vals = [value_fn(bc_summaries[bc]) for bc in bc_list]
        normalized = [None] * len(bc_list)
        # Rows at/below the Dup-funnel block are plain by default — EXCEPT rows with
        # an explicit ('thresh', fn) gradient, which always color (e.g. the Illumina
        # cross-orientation "incorrect %" rows at the bottom of the table).
        _is_thresh = isinstance(mode, tuple) and mode[0] == 'thresh'
        # ILLUMINA 'grey' = keep the row, render it 'n/a' de-emphasised, never shade it.
        # NB: this must NOT `continue` — the row still has to be APPENDED to summary_data
        # below. Skipping the append would make greyed rows vanish from the PDF entirely AND
        # desync every row index beneath them (divider lines and threshold cell colors
        # are keyed on this same index), i.e. the opposite of "keep the layout".
        if mode == 'grey':
            summary_grey_rows.add(ri_off)
        plain = (mode == 'grey') or (
            (mode == 'none' or (_dup_start is not None and ri_off >= _dup_start)) and not _is_thresh)
        if plain:
            summary_plain_rows.add(ri_off)
        elif isinstance(mode, tuple) and mode[0] == 'thresh':
            # Fixed threshold coloring: raw value -> explicit hex via mode[1].
            summary_plain_rows.add(ri_off)  # skip default ±30% highlighting
            color_fn = mode[1]
            for ci, v in enumerate(raw_vals):
                c = color_fn(v)
                if c is not None:
                    summary_cell_colors[(ri_off, ci + 1)] = c
        elif mode == 'abs_spec':
            for ci, v in enumerate(raw_vals):
                normalized[ci] = v
        elif mode == 'pct33_100':
            for ci, v in enumerate(raw_vals):
                if v is not None:
                    normalized[ci] = max(0.0, min(1.0, (v - 33.0) / 67.0))
        elif mode in ('minmax', 'minmax_inv'):
            nums = [v for v in raw_vals if v is not None]
            if nums:
                vmin, vmax = min(nums), max(nums)
                if vmax > vmin:
                    for ci, v in enumerate(raw_vals):
                        if v is None:
                            continue
                        x = (v - vmin) / (vmax - vmin)
                        normalized[ci] = (1.0 - x) if mode == 'minmax_inv' else x
                else:
                    for ci, v in enumerate(raw_vals):
                        if v is not None:
                            normalized[ci] = 0.5
        row = [label]
        for ci, bc in enumerate(bc_list):
            row.append(formatter(bc_summaries[bc]))
            if normalized[ci] is not None:
                summary_cell_values[(ri_off, ci + 1)] = normalized[ci]
        summary_data.append(row)
        if not plain and not (isinstance(mode, tuple) and mode[0] == 'thresh'):
            summary_grad_rows.add(ri_off)

    # ---- machine-readable per-run summary dump (mirrors THIS Master Summary) ----
    # Writes <run>_master_summary.{json,tsv} next to the PDF so cross-run
    # aggregation reads exactly the numbers rendered above. This is the only place
    # the fully-derived values (unified_funnel accuracy, mRNA-filtered species %,
    # etc.) exist — they are NOT in the per-barcode *_speciesmix_results.json.
    # Best-effort: a dump error must never break the PDF build.
    try:
        import re as _re_ms
        def _ms_plain(x):
            if x is None:
                return ''
            _t = getattr(x, 'text', None)
            _s = _t if isinstance(_t, str) else str(x)
            return _re_ms.sub('<[^>]+>', '', _s).replace('&amp;', '&').replace('&nbsp;', ' ').strip()
        def _ms_num(v):
            if v is None or isinstance(v, bool):
                return None
            try:
                return round(float(v), 4)
            except (TypeError, ValueError):
                return None
        _ms_metrics = []
        for (_lbl, _fmtr, _mode, _vfn) in summary_rows:
            _cells = {}
            for _bc in bc_list:
                _s = bc_summaries[_bc]
                try:
                    _raw = _vfn(_s)
                except Exception:
                    _raw = None
                try:
                    _txt = _ms_plain(_fmtr(_s))
                except Exception:
                    _txt = ''
                _cells[_bc] = {'value': _ms_num(_raw), 'text': _txt}
            _ms_metrics.append({'metric': _ms_plain(_lbl), 'cells': _cells})
        _run_name = os.path.basename(os.path.normpath(output_dir))
        _ms_dump = {
            'run': _run_name,
            'chemistry': EXPERIMENT_CONFIG.get('chemistry') or EXPERIMENT_CONFIG.get('chem'),
            'experiment_title': EXPERIMENT_CONFIG.get('experiment_title'),
            'species_map': EXPERIMENT_CONFIG.get('primary_species_map'),
            'rt_primers_used': EXPERIMENT_CONFIG.get('rt_primers_used'),
            'condition_labels': EXPERIMENT_CONFIG.get('condition_labels'),
            'barcodes': list(bc_list),
            'barcode_header': {_bc: _ms_plain(bc_header.get(_bc, _bc)) for _bc in bc_list},
            'metrics': _ms_metrics,
        }
        _ms_json = os.path.join(output_dir, f"{_run_name}_master_summary.json")
        with open(_ms_json, 'w') as _mf:
            json.dump(_ms_dump, _mf, indent=2)
        _ms_tsv = os.path.join(output_dir, f"{_run_name}_master_summary.tsv")
        with open(_ms_tsv, 'w') as _mf:
            _mf.write('metric\t' + '\t'.join(bc_list) + '\n')
            for _m in _ms_metrics:
                _mf.write(_m['metric'] + '\t'
                          + '\t'.join((_m['cells'][_bc]['text'] or '') for _bc in bc_list) + '\n')
        print(f"Machine-readable summary: {os.path.basename(_ms_json)}")
    except Exception as _ms_e:
        print(f"WARN: master-summary dump failed: {_ms_e}")

    _tcat_bold = ParagraphStyle('TaxCatB', parent=_tcat_style, fontName='Helvetica-Bold')
    for _ri in range(1, len(summary_data)):
        if isinstance(summary_data[_ri][0], str):
            summary_data[_ri][0] = Paragraph(
                summary_data[_ri][0].replace('&', '&amp;'),
                _tcat_bold if _ri <= _n_head else _tcat_style)

    # Divider lines render as LINEBELOW row N. Keep the composition/structure group
    # lines, then bracket the Dup-funnel block with a thick line ABOVE (below the row
    # before it) and BELOW (below its last row), plus the thick line above the Illumina
    # block. Funnel row span is found by label so it survives row reordering.
    def _sri(prefix):
        return next((i for i, r in enumerate(summary_rows, start=1)
                     if isinstance(r[0], str) and r[0].startswith(prefix)), None)
    _grp = tuple(x for x in ((_sri(_last_funnel_label) if _last_funnel_label else None),  # #3: after last funnel step
                             _sri('Human % correct'),      # close the barcode-accuracy block
                             _sri('Mitochondrial %')) if x)  # close the mRNA/rRNA/MT composition group
    _fun = [i for i, r in enumerate(summary_rows, start=1)
            if isinstance(r[0], str) and r[0].startswith(('Dup funnel', 'Duplicate rate'))]
    _fun_div = (_fun[0] - 1, _fun[-1]) if _fun else ()
    _ill_div = ((len(summary_rows) - len(_illumina_rows),) if _illumina_rows else ())
    # thick line above and below the ONT error-rate row
    _err = _sri('ONT error rate')
    _err_div = ((_err - 1, _err) if _err else ())
    summary_table = make_heatmap_table(
        summary_data, col_widths=auto_col_widths(240 if n_bc <= 8 else 200, n_bc),
        font_size=7, specificity_gradient_rows=summary_grad_rows,
        cell_values=summary_cell_values, plain_rows=summary_plain_rows,
        cell_colors=summary_cell_colors, grey_rows=summary_grey_rows,
        divider_rows=tuple(_grp) + _fun_div + _ill_div + _err_div + ((_n_head,) if _n_head else ()))

    # Non-templated G-run distribution per BC (all BCs overlaid; the 0-run bin is
    # omitted and each curve renormalized over its reads with >=1 G, so shapes are
    # directly comparable). Data = filter_funnel.grun_hist (positional counts).
    def _fit_box(png, max_w, max_h):
        """{'width','height'} scaling `png` to fit max_w x max_h WITHOUT distortion.
        Reads the real pixel size, so changing a figure's panel count or figsize can
        never again silently squash it against a hardcoded aspect ratio."""
        try:
            from reportlab.lib.utils import ImageReader
            _w, _h = ImageReader(png).getSize()
            _s = min(max_w / _w, max_h / _h)
            return {'width': _w * _s, 'height': _h * _s}
        except Exception:
            return {'width': max_w, 'height': max_h}

    _grun_img = None
    try:
        _gh = {bc: {int(k): v for k, v in
                    ((all_results[bc].get('filter_funnel') or {}).get('grun_hist') or {}).items()}
               for bc in bc_list}
        if any(sum(v for k, v in _gh[bc].items() if k >= 1) for bc in bc_list):
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as _plt
            # ONE plot, LOG y, with the 0 bin included (0 = no template switch), rather
            # than a pair of panels (linear distribution + separate 0-run bar chart).
            # Log scale is the informative view here: it shows the 0 bin sitting ABOVE the
            # 1 and 2 bins, i.e. reads either fail to template-switch at all or add >=3 G —
            # a real dip that the linear axis flattens away.
            # TWO panels: barcode accuracy by G-run on top, the distribution below, on a
            # shared x-axis so a dip in accuracy can be read directly against how many
            # reads sit in that bin. Accuracy is mRNA-only and post-deduplication (see
            # star_taxonomy.star_grun_accuracy).
            _gacc = {bc: {int(k): v for k, v in
                          (all_results[bc].get('star_grun_accuracy') or {}).items()}
                     for bc in bc_list}
            _has_acc = any(sum(v[1] for v in _gacc[bc].values()) for bc in bc_list)
            if _has_acc:
                _fig, (_axa, _ax) = _plt.subplots(
                    2, 1, figsize=(9, 6.4), sharex=True,
                    gridspec_kw={'height_ratios': [1.15, 1], 'hspace': 0.10})
            else:
                _fig, _ax = _plt.subplots(figsize=(9, 3.6))
                _axa = None
            _cmap = _plt.get_cmap('tab10' if n_bc <= 10 else 'tab20')
            _XMAX = 10
            if _axa is not None:
                _MINN = 200      # bins thinner than this are noise, not signal
                for _i, bc in enumerate(bc_list):
                    _h = _gacc[bc]
                    _xs, _ys = [], []
                    for _k in range(0, _XMAX):
                        _c, _n = _h.get(_k, [0, 0])
                        if _n >= _MINN:
                            _xs.append(_k); _ys.append(100.0 * _c / _n)
                    _ct = sum(v[0] for k, v in _h.items() if k >= _XMAX)
                    _nt = sum(v[1] for k, v in _h.items() if k >= _XMAX)
                    if _nt >= _MINN:
                        _xs.append(_XMAX); _ys.append(100.0 * _ct / _nt)
                    if _xs:
                        _axa.plot(_xs, _ys, '-o', lw=1.4, ms=3,
                                  color=_cmap(_i % 20), label=bc_nick[bc])
                _axa.set_ylabel('barcode % correct\n(mRNA, deduplicated)')
                _axa.set_title('Barcode accuracy and template-switch G-run length',
                               fontsize=9, fontweight='bold')
                _axa.grid(alpha=0.25)
                for _s in ('top', 'right'):
                    _axa.spines[_s].set_visible(False)
            for _i, bc in enumerate(bc_list):
                h = _gh[bc]
                denom = sum(h.values())        # ALL barcode-callable reads, 0 bin included
                if not denom:
                    continue
                ys = [100 * h.get(k, 0) / denom for k in range(0, _XMAX)]
                ys.append(100 * sum(v for k, v in h.items() if k >= _XMAX) / denom)
                _ax.plot(range(0, _XMAX + 1), ys, '-o', lw=1.4, ms=3,
                         color=_cmap(_i % 20), label=bc_nick[bc])
            _ax.set_yscale('log')
            _ax.set_xticks(range(0, _XMAX + 1))
            _ax.set_xticklabels([str(k) for k in range(0, _XMAX)] + [f'{_XMAX}+'])
            _ax.set_xlabel('non-templated G-run length  (0 = no template switch)')
            _ax.set_ylabel('% of barcode-callable reads (log)')
            if _axa is None:
                _ax.set_title('Non-templated G-run distribution per BC (template switch)',
                              fontsize=9, fontweight='bold')
            _ax.grid(alpha=0.25, which='both')
            for _s in ('top', 'right'):
                _ax.spines[_s].set_visible(False)
            _ax.legend(fontsize=6, ncol=2, frameon=False)
            _plt.tight_layout()
            _grp_png = os.path.join(fig_dir, 'grun_distribution.png')
            _fig.savefig(_grp_png, dpi=150, bbox_inches='tight'); _plt.close(_fig)
            _grun_img = _grp_png
    except Exception as _e:
        print(f'  (G-run distribution plot skipped: {_e})')

    _ms_heading = Paragraph("Master Summary — All BCs at a Glance", heading_style)
    _methods_para = (
        Paragraph(
            "<i>Per-BC headline metrics. Species/swap rows are <b>pure-STAR</b> "
            "(combined-genome alignment); structure/adapter rows are "
            "from read-structure analysis. The top block uses <b>fixed-threshold</b> "
            "coloring (not a per-row min→max scale): <b>Barcode % correct</b> rows are "
            "yellow &lt;99%, light-green 99–99.9%, deep-green &gt;99.9%. <b>Barcode % "
            "correct</b> = of barcoded reads, the fraction whose STAR-aligned species "
            "matches the species barcode, reported at the deepest barcode-core funnel level (reads scored in "
            "parentheses); every level is in the Read Structure section (" + _ms_step_desc + " …). "
            "<b>Mouse/Human % correct</b> split the last barcode level "
            "by aligned species — of reads aligning to that species, the fraction "
            "correctly barcoded. <b>% reads ending in ONT adapter</b> = the fraction whose "
            "last 110 nt contain the ONT adapter block (green &gt;75%, yellow 50–75%, magenta "
            "&lt;50%). <b>mRNA %</b> (green &gt;75%, light-green 50–75%, yellow 25–50%, "
            "magenta &lt;25%), <b>rRNA %</b> (green &lt;5%, light-green 5–15%, yellow "
            "15–30%, red &gt;30%) and <b>Mitochondrial %</b> (green &lt;2%, light-green "
            "2–5%, yellow 5–10%, magenta &gt;10%) are of <b>mapped</b> reads (same "
            "denominator as the §3 composition table). Insert lengths are reported as "
            "<b>mean / median</b> (bp). "
            "<b>ONT error rate</b> = per-base error over STAR-aligned read cores "
            "(substitutions + indels ÷ aligned columns; soft-clipped ends "
            "excluded, so it is a lower bound). "
            "The <b>Dup-funnel</b> block and every row below it are shown uncolored. "
            "<b>Dup-funnel</b> rows detect PCR duplicates over protein-coding "
            "non-ribosomal <b>mRNA reads only</b> (with a callable 6-nt UMI): reads are "
            "collapsed into molecules by adding one identity criterion at a time — "
            "<b>same gene</b> (mRNA only), <b>+ UMI</b>, <b>+ start/end</b> within ±10 bp, "
            "<b>+ same strand</b> — and each row is the duplicate % surviving that step; "
            "the same-strand step isolates PCR-amplification duplicates from both-strand "
            "sequencing of one molecule.</i>", body_style))
    master_summary_elems = [
        _ms_heading,
        summary_table,
        Spacer(1, 12),
    ] + ([
        KeepTogether([
            Paragraph("<b>Non-templated G-run distribution</b> <i>— when the reverse transcriptase reaches the 5′ end of the RNA it "
                      "adds a few extra bases and switches to the TSO; that leaves a short run of G right after the barcode "
                      "(typically 3 to 5). Left: distribution of that G-run length per barcode, among reads with at least one G, "
                      "scaled so the shapes can be compared. Right: share of barcoded reads with no G-run at all, i.e. molecules "
                      "that carry a barcode without a template switch. This share is a lower bound: when a read has two possible "
                      "barcode matches, the one followed by a G-run is taken.</i>", body_style),
            Spacer(1, 4),
            # Size from the PNG's ACTUAL pixel aspect, fitted inside a box. A height
            # hardcoded to one aspect ratio (e.g. width*(3.4/13)) is correct only for a
            # single-panel figure and renders a taller two-panel version squashed.
            Image(_grun_img, **_fit_box(_grun_img, 9.6 * inch, 5.4 * inch)),
        ]),
        Spacer(1, 12),
    ] if _grun_img else []) + [
        PageBreak(),
    ]
    # The report OPENS with the Master Summary (heading + table): insert right
    # after the title and the how-to-read note (index 2). It ends with a PageBreak, so the detail sections
    # follow on the next page.
    for offset, e in enumerate(master_summary_elems):
        elems.insert(2 + offset, e)

    # APPENDIX at the very back: construct/provenance line, the auto-detected
    # bobcodes table, and the Master Summary metric definitions.
    # Only the barcode map is rendered (no construct description or metric definitions):
    # which codes were declared, which sample and species each
    # names, and its share of reads, since every accuracy number is scored against it.
    _appendix = [PageBreak(), Paragraph("Barcodes used in this run", heading_style)] + _bobcode_top_elems
    # ILLUMINA: appendix removed (construct/provenance line, auto-detected bobcodes
    # table, Master-Summary metric definitions) — the definitions describe ONT-only
    # rows that no longer exist on this path.
    if not IS_ILLUMINA:
        elems.extend(_appendix)

    doc.build(elems)
    print(f"\nCross-barcode summary saved to: {pdf_path}")
    return pdf_path



if __name__ == '__main__':
    # All experiment-specific values come from EXPERIMENT_CONFIG (loaded from guide)
    barcodes = EXPERIMENT_CONFIG['barcodes']
    # scripts/ is inside the experiment dir, so go up one level to get the experiment root
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    fastq_dir = os.path.join(base_dir, EXPERIMENT_CONFIG['fastq_dir'])
    indiv_dir = os.path.join(base_dir, 'individual-analyses')
    crossbarcode_dir = os.path.join(base_dir, 'crossbarcode')
    os.makedirs(indiv_dir, exist_ok=True)
    os.makedirs(crossbarcode_dir, exist_ok=True)

    # Create dataset-specific sequence list with derivative sequences
    ds_seq_path = create_dataset_sequence_list(
        SEQ_LIST_PATH, base_dir, EXPERIMENT_CONFIG, GROUND_TRUTH, rc)
    print(f"Dataset sequence list: {os.path.basename(ds_seq_path)}")


    bc_overrides = EXPERIMENT_CONFIG.get('barcode_species_overrides', {}) or {}
    if bc_overrides:
        print(f"\nPer-BC species-map overrides:")
        for ov_bc, ov_map in bc_overrides.items():
            print(f"  {ov_bc}: {ov_map}")


    # RENDER-ONLY (per-sample mode): the per-unit analyses were run as their own
    # tasks and left <unit>_speciesmix_results.json behind; this invocation only assembles
    # the cross-unit report from them (the deferred computations below are no-ops once
    # their keys are cached).
    RENDER_ONLY = os.environ.get('BOBSEQ_RENDER_ONLY', '') == '1'
    all_results = {}
    for bc in ([] if RENDER_ONLY else barcodes):
        fq_bc_dir = os.path.join(fastq_dir, bc)
        # Some runs keep the (already adapter-trimmed) reads directly in
        # <base>/<bc>/ with no fastq/ subdir — run_star reads those in place, so
        # analyze must too. Fall back to the barcode dir if fastq/<bc>/ is absent.
        if not os.path.isdir(fq_bc_dir):
            _alt = os.path.join(base_dir, bc)
            if os.path.isdir(_alt):
                fq_bc_dir = _alt
        out_bc_dir = os.path.join(indiv_dir, bc)
        os.makedirs(out_bc_dir, exist_ok=True)

        # File discovery: prefer adapter-trimmed ("noadapter_*.fastq") if present,
        # then fall back to any other *.fastq (skipping "filtered" variants which
        # come from a prior manual filter step via filter_reads.py).
        _all_fq = [f for f in os.listdir(fq_bc_dir) if f.endswith('.fastq') and 'filtered' not in f]
        _noadapter = [f for f in _all_fq if f.startswith('noadapter_')]
        if _noadapter:
            fq_files = _noadapter
        else:
            fq_files = [f for f in _all_fq if not f.startswith('noadapter_')]
        if not fq_files:
            print(f"Skipping {bc}: no FASTQ")
            continue

        # Run analysis: read from fastq/, write to individual-analyses/
        results = analyze_barcode(fq_bc_dir, fq_files[0], output_dir=out_bc_dir)
        results['fastq_dir'] = fq_bc_dir
        # Pre-deduplication metrics, if dedup_reads_illumina.py ran. Loaded HERE rather
        # than in the deferred compute_dup_funnel_all pass because 'Reads (total)' must
        # report reads EXAMINED (pre-filter) on the very first PDF render — otherwise
        # the first PDF shows the post-filter count and only the regenerated one is right.
        _side = os.path.join(out_bc_dir, f'{bc}_dup_prefilter.json')
        if os.path.exists(_side):
            try:
                with open(_side) as _sf:
                    _pf = json.load(_sf)
                if _pf:
                    results['dup_funnel'] = _pf
            except (OSError, ValueError):
                pass
        _ps = os.path.join(out_bc_dir, f'{bc}_prep_stats.json')
        if os.path.exists(_ps):
            try:
                results['prep_stats'] = json.load(open(_ps))
                _persist_result_key(indiv_dir, bc, 'prep_stats', results['prep_stats'])
            except Exception:
                pass
        # No individual per-barcode PDF report. Analysis is pure-STAR (no minimap2
        # deep-correlation step writing {bc}_deep_correlation.json); regen only reads
        # a deep-correlation cache if one happens to be on disk.
        all_results[bc] = results

    # Print summary
    bc_range = f"{barcodes[0]}–{barcodes[-1]}" if len(barcodes) > 1 else barcodes[0]
    print(f"\n{'='*70}")
    loc = EXPERIMENT_CONFIG.get('barcode_location', 'BOBcode')
    primary_map = EXPERIMENT_CONFIG.get('primary_species_map', {})
    primary_codes = list(primary_map.keys())
    bc_type_label = 'BOB' if loc == 'BOBcode' else 'TSO'
    bc_counts_key = 'bob_barcode_counts' if loc == 'BOBcode' else 'tso_barcode_counts'
    bc_detect_key = 'reads_with_bob_c' if loc == 'BOBcode' else 'reads_with_tso'

    # Build column headers from primary species map
    code_headers = [f"{c}(>{sp[:2]})" for c, sp in primary_map.items()]
    code_header_str = ''.join(f'{h:>10}' for h in code_headers)
    print(f"SUMMARY: {bc_type_label} Barcode detection across {bc_range}")
    print(f"{'='*70}")
    print(f"{'Barcode':<12} {'Reads':>6} {bc_type_label+'%':>5} {code_header_str}")
    print("-" * 65)
    for bc in barcodes:
        r = all_results.get(bc, {})
        t = r.get('total_reads', 0)
        detect_n = r.get(bc_detect_key, 0)
        bc_counts = r.get(bc_counts_key, {})
        code_vals = ''.join(f'{bc_counts.get(c, 0):>10}' for c in primary_codes)
        print(f"{bc:<12} {t:>6} {detect_n/max(t,1)*100:>4.0f}% {code_vals}")

    # Enrich every per-barcode JSON with the pure-STAR taxonomy + rRNA-excluded
    # swap so a plain `python3 scripts/analyze_speciesmix_illumina.py` produces the report.
    print("\nEnriching per-barcode JSONs with STAR taxonomy + swap..." if not RENDER_ONLY
          else "\nRender-only: loading the per-unit JSONs written by the per-unit tasks...")
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        # Pure-STAR per-read taxonomy + rRNA-excluded swap (star_taxonomy_illumina.py) ->
        # results['star_taxonomy'] / ['star_swap'] for the 3A-STAR section. Only
        # runs if STAR combined BAMs exist; parses the combined GTF once (heavy)
        # then loops barcodes. Independent of minimap2 outputs.
        try:
            import star_taxonomy_illumina as _stax
            _star_bcs = [bc for bc in ([] if RENDER_ONLY else barcodes)
                         if os.path.exists(os.path.join(
                             indiv_dir, bc, f'{bc}_star_combined.bam'))]
            if _star_bcs:
                print(f"  star_taxonomy: parsing combined GTF for "
                      f"{len(_star_bcs)} barcode(s)...")
                _sgtf = _stax.parse_gtf(_stax.GTF)          # cached after 1st run
                _srdna = _stax.load_rdna(_stax.RDNA_BED)
                # Pre-warm the terminal-exon + adapter caches in the PARENT so the
                # forked workers inherit them (copy-on-write) instead of each
                # re-parsing the GTF. Then fan the per-barcode BAM walks over cores.
                _stax._get_terminal(); _stax._get_adapter()
                _stax._SHARED_GTF = _sgtf; _stax._SHARED_RDNA = _srdna
                _nproc = max(1, min(len(_star_bcs), (os.cpu_count() or 2) - 1))
                if _nproc > 1:
                    import multiprocessing as _mp
                    os.environ.setdefault('OBJC_DISABLE_INITIALIZE_FORK_SAFETY', 'YES')
                    print(f"  star_taxonomy: {_nproc} parallel workers over "
                          f"{len(_star_bcs)} barcode(s)...")
                    try:
                        _ctx = _mp.get_context('fork')
                        with _ctx.Pool(_nproc) as _pool:
                            for _bc, _err in _pool.imap_unordered(
                                    _stax._process_bc_worker, _star_bcs):
                                if _err:
                                    print(f"  {_bc}: star_taxonomy skipped ({_err})")
                    except Exception as _pe:
                        print(f"  parallel pool failed ({_pe}); falling back to serial")
                        for bc in _star_bcs:
                            try: _stax.process_barcode(bc, _sgtf, _srdna)
                            except Exception as e: print(f"  {bc}: skipped ({e})")
                else:
                    for bc in _star_bcs:
                        try:
                            _stax.process_barcode(bc, _sgtf, _srdna)
                        except Exception as e:
                            print(f"  {bc}: star_taxonomy skipped ({e})")
            else:
                print("  star_taxonomy: no *_star_combined.bam found, skipping")
        except Exception as e:
            print(f"  star_taxonomy import skipped: {e}")
        # Reload enriched JSONs into all_results so generate_cross_barcode_summary
        # picks up the STAR fields written above.
        for bc in barcodes:
            json_path = os.path.join(indiv_dir, bc, f'{bc}_speciesmix_results.json')
            if os.path.exists(json_path):
                with open(json_path) as f:
                    all_results[bc] = json.load(f)
    except Exception as e:
        print(f"  STAR enrichment skipped at top level: {e}")

    # Generate the cross-barcode summary PDF (pure-STAR: read structure, STAR
    # RNA-class composition & species-barcode swap, sequence-ref detection, methods).
    # TWO-PHASE: render the PDF first (fast) with the Illumina QC + confident-dup-funnel
    # rows reading "not run yet", THEN run those slow analyses and regenerate the PDF
    # with the values filled in. (Both are no-ops for chemistries that lack P5/P7 or a
    # UMI, in which case the single first PDF is final and the rows read "n/a".)
    generate_cross_barcode_summary(all_results, base_dir)
    # ILLUMINA: skipped. Its P5/P7/i5/i7 page and all 13 of its Master-Summary rows are
    # removed on this path, so the whole FASTQ pass was dead work — measured at 13.3 s of
    # 57.6 s (23%) on an 84k-read subsample, i.e. ~2.6 h at full depth for output nobody
    # renders. Its metrics also need BOTH molecule ends inside one read, which a 150 bp
    # mate never has.
    _did_iq = False if IS_ILLUMINA else compute_illumina_qc_all(all_results, base_dir)
    _did_df = compute_dup_funnel_all(all_results, base_dir)
    _did_rc = compute_read_composition_all(all_results, base_dir)
    _did_rt = compute_rt_extent_all(all_results, base_dir)
    if _did_iq or _did_df or _did_rc or _did_rt:
        print("  regenerating cross-barcode PDF with deferred sections filled in...")
        generate_cross_barcode_summary(all_results, base_dir)

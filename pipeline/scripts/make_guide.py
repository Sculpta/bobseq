#!/usr/bin/env python3
"""Render a run's machine_readable_guide/<run>_analysisguide.txt from its run JSON.

`analyze_speciesmix.load_experiment_config()` parses the guide's BEGIN_CONFIG/END_CONFIG
block; the SOURCE of that block is one JSON per run (config/run_24plex.json for the 24-plex library)
rather than one central registry holding every run.

Why one JSON per run:
  * one file per run diffs cleanly in git; a shared registry makes every run a
    merge conflict
  * JSON is a real type system: `barcodes` is a list, `grun_offsets` is a map of
    int, and a typo fails here rather than silently producing a guide with a
    missing field
  * one config shape with a `chemistry` field covers every chemistry, so there is
    no need for divergent per-chemistry registries

Usage:
    python3 make_guide.py <run.json> [-o <outdir>]

Writes <outdir>/machine_readable_guide/<run>_analysisguide.txt (default outdir: cwd).
"""
import argparse
import json
import os
import sys

# Field order in the emitted guide (fixed so guides are byte-comparable across runs).
_ORDER = ['barcodes', 'fastq_dir', 'experiment_title', 'construct_description',
          'barcode_location', 'tso_arch', 'tso_species_map', 'bob_species_map',
          'grun_offsets', 'unmeasurable', 'rt_primers_used', 'rt_end_label', 'species_design',
          'condition_labels', 'bobcode_labels', 'bobcode_show_max',
          'bobcode_show_frac', 'transcriptome_ref', 'genome_ref', 'mouse_gtf',
          'human_gtf']

# Keys the guide expects as "k=v,k=v" strings.
# 'bobcode_labels' maps a 7-mer to the SAMPLE it came from. condition_labels is keyed by
# ONT barcode, which is enough when one ONT barcode is one sample, but a 24plex pool has
# many samples behind one ONT barcode and they were otherwise nameless in the report.
_MAP_FIELDS = ('tso_species_map', 'bob_species_map', 'grun_offsets', 'condition_labels',
               'bobcode_labels')
# Keys the guide expects as "a,b,c" strings.
# 'unmeasurable' = geometry-dependent elements this run's READ LENGTH cannot observe.
# prep_reads derives the same set from the reads and refuses to run
# unless they agree, so this is an acknowledgement, never a source of truth.
_LIST_FIELDS = ('barcodes', 'rt_primers_used', 'unmeasurable')

REQUIRED = ('run', 'barcodes', 'barcode_location')


def _fmt(key, value):
    """Render a JSON value into the guide's flat string form."""
    if key in _MAP_FIELDS and isinstance(value, dict):
        return ','.join(f'{k}={v}' for k, v in value.items())
    if key in _LIST_FIELDS and isinstance(value, list):
        return ','.join(str(v) for v in value)
    if isinstance(value, bool):
        return 'true' if value else 'false'
    return str(value)


def load_v1_barcodes(sheet):
    """{sequence} for the V1 24plex set, or an empty set if the sheet is unavailable."""
    if not sheet or not os.path.exists(sheet):
        return set()
    out = set()
    with open(sheet) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        try:
            col = header.index('Barcode')
        except ValueError:
            return set()
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            if len(f) > col and f[col]:
                out.add(f[col].strip().upper())
    return out


def check_chemistry(cfg, path, sheet):
    """Guard the one config mistake that silently produces wrong numbers.

    `barcode_adapter.is_24plex` prefers the explicit `tso_arch` field and otherwise
    falls back to a 6-member allowlist of 24plex 7-mers -- a list that predates the
    field and covers only a quarter of the V1 set. So a 24plex run whose barcodes are
    outside those six (e.g. TGGATGA/TGTCCAC) and which omits `tso_arch` is parsed with
    the CLASSIC tso-bob adapter instead: wrong structure model, wrong barcode accuracy,
    no error anywhere. Copying a run config and dropping one field is enough.

    So: if the species map contains a known V1 barcode, `tso_arch` is mandatory.
    """
    if cfg.get('barcode_location') != 'TSO':
        return
    smap = cfg.get('tso_species_map') or {}
    codes = {str(c).upper() for c in smap}
    v1 = load_v1_barcodes(sheet)
    arch = (cfg.get('tso_arch') or '').strip()

    if not arch:
        hits = sorted(codes & v1)
        if hits:
            sys.exit(
                f"make_guide: {path} uses V1 24plex barcode(s) {', '.join(hits)} but does "
                f"not set 'tso_arch'.\n"
                f"  Without it the analyzer falls back to a 6-member allowlist that does "
                f"not cover the whole V1 set,\n"
                f"  silently parsing this run with the classic tso-bob chemistry and "
                f"producing wrong barcode accuracy.\n"
                f'  Add:  "tso_arch": "24plex"')
        print(f"make_guide: NOTE no 'tso_arch' set; chemistry will be inferred from the "
              f"barcode sequences. Set it explicitly if this is a 24plex run.")
        return

    if arch.lower() in ('24plex', '24-plex') and v1:
        unknown = sorted(codes - v1)
        if unknown:
            print(f"make_guide: WARNING barcode(s) {', '.join(unknown)} are not in the V1 "
                  f"24plex registry (assets/24plex_bobcode_set_v1.tsv).")
            print(f"  Alignment and barcode accuracy are driven by this config and are "
                  f"unaffected, but metrics that SCAN for known V1 motifs -- "
                  f"'Empty products %' and 'Poly(dT)-less artifacts %' -- only count "
                  f"reads carrying a registered barcode, so they will under-report.")
            print(f"  If this is a new bobcode set, add it to the registry.")


_PRIMING_TOKENS = {
    'polydt':      ('polydt', 'polyt', '_dt', '-dt'),
    # '-6n'/'_6n' rather than a bare '6n': the polydT primer is named R1_26N_polydT, which
    # CONTAINS '6n', so a bare token would make it match both chemistries, _priming_from_
    # rt_primers would see two hits and return None, and the consistency guard would stop
    # recognising every polydT run. The separator makes the 6N random primer matchable
    # without that collision.
    'randomprime': ('9n', '-6n', '_6n', 'random', 'hexamer'),
}


def _priming_from_rt_primers(rt_primers):
    """Which priming chemistry the RT primer names imply, or None if not recognised."""
    joined = ' '.join(str(p) for p in (rt_primers or [])).lower()
    hits = {name for name, toks in _PRIMING_TOKENS.items() if any(t in joined for t in toks)}
    return hits.pop() if len(hits) == 1 else None


def check_species_sheet(cfg, path):
    """Every bobcode must carry an explicit species from THIS run's sample sheet, and the
    species map and the sample labels must name the same codes. A code that appears in one
    and not the other is a sheet/config mismatch, and a code without a species would otherwise
    silently inherit one from an earlier run (the same code can be mouse in one run and
    human in another)."""
    if cfg.get('barcode_location') != 'TSO':
        return
    smap = {str(k).upper(): str(v).lower() for k, v in (cfg.get('tso_species_map') or {}).items()}
    labels = {str(k).upper() for k in (cfg.get('bobcode_labels') or {})}
    bad = [k for k, v in smap.items() if v not in ('human', 'mouse')]
    if bad:
        sys.exit(f"make_guide: {path}: tso_species_map has non human/mouse species for {bad}")
    if labels:
        only_labels = sorted(labels - set(smap)); only_map = sorted(set(smap) - labels)
        if only_labels or only_map:
            sys.exit(f"make_guide: {path}: sample sheet mismatch -- bobcodes in bobcode_labels but not "
                     f"tso_species_map: {only_labels}; in tso_species_map but not bobcode_labels: {only_map}. "
                     f"Every sample needs an explicit species from this run's sheet.")
    print(f"make_guide: species sheet ok: {len(smap)} bobcodes, "
          f"{sum(1 for v in smap.values() if v == 'human')} human / {sum(1 for v in smap.values() if v == 'mouse')} mouse")


def check_chemistry_consistency(cfg, path):
    """`chemistry` is DECLARATIVE -- nothing selects behaviour from it.

    The adapter is chosen by barcode_location + tso_arch (barcode_adapter.get_adapter),
    and the priming chemistry shows up only through rt_primers_used. So a config can say
    "24plex-polydt" while actually being analysed as random-priming, and no stage will
    disagree. Nothing crashes; the run is simply filed under the wrong chemistry.

    That matters because `chemistry` is what the cross-run matrix groups by
    (aggregate_summaries) and what comparable_funnel stamps into its output. A mislabelled
    run is compared against the wrong population, which is precisely the error the
    chemistry-general funnel exists to prevent.

    Only a CONFIDENT contradiction fails. An unrecognised chemistry string or an
    unfamiliar RT primer is a note, not an error -- new chemistries are expected, and this
    guard must not be the thing blocking one.
    """
    declared = (cfg.get('chemistry') or '').strip()
    if not declared:
        return

    arch_part, _, priming_part = declared.lower().partition('-')

    # arch half vs the field that actually drives adapter selection
    tso_arch = (cfg.get('tso_arch') or '').strip().lower().replace('-', '')
    if tso_arch and arch_part.replace('-', '') and arch_part.replace('-', '') != tso_arch:
        sys.exit(
            f"make_guide: {path} declares chemistry={declared!r} but tso_arch={cfg['tso_arch']!r}.\n"
            f"  tso_arch is what selects the analysis adapter; chemistry is only a label.\n"
            f"  They disagree, so the run would be analysed as one chemistry and reported "
            f"as another.")

    # priming half vs the RT primers. Both sides go through the SAME matcher, so a
    # contradiction is only ever reported between two things this guard actually
    # recognises -- an unfamiliar chemistry name is unknown, not wrong, and must not be
    # the thing that blocks a new chemistry.
    implied = _priming_from_rt_primers(cfg.get('rt_primers_used'))
    declared_priming = _priming_from_rt_primers([priming_part]) if priming_part else None
    if declared_priming and implied and declared_priming != implied:
        sys.exit(
            f"make_guide: {path} declares chemistry={declared!r} but "
            f"rt_primers_used={cfg.get('rt_primers_used')!r} implies {implied!r}.\n"
            f"  Nothing downstream reconciles these: the analysis follows the primers, the "
            f"cross-run matrix groups by chemistry.\n"
            f"  Fix whichever is wrong -- a mislabelled run is compared against the wrong "
            f"population.")

    if priming_part and not (declared_priming and implied):
        print(f"make_guide: NOTE chemistry={declared!r} declares priming {priming_part!r} "
              f"with rt_primers_used={cfg.get('rt_primers_used')!r}; this guard recognises "
              f"one side but not both, so the pair was NOT cross-checked. Not an error -- "
              f"add the primer token to _PRIMING_TOKENS if this chemistry is here to stay.")


def validate(cfg, path, sheet=None):
    missing = [k for k in REQUIRED if not cfg.get(k)]
    if missing:
        sys.exit(f"make_guide: {path} is missing required field(s): {', '.join(missing)}")

    if cfg['barcode_location'] not in ('TSO', 'BOBcode'):
        sys.exit(f"make_guide: barcode_location must be 'TSO' or 'BOBcode', "
                 f"got {cfg['barcode_location']!r}")

    smap_key = 'tso_species_map' if cfg['barcode_location'] == 'TSO' else 'bob_species_map'
    smap = cfg.get(smap_key) or {}
    if not smap:
        sys.exit(f"make_guide: {path} has barcode_location={cfg['barcode_location']} "
                 f"but no {smap_key}")
    for code, species in smap.items():
        if species not in ('human', 'mouse'):
            sys.exit(f"make_guide: {smap_key}[{code}] must be 'human' or 'mouse', "
                     f"got {species!r}")
        if not set(code) <= set('ACGT'):
            sys.exit(f"make_guide: {smap_key} key {code!r} is not a DNA sequence")

    # Every barcode in the species map should carry a GGG offset for the 24plex
    # distance step; a missing one silently falls back to 8 and skews that funnel level.
    offsets = cfg.get('grun_offsets') or {}
    if cfg.get('tso_arch') == '24plex':
        no_off = [c for c in smap if c not in offsets]
        if no_off:
            sys.exit(f"make_guide: tso_arch=24plex requires grun_offsets for every "
                     f"species barcode; missing: {', '.join(no_off)}")

    # Condition labels are per-barcode; a stale one usually means a copy-paste run config.
    labels = cfg.get('condition_labels') or {}
    stray = [b for b in labels if b not in cfg['barcodes']]
    if stray:
        sys.exit(f"make_guide: condition_labels refer to barcodes not in this run: "
                 f"{', '.join(stray)}")

    check_chemistry(cfg, path, sheet)

    check_species_sheet(cfg, path)
    check_chemistry_consistency(cfg, path)


def build(cfg):
    """Return the guide's key->string mapping from the run config."""
    out = {}
    for k, v in cfg.items():
        if k in ('run', 'chemistry', 'notes'):     # run-config-only keys, not guide fields
            continue
        out[k] = _fmt(k, v)
    out.setdefault('fastq_dir', 'fastq')
    # Both maps are always set to the same value so code paths reading either key
    # behave identically.
    if cfg.get('barcode_location') == 'TSO' and 'bob_species_map' not in out:
        out['bob_species_map'] = out.get('tso_species_map', '')
    return out


def render(run, fields):
    lines = [
        f"# Auto-generated by make_guide.py for run '{run}'.",
        "# DO NOT EDIT here -- edit runs/<run>.json and re-run the pipeline.",
        "# The master sequence list is config/sequence_list.txt.",
        "",
        "BEGIN_CONFIG",
    ]
    ordered = [k for k in _ORDER if k in fields] + [k for k in fields if k not in _ORDER]
    for k in ordered:
        lines.append(f"{k}\t{fields[k]}")
    lines.append("END_CONFIG")
    return '\n'.join(lines) + '\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_json')
    ap.add_argument('-o', '--outdir', default='.')
    ap.add_argument('--bobcodes', default=None,
                    help='24plex bobcode registry TSV (default: <repo>/assets/24plex_bobcode_set_v1.tsv)')
    args = ap.parse_args()

    with open(args.run_json) as fh:
        cfg = json.load(fh)
    sheet = args.bobcodes or os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'assets', '24plex_bobcode_set_v1.tsv')
    validate(cfg, args.run_json, sheet)

    run = cfg['run']
    mrg = os.path.join(args.outdir, 'machine_readable_guide')
    os.makedirs(mrg, exist_ok=True)
    out = os.path.join(mrg, f'{run}_analysisguide.txt')
    with open(out, 'w') as fh:
        fh.write(render(run, build(cfg)))

    smap = cfg.get('tso_species_map') or cfg.get('bob_species_map') or {}
    print(f"make_guide: wrote {out}")
    print(f"  run={run}  barcodes={len(cfg['barcodes'])}  "
          f"location={cfg['barcode_location']}  arch={cfg.get('tso_arch', '(unset)')}")
    print(f"  species map: {smap}  (PROVISIONAL -- confirmed against the alignment "
          f"by detect_species_map)")


if __name__ == '__main__':
    main()

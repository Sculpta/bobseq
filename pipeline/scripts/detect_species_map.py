#!/usr/bin/env python3
"""Auto-detect the 7-mer -> species DIRECTION from the STAR alignment and check it
against the run's declared map (from the registry/guide).

WHY this is auto-detected, not trusted: the 7-mer->species direction INVERTS from
run to run (e.g. GATATGG=mouse in one run, =human in another). A wrong direction
silently drops barcode accuracy to ~1-3% while everything else looks fine. So the
registry's `tso_species_map` is only a PROVISIONAL/expected value; the alignment is
ground truth. Each barcode read carries a 7-mer (in the sequence) and aligns to a
species (chromosome prefix) — the per-code human/mouse tally (persisted by
star_taxonomy, independent of the declared direction) tells us the true direction.

  *** APPLIED TO ALL CHEMISTRIES ***
  This calls the direction from `star_crosstab_mtrep` (rRNA + mitochondrial
  reads EXCLUDED), NOT the all-reads `star_crosstab`. Random priming has no polyA
  enrichment, so 42-80% of mapped reads are rRNA; rRNA is highly conserved and
  cross-maps to human on the combined index, which drags the all-reads tally toward
  human and falsely flags the mouse 7-mer as LOW-CONFIDENCE (measured on one
  random-priming library: all-reads 66% human -> WRONG; mtrep 97% mouse ->
  correct). Excluding rRNA/MT is strictly safer for every chemistry (polydT rRNA is
  low, so mtrep ~ all there), so it is applied everywhere: calling from the
  all-reads tally is what would let a confident-MISMATCH auto-fix INVERT a correct
  species map on an rRNA-heavy run. Per code, if mtrep has < N_MIN reads we fall
  back to the all-reads tally and say so, so low-depth runs still get a call.

Reads run_dir/individual-analyses/*/*_results.json, aggregates the crosstab across
barcodes, and for each code reports the dominant aligned species + confidence.
Compares to the guide's tso_species_map. Writes species_map_check.txt in the run
dir. Exit codes: 0 = matches declared (confirmed), 2 = MISMATCH (declared direction
is wrong -> fix the guide/registry), 3 = low-confidence/ambiguous (needs a human
look). With --fix-guide, on a confident mismatch it rewrites the guide's
tso_species_map/bob_species_map to the detected direction (so a re-run scores right).

Usage:  python3 detect_species_map.py <run_dir> [--fix-guide]
"""
import glob
import json
import os
import sys

CONF_MIN = 0.70   # dominant-species fraction to call a direction confidently
N_MIN = 50        # min human+mouse reads for a code to be callable


def _guide_path(run_dir):
    mrg = os.path.join(run_dir, 'machine_readable_guide')
    if not os.path.isdir(mrg):
        return None
    for f in sorted(os.listdir(mrg)):
        if f.endswith('_analysisguide.txt') or f.endswith('_guide.txt'):
            return os.path.join(mrg, f)
    return None


def _parse_map(guide_path, key='tso_species_map'):
    if not guide_path or not os.path.exists(guide_path):
        return {}
    infl = False
    for ln in open(guide_path):
        s = ln.strip()
        if s == 'BEGIN_CONFIG':
            infl = True
            continue
        if s == 'END_CONFIG':
            break
        if infl and '\t' in ln:
            k, v = ln.rstrip('\n').split('\t', 1)
            if k.strip() == key:
                out = {}
                for pair in v.split(','):
                    if '=' in pair:
                        code, sp = pair.split('=', 1)
                        out[code.strip()] = sp.strip()
                return out
    return {}


def _agg(run_dir, key):
    """Aggregate a persisted crosstab (by key) -> {code: {'human':n,'mouse':n}}."""
    agg = {}
    for jp in sorted(glob.glob(os.path.join(run_dir, 'individual-analyses', '*', '*_results.json'))):
        ct = (json.load(open(jp)) or {}).get(key) or {}
        for code, sp in ct.items():
            d = agg.setdefault(code, {'human': 0, 'mouse': 0})
            d['human'] += sp.get('human', 0)
            d['mouse'] += sp.get('mouse', 0)
    return agg


def detect(run_dir):
    """RANDOM-PRIMING: call direction from the rRNA+MT-excluded crosstab.

    Returns {code: {'human':n,'mouse':n,'src':'mtrep'|'all'}}. Per code we prefer
    the rRNA/MT-excluded tally (`star_crosstab_mtrep`); if that code has fewer than
    N_MIN reads there we fall back to the all-reads tally (`star_crosstab`) so a
    low-depth code is still callable — the source used is recorded in 'src'.
    """
    mtrep = _agg(run_dir, 'star_crosstab_mtrep')
    allx = _agg(run_dir, 'star_crosstab')
    out = {}
    for code in set(mtrep) | set(allx):
        mt = mtrep.get(code, {'human': 0, 'mouse': 0})
        if mt['human'] + mt['mouse'] >= N_MIN:
            out[code] = {**mt, 'src': 'mtrep'}
        else:
            al = allx.get(code, {'human': 0, 'mouse': 0})
            out[code] = {**al, 'src': 'all'}
    return out


def main():
    if len(sys.argv) < 2:
        sys.exit('usage: detect_species_map.py <run_dir> [--fix-guide]')
    run_dir = os.path.abspath(sys.argv[1])
    fix = '--fix-guide' in sys.argv[2:]
    guide = _guide_path(run_dir)
    declared = _parse_map(guide, 'tso_species_map') or _parse_map(guide, 'bob_species_map')
    agg = detect(run_dir)
    if not agg:
        sys.exit('detect_species_map: no star_crosstab found — run analyze first')

    # Two independent conditions, tracked separately and combined at the END:
    #   mismatch  -> a code whose declared direction is confidently WRONG (exit 2)
    #   lowconf   -> a code that is thin, low-confidence or undeclared (exit 3)
    # Folding them into one `status = max(status, ...)` would let 3 outrank 2: a run
    # with ONE inverted code and ONE thin code would exit 3, the `--fix-guide` branch
    # (gated on status == 2) would never run, and the run would publish with the
    # wrong direction behind a "needs a human look" banner. A confident mismatch is
    # the more serious finding and must win.
    lines, mismatch, lowconf = [], False, False
    detected = {}
    hdr = (f"{'7-mer':<12}{'aligned human/mouse':>22}{'src':>7}{'detected':>10}"
           f"{'conf':>7}{'declared':>10}   verdict")
    lines.append(hdr)
    for code in sorted(agg):
        h, m = agg[code]['human'], agg[code]['mouse']
        src = agg[code].get('src', 'mtrep')   # 'mtrep' = rRNA/MT-excluded, 'all' = fallback
        n = h + m
        if n < N_MIN:
            det, conf, verdict = '?', 0.0, f'too few reads (n={n})'
            lowconf = True
        else:
            det = 'human' if h >= m else 'mouse'
            conf = max(h, m) / n
            decl = declared.get(code)
            if conf < CONF_MIN:
                # NOT recorded in `detected`: a low-confidence call must never be
                # merged into the guide by --fix-guide.
                verdict = f'LOW CONFIDENCE ({conf:.0%}) — check'
                lowconf = True
            elif decl is None:
                detected[code] = det
                verdict = 'no declared value'
                lowconf = True
            elif det == decl:
                detected[code] = det
                verdict = 'OK — matches declared'
            else:
                detected[code] = det
                verdict = f'MISMATCH — declared {decl}, alignment says {det}'
                mismatch = True
            if src == 'all':
                verdict += ' [rRNA-excluded tally too thin — fell back to all-reads]'
        lines.append(f"{code:<12}{f'{h:,}/{m:,}':>22}{src:>7}{det:>10}{conf:>6.0%}"
                     f"{declared.get(code, '-'):>10}   {verdict}")

    # A confident mismatch takes precedence over any low-confidence code.
    status = 2 if mismatch else (3 if lowconf else 0)
    banner = {0: '✓ species direction CONFIRMED by alignment (matches the registry).',
              2: '⚠ species direction MISMATCH — the registry/guide direction is WRONG for this run.',
              3: '⚠ species direction needs a human look (low confidence / missing declared value).'}[status]
    if mismatch and lowconf:
        banner += '\n  (some codes are ALSO low-confidence or thin; see rows — they are not auto-corrected)'
    report = ['SPECIES-DIRECTION AUTO-DETECTION (7-mer -> species, from STAR alignment)',
              '  [called from rRNA+MT-EXCLUDED reads (star_crosstab_mtrep);',
              '   conserved rRNA would otherwise cross-map to human and muddy the call]',
              '='*78, banner, '',
              'NOTE: on 24-plex runs the rows below are species CLASSES (all human codes /',
              '  all mouse codes) labelled by the legacy representative 7-mer, not individual',
              '  bobcodes. The per-bobcode check against the sample sheet is validate_outputs',
              '  invariant 12 (per-sample species direction).', '',
              'The direction is auto-detected here and should be CONFIRMED by you:',
              '  - if it matches the registry, nothing to do;',
              '  - if it MISMATCHES, update run_registry.txt [<run>] tso_species_map to the',
              '    detected direction (this run\'s barcode accuracy is only trustworthy once',
              '    the declared direction matches the alignment).', '',
              *lines, '']
    out = '\n'.join(report)
    print(out)
    with open(os.path.join(run_dir, 'species_map_check.txt'), 'w') as fh:
        fh.write(out + '\n')

    if fix and status == 2 and detected and guide:
        # MERGE the detected direction into the declared map. Replacing the map with
        # only the detected codes silently dropped every code the library-level
        # crosstab does not tally (24-plex runs collapse to two class rows), so a
        # re-run lost 22 of 24 samples. Codes without a detection keep their sheet.
        merged = {**declared, **detected}
        new = ','.join(f'{c}={merged[c]}' for c in sorted(merged))
        # Rewrite ONLY inside the BEGIN_CONFIG/END_CONFIG block that _parse_map
        # reads; a bare whole-file substitution also hit any prose line that
        # happened to start with the key.
        out_lines, in_cfg, n_sub = [], False, 0
        for ln in open(guide):
            s = ln.strip()
            if s == 'BEGIN_CONFIG':
                in_cfg = True
            elif s == 'END_CONFIG':
                in_cfg = False
            elif in_cfg:
                key = ln.split('\t', 1)[0]
                if key in ('tso_species_map', 'bob_species_map'):
                    ln = f'{key}\t{new}\n'
                    n_sub += 1
            out_lines.append(ln)
        if n_sub == 0:
            sys.exit('detect_species_map: --fix-guide found no species-map line inside '
                     'the BEGIN_CONFIG/END_CONFIG block; guide NOT rewritten')
        open(guide, 'w').writelines(out_lines)
        print(f"detect_species_map: --fix-guide rewrote the guide direction -> {new}")
    sys.exit(status)


if __name__ == '__main__':
    main()

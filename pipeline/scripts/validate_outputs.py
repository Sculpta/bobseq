#!/usr/bin/env python3
"""Post-run invariants for the pipeline. Exits non-zero when an output is internally
inconsistent, so a run FAILS instead of publishing quietly.

Each check guards a failure mode that would otherwise ship silently:

  sample-count   a whole library dropped by a pipeline wiring mistake while the
                 run reports SUCCESS.
  libs-sum       the per-library split must reconcile with the pooled table; a
                 dropped library is invisible in the pooled numbers alone.
  composition    a pooled composition with a mapped-only denominator, which
                 zeroes Unmapped and inflates every other class.
  conservation   raw == assigned + unclassified + undeclared, per ONT barcode.

TWO MODES, so single-library runs are checked as well as multiplexed ones:

  --mode library     runs on EVERY run, after the analysis and before the run
                     is published:
                     L1 the Master Summary carries exactly the declared barcodes
                     L2 the species-direction check was produced and is readable
                     L3 (Illumina) prep read accounting conserves: raw = written +
                        dropped, and the sidecar / published FASTQ agree with it
                     9/10/11 on the library sidecar and published FASTQ (Illumina)
  --mode per-sample  invariants 1-12 over the per-sample tables.

INPUTS ARE REQUIRED PER MODE AND PLATFORM. An input that changes the meaning of a
check is never optional ("absent = not checked" lets a missing file degrade a
check instead of failing it). A missing required input FAILS unless
--allow-missing is given, which prints what was skipped in capitals.

CLI: validate_outputs.py --mode library    --run-json R --master-summary-tsv T
                                            --species-check S [--prep-stats P
                                            --library-sidecar L --prepped-fastq F]
     validate_outputs.py --mode per-sample --run-json R --composition C
                                            --demux-stats D... --per-sample P
                                            [--library-composition L...]
                                            [--prep-stats P --sidecars S...
                                            --library-sidecar L --prepped-fastq F]
"""
import argparse, csv, gzip, json, os, sys

CATS = ['mRNA (protein-coding)', 'Ribosomal-protein mRNA',
        'Other ncRNA (lncRNA/NMD/ret-intron/pseudo)', 'rRNA', 'Mitochondrial',
        'pre-mRNA / intronic', 'Intergenic / genomic', 'Unmapped']


def load(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith('#')), delimiter='\t'))


def _fastq_reads(path):
    n = 0
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as fh:
        for n, _ in enumerate(fh, 1):
            pass
    return n // 4


def _undeclared_map(s):
    """The undeclared-bobcode tally. The ONT demux writes `undeclared_bobcodes`,
    the Illumina demux (demux_qname_illumina.py) writes `undeclared`; both keys
    are read so the check works on either platform. Returns (mapping, present)."""
    for k in ('undeclared', 'undeclared_bobcodes'):
        if isinstance(s.get(k), dict):
            return s[k], True
    return {}, False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', choices=('per-sample', 'library'), default='per-sample')
    ap.add_argument('--run-json', required=True)
    # library mode
    ap.add_argument('--master-summary-tsv', default=None)
    ap.add_argument('--species-check', default=None)
    # per-sample mode
    ap.add_argument('--composition', default=None)
    ap.add_argument('--library-composition', nargs='*', default=[])
    ap.add_argument('--demux-stats', nargs='*', default=[])
    ap.add_argument('--per-sample', default=None)
    # Illumina extras (both modes)
    ap.add_argument('--prep-stats', default=None)
    ap.add_argument('--sidecars', nargs='*', default=[])
    ap.add_argument('--library-sidecar', default=None)
    ap.add_argument('--prepped-fastq', default=None)
    ap.add_argument('--allow-missing', action='store_true',
                    help='downgrade a missing REQUIRED input to a loud note (escape hatch)')
    # Read floor per declared bobcode. PROPORTIONAL by default: a declared code with less
    # than this fraction of the library's assigned reads is "absent" (the run json names a
    # sample that is not in the pool, or a lane was dropped). An absolute floor (e.g.
    # 1000) would fail every shallow fixture and any low-yield sample outright; pass
    # --min-reads-per-sample to add one on top when a run warrants it.
    ap.add_argument('--min-frac-per-sample', type=float, default=0.001)
    ap.add_argument('--min-reads-per-sample', type=int, default=0)
    ap.add_argument('--max-undeclared-frac', type=float, default=0.005)
    ap.add_argument('--dedup-tolerance', type=float, default=0.02)
    ap.add_argument('--tolerance', type=float, default=0.05,
                    help='max |100 - sum(categories)| in percentage points')
    args = ap.parse_args()

    fail = []
    note = lambda m: sys.stderr.write(f'  {m}\n')
    cfg = json.load(open(args.run_json))
    illumina = str(cfg.get('platform') or '').lower().startswith('illumina') \
        or cfg.get('read_layout') == 'r2_positional'
    note(f'mode={args.mode} platform={"illumina" if illumina else "ont"}')

    # ---- required inputs per mode/platform -------------------------------------
    req = {'library': ['master_summary_tsv', 'species_check'],
           'per-sample': ['composition', 'demux_stats', 'per_sample']}[args.mode]
    if illumina:
        req += ['prep_stats', 'library_sidecar', 'prepped_fastq']
        if args.mode == 'per-sample':
            req.append('sidecars')
    missing = [r for r in req if not getattr(args, r)]
    present = lambda r: bool(getattr(args, r)) and r not in missing
    if missing:
        msg = f'REQUIRED input(s) not provided for --mode {args.mode}: {missing}'
        if args.allow_missing:
            note('SKIPPED (--allow-missing): ' + msg)
        else:
            fail.append(msg + ' (pass --allow-missing to skip, loudly)')
    for r in req:
        v = getattr(args, r)
        for p in (v if isinstance(v, list) else [v] if v else []):
            if not os.path.exists(p):
                fail.append(f'input {r}: file not found: {p}')

    ps = json.load(open(args.prep_stats)) if present('prep_stats') and os.path.exists(args.prep_stats) else None
    lib_side = json.load(open(args.library_sidecar)) if present('library_sidecar') and os.path.exists(args.library_sidecar) else None

    # =========================================================================
    # LIBRARY MODE
    # =========================================================================
    if args.mode == 'library':
        declared_bcs = list(cfg.get('barcodes') or [])
        # L1. the Master Summary carries exactly the declared units. A dropped
        # library shows here as a missing column.
        if present('master_summary_tsv') and os.path.exists(args.master_summary_tsv):
            with open(args.master_summary_tsv) as fh:
                hdr = fh.readline().rstrip('\n').split('\t')
                nrows = sum(1 for _ in fh)
            cols = hdr[1:]
            if sorted(cols) != sorted(declared_bcs):
                fail.append(f'L1 Master Summary columns {cols} != declared barcodes {declared_bcs}')
            if nrows == 0:
                fail.append('L1 Master Summary has no metric rows')
            note(f'L1 master summary: {len(cols)} column(s), {nrows} rows, declared {len(declared_bcs)}')

        # L2. the species-direction sign-off exists and was actually produced.
        if present('species_check') and os.path.exists(args.species_check):
            txt = open(args.species_check).read()
            if 'not produced' in txt.lower():
                fail.append('L2 species_map_check.txt says the check was NOT produced: '
                            'barcode accuracy is unconfirmed')
            banner = next((l for l in txt.splitlines() if 'species direction' in l), '(no banner)')
            if 'MISMATCH' in banner:
                fail.append('L2 species direction still reports MISMATCH after analysis: '
                            'the self-correction did not take')
            note(f'L2 species check: {banner.strip()}')

        # L3. Illumina prep accounting conserves, and the downstream files agree with it.
        if ps:
            n = ps.get('n', 0)
            parts = (ps.get('written', 0) + ps.get('no_bobcode', 0)
                     + ps.get('too_short_after_trim', 0) + ps.get('strided_out', 0))
            if n != parts:
                fail.append(f'L3 prep accounting: n={n} != written+no_bobcode+too_short+strided_out={parts}')
            called = ps.get('bobcode_called', 0)
            if called != ps.get('written', 0) + ps.get('too_short_after_trim', 0):
                fail.append(f'L3 prep accounting: bobcode_called={called} != written+too_short='
                            f'{ps.get("written", 0) + ps.get("too_short_after_trim", 0)}')
            if ps.get('strided_out'):
                fail.append(f'L3 prep ran with a stride ({ps["strided_out"]} pairs skipped): '
                            'every percentage in prep_stats is computed over ALL pairs and is wrong')
            note(f'L3 prep accounting: n={n:,} written={ps.get("written", 0):,} '
                 f'no_bobcode={ps.get("no_bobcode", 0):,} too_short={ps.get("too_short_after_trim", 0):,}')
        if present('prepped_fastq') and ps and os.path.exists(args.prepped_fastq):
            nfq = _fastq_reads(args.prepped_fastq)
            if nfq != ps.get('written'):
                fail.append(f'11 published prepped FASTQ has {nfq} reads, prep wrote {ps.get("written")}')
            note(f'11 published input: {nfq:,} reads == prep written')
        if lib_side:
            _check_sidecar('library', lib_side, cfg, fail, note)
            n_in = lib_side.get('fastq_reads', lib_side.get('n_reads', 0))
            kept = lib_side.get('dedup_kept', 0)
            if kept > n_in:
                fail.append(f'10 library: dedup kept {kept} > input {n_in}')
            if ps and n_in > ps.get('written', 0):
                fail.append(f'10 library: sidecar input {n_in} > prep written {ps.get("written")}')
            note(f'10 dedup sanity: kept {kept:,} of {n_in:,}')
        _finish(fail)
        return

    # =========================================================================
    # PER-SAMPLE MODE (invariants 1-12)
    # =========================================================================
    # unlabeled run: the sample IS the code (the same rule as demux and the per-sample tables)
    declared = cfg.get('bobcode_labels') or {c: c for c in (cfg.get('tso_species_map') or {})}
    expect_n = len(declared)
    comp = load(args.composition) if present('composition') and os.path.exists(args.composition) else []

    # 1. every declared sample present exactly once
    names = [r['sample'] for r in comp]
    want = set(declared.values())
    if len(comp) != expect_n:
        fail.append(f'sample count {len(comp)} != {expect_n} declared in the run config')
    missing_s = sorted(want - set(names))
    if missing_s:
        fail.append(f'declared samples absent from the composition table: {missing_s}')
    dupes = sorted({n for n in names if names.count(n) > 1})
    if dupes:
        fail.append(f'duplicate sample rows: {dupes}')
    note(f'samples: {len(comp)} rows, {expect_n} declared, missing={len(missing_s)}')

    # 2. composition is exhaustive (scan every row; report the worst)
    worst = 0.0; worst_s = None; miss_rows = []
    for r in comp:
        miss = [c for c in CATS if c not in r]
        if miss:
            miss_rows.append((r['sample'], miss)); continue
        d = abs(100.0 - sum(float(r[c]) for c in CATS))
        if d > worst: worst, worst_s = d, r['sample']
    if miss_rows:
        fail.append(f'composition missing categories in {len(miss_rows)} row(s), e.g. {miss_rows[0]}')
    if worst > args.tolerance:
        fail.append(f'composition does not sum to 100 (worst {worst:.3f} pts, {worst_s})')
    note(f'composition: max |100-sum| = {worst:.4f} pts')

    # 3. the mapped-only-denominator regression. ANY sample with assigned == classified
    #    and exactly zero unmapped is degenerate (real data never has exactly 0%); it is
    #    checked per sample so one bad lane cannot pass.
    degenerate = [r['sample'] for r in comp
                  if r.get('reads_assigned') == r.get('reads_classified')
                  and float(r.get('Unmapped', 0)) == 0.0]
    if degenerate:
        fail.append(f'{len(degenerate)}/{len(comp)} sample(s) have reads_assigned == '
                    f'reads_classified and Unmapped = 0 (mapped-only denominator): '
                    f'{degenerate[:5]}')
    note(f'denominator: {len(comp)-len(degenerate)}/{len(comp)} samples carry an unmapped fraction')

    # 4. per-library tables reconcile with the pooled one
    if args.library_composition:
        pooled = {r['sample']: int(r['reads_classified']) for r in comp}
        summed = {}
        for p in args.library_composition:
            for r in load(p):
                summed[r['sample']] = summed.get(r['sample'], 0) + int(r['reads_classified'])
        bad = {s: (summed.get(s, 0), pooled[s]) for s in pooled if summed.get(s, 0) != pooled[s]}
        if bad:
            ex = list(bad.items())[:3]
            fail.append(f'per-library classified reads do not sum to pooled for '
                        f'{len(bad)} sample(s), e.g. {ex}')
        same = [p for p in args.library_composition
                if os.path.realpath(p) == os.path.realpath(args.composition or '')]
        note(f'libs-sum: {len(pooled)-len(bad)}/{len(pooled)} samples reconcile'
             + (' (single library: pooled == per-library by construction)' if same else ''))

    # 5. demux conservation
    dm = []
    for d in args.demux_stats:
        if not os.path.exists(d):
            fail.append(f'demux stats not found: {d}'); continue
        j = json.load(open(d)); dm.append(j)
        tot = j['reads_total']
        parts = j['reads_assigned'] + j['reads_unclassified'] + j.get('reads_undeclared', 0)
        if tot != parts:
            fail.append(f'{os.path.basename(d)}: {tot} != assigned+unclassified+undeclared '
                        f'({parts}), diff {tot-parts}')
        note(f'conservation {j.get("ont_barcode") or os.path.basename(d)}: '
             f'{tot:,} == {parts:,}  {"OK" if tot==parts else "FAIL"}')

    # 6. the species-accuracy column must be the rDNA-EXCLUDED one on a two-species run.
    if present('per_sample') and os.path.exists(args.per_sample):
        psr = load(args.per_sample)
        cols = set(psr[0]) if psr else set()
        need = {'pct_expected_species_no_rdna',
                'pct_expected_species_ALLREADS_rRNA_INFLATED', 'human_rdna', 'mouse_rdna'}
        absent = sorted(need - cols)
        if absent:
            fail.append(f'per_sample table missing rDNA-aware columns {absent}: '
                        f'per_sample_table.py ran without --rdna-bed')
        elif len({r['species'] for r in psr if r.get('species')}) > 1:
            rd = sum(int(r['human_rdna'] or 0) + int(r['mouse_rdna'] or 0) for r in psr)
            if rd == 0:
                fail.append('two-species run but ZERO rDNA reads counted: the rDNA BED did '
                            'not match the reference contig names, so the corrected '
                            'species-accuracy column is silently the inflated one')
            note(f'species-accuracy: rDNA-aware, {rd:,} rDNA reads excluded')

    side = {os.path.basename(f).replace('_dup_prefilter.json', ''): json.load(open(f))
            for f in (args.sidecars or []) if os.path.exists(f)}

    # 7. conservation through the split
    for s in dm:
        parts = s['reads_assigned'] + s.get('reads_unclassified', 0) + s.get('reads_undeclared', 0)
        if parts != s['reads_total']:
            fail.append(f'7 split conservation: reads_total {s["reads_total"]} != assigned+unclassified+undeclared {parts}')
        if ps and ps.get('written') is not None and s['reads_total'] != ps['written']:
            fail.append(f'7 split conservation: demux saw {s["reads_total"]} reads but prep wrote {ps["written"]}')
        per_sum = sum(r['reads'] for r in s['per_bobcode'].values())
        if per_sum != s['reads_assigned']:
            fail.append(f'7 split conservation: sum(per_bobcode.reads) {per_sum} != reads_assigned {s["reads_assigned"]}')
        if 'bam_records_total' in s:
            bsum = sum(r.get('bam_records', 0) for r in s['per_bobcode'].values())
            if bsum + s.get('bam_records_unassigned', 0) != s['bam_records_total']:
                fail.append(f'7 split conservation (BAM): per-sample records {bsum} + unassigned {s.get("bam_records_unassigned",0)} != pooled {s["bam_records_total"]}')
        else:
            note('7 BAM conservation not checked: demux stats carry no bam_records_total')
    if dm:
        note('7 split conservation checked')

    # 8. every declared bobcode present with a read floor; undeclared bobcodes below noise
    for s in dm:
        tot = max(1, s['reads_total'])
        # ceil, and never below 1: a declared code with ZERO reads must always fail
        floor = max(args.min_reads_per_sample, 1,
                    -(-int(args.min_frac_per_sample * 1e6 * s.get('reads_assigned', tot)) // 1000000))
        for code, r in s['per_bobcode'].items():
            if r['reads'] < floor:
                fail.append(f'8 declared bobcode {code} ({r.get("sample")}) has {r["reads"]} reads, '
                            f'floor {floor} ({100*args.min_frac_per_sample:.2f}% of assigned): '
                            f'a declared sample that is not in the pool, or a dropped lane')
        und, has_map = _undeclared_map(s)
        if s.get('reads_undeclared', 0) > 0 and not has_map:
            fail.append(f'8 demux stats report {s["reads_undeclared"]} undeclared reads but carry no '
                        f'per-code tally (undeclared / undeclared_bobcodes): noise check impossible')
        for code, n in und.items():
            n = n.get('reads', n) if isinstance(n, dict) else n
            if n / tot > args.max_undeclared_frac:
                fail.append(f'8 UNDECLARED bobcode {code} carries {n} reads ({100*n/tot:.2f}% of the library): the run json is missing a sample')
    if dm:
        note(f'8 bobcode presence checked ({sum(len(_undeclared_map(s)[0]) for s in dm)} undeclared codes tallied)')

    # 9. chemistry consistency: every sidecar's UMI matches what the RT primers imply
    for name, d in list(side.items()) + ([('library', lib_side)] if lib_side else []):
        if d.get('empty_sample'):
            continue
        _check_sidecar(name, d, cfg, fail, note)
    if side or lib_side:
        note('9 chemistry/UMI consistency checked')

    # 10. dedup sanity
    for name, d in side.items():
        n_in = d.get('fastq_reads', d.get('n_reads', 0)); kept = d.get('dedup_kept', 0)
        if kept > n_in:
            fail.append(f'10 {name}: dedup kept {kept} > input {n_in}')
    if side and lib_side and not lib_side.get('empty_sample'):
        ks = sum(d.get('dedup_kept', 0) for d in side.values()); kl = lib_side.get('dedup_kept', 0)
        if kl and abs(ks - kl) / kl > args.dedup_tolerance:
            fail.append(f'10 sum(per-sample kept) {ks} vs library kept {kl}: differ by {100*abs(ks-kl)/kl:.2f}% > {100*args.dedup_tolerance:.0f}%')
        note('10 dedup sanity checked')

    # 11. the published pre-dedup FASTQ carries exactly what prep wrote
    if present('prepped_fastq') and ps and os.path.exists(args.prepped_fastq):
        n = _fastq_reads(args.prepped_fastq)
        if n != ps.get('written'):
            fail.append(f'11 published prepped FASTQ has {n} reads, prep wrote {ps.get("written")}')
        note('11 published input checked')

    # 12. species direction PER BOBCODE, against the sample sheet.
    smap = cfg.get('tso_species_map') or cfg.get('bob_species_map') or {}
    label_code = {v: k for k, v in declared.items()}
    n12 = 0
    for r in comp:
        code = label_code.get(r.get('sample'))
        try:
            pct = float(r.get('pct_correct_mrna_only')); n = int(float(r.get('n_scorable_mrna') or 0))
        except (TypeError, ValueError):
            continue
        exp = r.get('expected_species'); decl = smap.get(code) if code else None
        if decl and exp and decl != exp:
            fail.append(f'12 {r["sample"]} ({code}): composition scored against {exp} but the sheet declares {decl}')
        if n < N12_MIN:
            note(f'12 {r["sample"]} ({code}): only {n} scorable mRNA reads, species direction not callable')
            continue
        n12 += 1
        if pct < 100 * CONF12_FAIL:
            fail.append(f'12 {r["sample"]} ({code}): {pct:.1f}% of {n} scorable mRNA reads align to the declared '
                        f'species {exp}: the sheet direction is INVERTED for this bobcode')
        elif pct < 100 * CONF12_LOW:
            note(f'12 {r["sample"]} ({code}): LOW-CONFIDENCE species direction ({pct:.1f}% of {n})')
    note(f'12 species direction per bobcode: {n12} callable of {len(comp)} samples')
    _finish(fail)


# Invariant-12 thresholds. Same meaning as detect_species_map.N_MIN / CONF_MIN; kept
# in one place here and imported from there when the ONT tree is on the path so the
# two cannot drift apart.
N12_MIN, CONF12_FAIL, CONF12_LOW = 50, 0.50, 0.70
try:
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'run'))
    import detect_species_map as _dsm
    N12_MIN, CONF12_LOW = _dsm.N_MIN, _dsm.CONF_MIN
except Exception:
    pass


def _check_sidecar(name, d, cfg, fail, note):
    """9. the sidecar's UMI source matches what the run's RT primers imply, and the
    dedup key is the same-position definition (D)."""
    rt = ' '.join(str(x) for x in (cfg.get('rt_primers_used') or [])).lower()
    # random priming: the whole degenerate stretch after the bobcode (6N + spacer)
    want = '32N' if '26n' in rt else ('6N+HH' if any(t in rt for t in ('6n', '9n', 'random', 'hexamer')) else None)
    if want and d.get('umi_source') != want:
        fail.append(f'9 {name}: sidecar umi_source={d.get("umi_source")} but rt_primers_used={rt!r} implies {want}')
    if d.get('dedup_key') and 'position' not in d['dedup_key']:
        fail.append(f'9 {name}: dedup_key={d["dedup_key"]} is not the same-position definition')


def _finish(fail):
    if fail:
        sys.stderr.write('\nvalidate_outputs: FAILED\n')
        for f in fail:
            sys.stderr.write(f'  - {f}\n')
        sys.exit(1)
    sys.stderr.write('validate_outputs: all invariants hold\n')


if __name__ == '__main__':
    main()

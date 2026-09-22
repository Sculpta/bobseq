#!/usr/bin/env python3
"""prep_reads_illumina.py — ILLUMINA PE150 read prep for the 3'TSO species-mixing pipeline.

Replaces the ONT step (trim_ont_adapters_illumina.py), which is meaningless here: Illumina
reads carry no ONT LSK adapter / native barcode, and the TruSeq handles are the SEQUENCING
PRIMERS, so they are not present in the read at all.

WHAT THIS DOES (read layout verified on a 2x150 run):

  READ 2 — the informative mate. Sequencing starts immediately after the TruSeq2 primer, so
  the read BEGINS with the species bobcode. There is NO upstream anchor to search for; the
  ONT code found the 7-mer via strand.find('CTCTTCCGATCT')+12, which cannot work here. This
  module calls the bobcode POSITIONALLY at seq[0:7] instead.

      [bobcode 7nt][UMI 6nt][spacer][GGG...][cDNA insert]------->[polyA read-through]
      spacer length = grun_offset - 6   (2 nt "HH" for a +8 code, 5 nt "HHWMW" for a +11 code)
      G-run start   = 7 + grun_offset   (pos 15 / pos 18 for a +8 / +11 pair)

  READ 1 — depends on the RT primer (rt_primers_used in the guide):
    R1_26N_polydT : [26N UMI][polyT ~30nt][degraded sequence]. Quality collapses after the
                    homopolymer, so R1 is never aligned; only R1[0:26] is taken (the UMI).
    R1long-6N     : [priming hexamer][cDNA toward the RNA 5' end]. Real cDNA: written
                    UNTOUCHED to --r1-out with the same read name, and aligned as the mate
                    of R2 (paired-end). On this
                    chemistry the UMI is the WHOLE degenerate stretch after the bobcode up to
                    the G-run (6N + spacer, 8/11/14/7 nt by code set; every position measured
                    degenerate), not the first 6 nt.

OUTPUT: the trimmed R2 FASTQ (mate 1 of the alignment), optionally the R1 mate FASTQ, with the
UMI(s) and the called bobcode carried in the read name so no downstream module has to re-parse
the raw read:

    @<readid>_<UMI1_26nt>_<UMI2> BC:<bobcode> MM:<mismatches> GR:<grun_len>

The umi_tools "_UMI" name convention is used so the read id stays the join key (everything
downstream keys on `rid = qname.split()[0]`, and appending after an underscore keeps that
stable while making the UMIs recoverable).

NOTHING IS HARDCODED: the bobcodes and G-run offsets come from the run's analysisguide
(tso_species_map + grun_offsets). This is deliberate — the ONT modules hardcode 7-mers
(override_24plex BC_HUMAN/BC_MOUSE, filter_funnel _both_bc), which silently makes those
filters inert whenever a run uses a different bobcode pair.

USAGE:
  python3 prep_reads_illumina.py <R1.fq[.gz]> <R2.fq[.gz]> <out_R2.fq.gz> \
      [--guide <analysisguide.txt>] [--bobcodes TGGATGA=8,TGTCCAC=11] \
      [--max-n N] [--min-insert 30] [--stats out.json] [--r1-out <mate_R1.fq.gz>]
  --r1-out is written for every chemistry so the pipeline's file layout is static: on random
  priming it holds the R1 mate of every written R2 (same order, same names); on poly-dT it is
  an EMPTY gzip, which the alignment step treats as "single-end".
"""
import argparse
import gzip
import json
import os
import re
import sys

BOBCODE_LEN = 7
UMI2_LEN = 6          # bobcode-side UMI core, immediately 3' of the bobcode (poly-dT runs:
                      # exactly this; random priming: the whole stretch to the G-run, see main)
UMI1_LEN = 26         # DEFAULT ONLY. Set from the guide's rt_primers_used in main()
                      # via umi_dedup_illumina.configure_umi; 0 for random priming.
_POLYA = re.compile(r'A{10,}')
_POLYT = re.compile(r'T{10,}')      # R1 polydT tract (RT-end evidence)


def _open(p, mode='rt'):
    return gzip.open(p, mode) if p.endswith('.gz') else open(p, mode)


def parse_guide(path):
    """-> ({7mer: species}, {7mer: grun_offset}) from a materialized analysisguide."""
    smap, offs = {}, {}
    infl = False
    for ln in open(path):
        s = ln.strip()
        if s == 'BEGIN_CONFIG':
            infl = True
            continue
        if s == 'END_CONFIG':
            break
        if not infl or '\t' not in ln:
            continue
        k, v = ln.rstrip('\n').split('\t', 1)
        k = k.strip()
        if k in ('tso_species_map', 'bob_species_map'):
            for pair in v.split(','):
                if '=' in pair:
                    a, b = pair.split('=', 1)
                    smap.setdefault(a.strip(), b.strip())
        elif k == 'grun_offsets':
            for pair in v.split(','):
                if '=' in pair:
                    a, b = pair.split('=', 1)
                    try:
                        offs[a.strip()] = int(b.strip())
                    except ValueError:
                        pass
    return smap, offs


def guide_field(path, key):
    """Single scalar field from a materialized analysisguide, or ''."""
    if not path:
        return ''
    infl = False
    for ln in open(path):
        s = ln.strip()
        if s == 'BEGIN_CONFIG':
            infl = True; continue
        if s == 'END_CONFIG':
            break
        if infl and '\t' in ln:
            k, v = ln.rstrip('\n').split('\t', 1)
            if k.strip() == key:
                return v.strip()
    return ''


# --- READ GEOMETRY -> MEASURABILITY -----------------------------------------
# The analysis assumes the reads are long enough to CONTAIN the elements it
# measures. Read length is a property of the RUN, not of the library, and every
# provider picks it independently: one run may be 2x150, another 28/94 with R1
# covering only the 26N UMI. When an element is not in the read, every check for
# it returns 0 and the report would show TOTAL STRUCTURAL FAILURE instead of an
# inapplicable filter.
#
# This CANNOT be inferred from "no reads passed": a genuinely failed library is
# indistinguishable from an unobservable element by that test. So geometry is
# MEASURED on every run, capabilities are DERIVED from it, and anything
# unavailable must be ACKNOWLEDGED in the run JSON as "unmeasurable": [...].
# A normal run needs no config; an unusual one stops until a human agrees, so a
# provider silently shortening a read can never silently change what we report.
HARD_CAPABILITIES = ('umi1',)     # never waivable: silently corrupts, not just blanks


def derive_capabilities(r1_len, r2_len, offsets):
    """{name: (available, detail)} for every geometry-dependent element.

    Each rule states what the READ must contain, so the failure message can say
    what was found and what is needed rather than just refusing.
    """
    room1 = r1_len - UMI1_LEN                       # bases past the 26N UMI
    need_grun = (max(len(c) + o for c, o in offsets.items()) if offsets else BOBCODE_LEN) + 3
    return {
        # A short R1 truncates the UMI itself. Dedup would still run and still
        # produce a number, on a shorter key -- silent collision inflation. Hard.
        'umi1':  (r1_len >= UMI1_LEN,
                  f'R1 is {r1_len} nt, the UMI alone needs {UMI1_LEN}'),
        # The RT-end polydT tract: >=10 T starting right after the UMI.
        'polyt': (room1 >= 10,
                  f'R1 is {r1_len} nt = {room1} past the {UMI1_LEN}N UMI, '
                  f'a >=10 nt polydT tract needs 10'),
        # Non-templated G-run, at 7 + the largest per-bobcode offset.
        'grun':  (r2_len >= need_grun,
                  f'R2 is {r2_len} nt, the G-run sits at {need_grun}'),
    }


def enforce_geometry(st, guide_path, offsets, n):
    """Measure -> derive -> require acknowledgement. Returns the capability map."""
    r1_len = st.get('r1_len_max', 0)
    r2_len = st.get('r2_len_max', 0)
    caps = derive_capabilities(r1_len, r2_len, offsets)
    acked = {x.strip() for x in (guide_field(guide_path, 'unmeasurable') or '').split(',')
             if x.strip()}
    st['read_geometry'] = {'r1_len': r1_len, 'r2_len': r2_len}
    st['capabilities'] = {k: ok for k, (ok, _) in caps.items()}
    st['unmeasurable_acked'] = sorted(acked)
    # RT-end is measured only on bobcode-called reads (the loop `continue`s before it on
    # a failed call), so that is the denominator; over all pairs it would be diluted by
    # the no-bobcode fraction.
    _called = st.get('bobcode_called', 0) or 1
    st['pct_rt_end'] = round(100 * st.get('rt_end', 0) / _called, 3)
    st['pct_rt_end_denominator'] = 'bobcode_called'

    print(f'  read geometry: R1 {r1_len} nt, R2 {r2_len} nt')
    for k, (ok, why) in sorted(caps.items()):
        print(f'    {k:<6} {"available" if ok else "NOT OBSERVABLE"}  ({why})')

    unknown = sorted(acked - set(caps))
    if unknown:
        sys.exit(f'prep_reads: FATAL "unmeasurable" names {unknown}, which are not '
                 f'geometry capabilities. Known: {sorted(caps)}.')
    for k in sorted(k for k, (ok, _) in caps.items() if not ok):
        if k in HARD_CAPABILITIES:
            sys.exit(f'prep_reads: FATAL {k} is not observable at this read length '
                     f'({caps[k][1]}). This cannot be waived with "unmeasurable": a '
                     f'truncated UMI does not blank the metric, it silently shortens '
                     f'the dedup key and inflates collisions.')
        if k not in acked:
            sys.exit(f'prep_reads: FATAL {k} is not observable at this read length '
                     f'({caps[k][1]}), and the run does not acknowledge it. Every '
                     f'check for {k} would return 0, and the report would show a '
                     f'total structural failure rather than an inapplicable filter. '
                     f'If that is expected for this run, add "unmeasurable": '
                     f'["{k}"] to the run JSON.')
    stale = sorted(k for k in acked if caps.get(k, (False,))[0])
    if stale:
        sys.exit(f'prep_reads: FATAL "unmeasurable" lists {stale}, but at this read '
                 f'length {"they are" if len(stale) > 1 else "it is"} observable '
                 f'({"; ".join(caps[k][1] for k in stale)}). Dropping the step would '
                 f'discard a working filter -- remove it from "unmeasurable".')
    return caps


def call_bobcode(seq, codes, max_mm=1):
    """POSITIONAL bobcode call at seq[0:len(code)], each code at ITS OWN length, longest
    codes first (the V1 set has 24 seven-mers and the 11-mer BOB25D), at most max_mm
    substitutions, ties rejected as ambiguous. Returns (bobcode, mismatches) or (None, None).
    Measured on ONT reads of the same bobcode set: the G-run offset counts from the END of
    the code whatever its length, so nothing here may assume 7."""
    for L in sorted({len(b) for b in codes}, reverse=True):
        w = seq[:L]
        if len(w) < L:
            continue
        best, bestd, tie = None, max_mm + 1, False
        for b in codes:
            if len(b) != L:
                continue
            d = 0
            for x, y in zip(w, b):
                if x != y:
                    d += 1
                    if d > bestd:
                        break
            if d < bestd:
                best, bestd, tie = b, d, False
            elif d == bestd and d <= max_mm:
                tie = True
        if best is not None and bestd <= max_mm:
            return (None, None) if tie else (best, bestd)
    return None, None


def grun_count(seq, start):
    """Length of the G-run at the EXACT expected position. Not capped: template
    switching can add well over 3 non-templated G. 0 = no template switch, which is a
    real QC signal rather than a detection miss."""
    n = 0
    while start + n < len(seq) and seq[start + n] == 'G':
        n += 1
    return n


def trim_polya(seq, qual):
    """Cut 3' polyA read-through (measured on a 24-plex library: 29% of reads run off the
    cDNA end into the polydT). Everything from the first >=10A run onward is removed."""
    m = _POLYA.search(seq)
    if m:
        return seq[:m.start()], qual[:m.start()]
    return seq, qual



# ---- per-sample mode: merge the statistics of chunked preps ------------------------
# Prep can run on chunks of the library in parallel; every per-read decision
# is identical (prep has no cross-read state), so the chunk outputs concatenate and the
# statistics ADD. Counts and count-dicts are summed, the read-length maxima are maxed, the
# constants must agree, and the derived fields are recomputed by the same code that
# computes them for a single run (pct_* here, the geometry block by enforce_geometry).
_STAT_COUNTS = ('n', 'bobcode_called', 'no_bobcode', 'mm0', 'mm1', 'grun_missing',
                'polya_trimmed', 'too_short_after_trim', 'written', 'rt_end',
                'polya_tail_ge40_other', 'polya_tail_lt40',
                'polya_tail_readthrough_adapter', 'polya_tail_lowcx', 'mate_r1_written')
_STAT_COUNT_DICTS = ('by_bobcode', 'grun_len')
_STAT_MAXIMA = ('r1_len_max', 'r2_len_max')
_STAT_CONSTANTS = ('species_map', 'grun_offsets', 'umi_source', 'aligned_mode')


def merge_stats(parts, guide_path):
    """One prep_stats dict from the per-chunk dicts, equal (as a dict) to the stats a
    single prep of the concatenated input would have written."""
    if not parts:
        sys.exit('prep_reads merge-stats: no chunk statistics given')
    st = {}
    for k in ('n', 'bobcode_called', 'no_bobcode'):
        st[k] = sum(int(p.get(k, 0)) for p in parts)
    st['by_bobcode'] = {}
    for p in parts:
        for c, v in (p.get('by_bobcode') or {}).items():
            st['by_bobcode'][c] = st['by_bobcode'].get(c, 0) + int(v)
    for k in ('mm0', 'mm1', 'grun_missing', 'polya_trimmed', 'too_short_after_trim', 'written'):
        st[k] = sum(int(p.get(k, 0)) for p in parts)
    gl = {}
    for p in parts:
        for c, v in (p.get('grun_len') or {}).items():
            gl[str(c)] = gl.get(str(c), 0) + int(v)
    st['grun_len'] = gl
    for k in _STAT_CONSTANTS:
        vals = [p.get(k) for p in parts]
        if any(v != vals[0] for v in vals):
            sys.exit(f'prep_reads merge-stats: chunks disagree on {k}: {vals}')
        if vals[0] is not None:
            st[k] = vals[0]
    for k in _STAT_MAXIMA:
        st[k] = max(int(p.get(k, 0)) for p in parts)
    st['rt_end'] = sum(int(p.get('rt_end', 0)) for p in parts)
    for k in ('polya_tail_ge40_other', 'polya_tail_lt40', 'polya_tail_readthrough_adapter',
              'polya_tail_lowcx', 'mate_r1_written'):
        v = sum(int(p.get(k, 0)) for p in parts)
        if v or any(k in p for p in parts):
            st[k] = v
    n = st['n'] or 1
    st['pct_bobcode_called'] = round(100 * st['bobcode_called'] / n, 2)
    st['pct_written'] = round(100 * st['written'] / n, 2)
    enforce_geometry(st, guide_path, st['grun_offsets'], n)
    return st


def merge_main(argv):
    """prep_reads_illumina.py merge-stats --guide G --out merged.json chunk1.json chunk2.json ..."""
    ap = argparse.ArgumentParser(prog='prep_reads_illumina.py merge-stats')
    ap.add_argument('parts', nargs='+')
    ap.add_argument('--guide', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args(argv)
    st = merge_stats([json.load(open(p)) for p in a.parts], a.guide)
    json.dump(st, open(a.out, 'w'), indent=1)
    print(f'prep_reads merge-stats: {len(a.parts)} chunk(s) -> {a.out}: n={st["n"]:,} '
          f'called={st["bobcode_called"]:,} written={st["written"]:,}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('r1'); ap.add_argument('r2'); ap.add_argument('out')
    ap.add_argument('--guide', default=None)
    ap.add_argument('--bobcodes', default=None,
                    help='fallback if no guide: "TGGATGA=8,TGTCCAC=11" (7mer=grun_offset)')
    ap.add_argument('--max-mm', type=int, default=1)
    ap.add_argument('--min-insert', type=int, default=30)
    ap.add_argument('--max-n', type=int, default=0, help='0 = all reads')
    ap.add_argument('--stride', type=int, default=1,
                    help='keep every Nth read pair (1 = all). Use this INSTEAD of '
                         'pre-subsampling with awk/head: prep already streams both '
                         'FASTQs, so striding here costs one pass instead of '
                         'decompress->awk->recompress->decompress. Striding also '
                         'beats head-sampling for QC, because an Illumina FASTQ is in '
                         'flowcell order — the first N reads come from one or two '
                         'tiles and are enriched for optical/ExAmp (spatially '
                         'adjacent) duplicates, which inflates the duplicate rate.')
    ap.add_argument('--stats', default=None)
    ap.add_argument('--r1-out', default=None,
                    help='R1 mate FASTQ (gz) of every written read, same names and order; '
                         'empty on chemistries whose R1 is not cDNA (poly-dT)')
    a = ap.parse_args()
    global UMI1_LEN
    _rt = guide_field(a.guide, 'rt_primers_used') if a.guide else ''
    if not _rt:
        sys.exit('prep_reads: FATAL the guide carries no rt_primers_used; the UMI length '
                 'is derived from it and will not be guessed.')
    import umi_dedup_illumina as _ud
    _umi = _ud.configure_umi(_rt)
    UMI1_LEN = _umi['umi1_len']
    # Random priming: R1 is cDNA -> paired-end alignment, and the TSO UMI is the whole
    # degenerate stretch up to the G-run (length = the code's G-run offset).
    PAIRED = _umi['umi_source'] in ('6N+HH', '6N+spacer')
    print(f"  UMI from chemistry: rt_primers_used={_rt} -> {_umi['umi_source']} "
          f"(R1 UMI {UMI1_LEN} nt, TSO UMI {'written in full to the read name, dedup uses ' + _umi['umi_source'] if PAIRED else f'{UMI2_LEN} nt'}, "
          f"ed<={_umi['umi_max_ed']}); alignment {'paired-end (R1 mate written)' if PAIRED else 'R2 only'}")

    smap, offs = ({}, {})
    if a.guide and os.path.exists(a.guide):
        smap, offs = parse_guide(a.guide)
    if a.bobcodes:
        for pair in a.bobcodes.split(','):
            if '=' in pair:
                k, v = pair.split('=', 1)
                offs[k.strip()] = int(v.strip())
    if not offs:
        sys.exit('prep_reads: no bobcodes/grun_offsets (pass --guide or --bobcodes)')
    codes = list(offs)

    st = {'n': 0, 'bobcode_called': 0, 'no_bobcode': 0, 'by_bobcode': {c: 0 for c in codes},
          'mm0': 0, 'mm1': 0, 'grun_missing': 0, 'polya_trimmed': 0,
          'too_short_after_trim': 0, 'written': 0, 'grun_len': {}, 'species_map': smap,
          'grun_offsets': offs, 'umi_source': _umi['umi_source'],
          'aligned_mode': 'paired' if PAIRED else 'single', 'mate_r1_written': 0}

    r1out = gzip.open(a.r1_out, 'wt') if a.r1_out else None
    with _open(a.r1) as f1, _open(a.r2) as f2, gzip.open(a.out, 'wt') as out:
        while True:
            h1 = f1.readline()
            h2 = f2.readline()
            if not h1 or not h2:
                break
            s1 = f1.readline().rstrip('\n'); f1.readline(); q1 = f1.readline().rstrip('\n')
            s2 = f2.readline().rstrip('\n'); f2.readline(); q2 = f2.readline().rstrip('\n')
            st['n'] += 1
            if a.max_n and st['n'] > a.max_n:
                st['n'] -= 1
                break
            if a.stride > 1 and (st['n'] - 1) % a.stride:
                st['strided_out'] = st.get('strided_out', 0) + 1
                continue
            rid1 = h1[1:].split()[0]
            rid2 = h2[1:].split()[0]
            if rid1 != rid2:
                sys.exit(f'prep_reads: R1/R2 out of sync at record {st["n"]}: {rid1} vs {rid2}')

            bc, mm = call_bobcode(s2, codes, a.max_mm)
            if bc is None:
                st['no_bobcode'] += 1
                continue
            st['bobcode_called'] += 1
            st['by_bobcode'][bc] += 1
            st['mm0' if mm == 0 else 'mm1'] += 1

            gstart = len(bc) + offs[bc]              # code end + (UMI 6 + spacer): measured, not 7
            glen = grun_count(s2, gstart)            # full run at the exact position
            # G-RUN IS NOT TRIMMED. The cut ends just BEFORE the
            # G-run, so the emitted read begins with it and STAR soft-clips it
            # (--alignEndsType Local). Removing the whole run would also remove any
            # TEMPLATED G at a transcript's true 5' end -- non-templated and
            # templated G are indistinguishable by sequence, so a greedy trim would silently
            # eat 1-2 real bases whenever a transcript started with G. Leaving the run in
            # lets the aligner decide where the genome match actually begins.
            # glen is MEASURED here (exact expected offset) and carried in the
            # QNAME, so every downstream G-run metric is unchanged -- it is deliberately
            # not re-derived from the soft-clip, which would conflate non-templated G,
            # templated G, and ordinary end-quality clipping.
            gend = gstart
            if glen == 0:
                st['grun_missing'] += 1
            st['grun_len'][glen] = st['grun_len'].get(glen, 0) + 1

            # poly-dT: the 6N core (the 26N on R1 carries the rest of the 32N composite).
            # random priming: everything up to the G-run, 8/11/14/7 nt by code (measured
            # degenerate at every position).
            umi2 = s2[len(bc):len(bc) + (offs[bc] if PAIRED else UMI2_LEN)]
            umi1 = s1[:UMI1_LEN]
            # RT-end evidence is on R1 ONLY: a polyT tract starting right after the 26N
            # UMI. Measured here (not downstream) because R1 is not aligned and so is
            # unavailable from the BAM. On R2 an A/T run means 3' polyA READ-THROUGH —
            # i.e. a short insert — so testing R2 would invert the metric's meaning.
            if len(s1) > st.get('r1_len_max', 0):
                st['r1_len_max'] = len(s1)
            if len(s2) > st.get('r2_len_max', 0):
                st['r2_len_max'] = len(s2)
            _pt = _POLYT.search(s1)
            rt_end = bool(_pt and 18 <= _pt.start() <= 40)
            st['rt_end'] = st.get('rt_end', 0) + (1 if rt_end else 0)

            ins, insq = s2[gend:], q2[gend:]
            n0 = len(ins)
            ins, insq = trim_polya(ins, insq)
            if len(ins) < n0:
                st['polya_trimmed'] += 1
                # What sat after the A-run. The trim fires on the FIRST
                # A{10,} anywhere in the insert; a genomic A-tract inside a real insert
                # would be cut too. Measured on two libraries (100k reads each):
                # tails >=40 nt beyond the run are 6% / 37% of reads, and of those the
                # Read-1 primer rc (read-through) appears at ~26 nt (the 26N) in 13% / 7%;
                # the rest are mostly further A-runs and low-complexity sequence, not
                # alignable insert. So the rule stands; this tally makes the potential
                # internal-tract fraction visible per run instead of assumed.
                _tail = s2[gend + len(ins):]
                _m = _POLYA.match(_tail)
                _after = _tail[_m.end():] if _m else _tail
                if len(_after) < 40:
                    _k = 'polya_tail_lt40'
                elif 'AGATCGGAAGAGC' in _after:
                    _k = 'polya_tail_readthrough_adapter'
                elif max(_after.count(b) for b in 'ACGT') / len(_after) > 0.6:
                    _k = 'polya_tail_lowcx'
                else:
                    _k = 'polya_tail_ge40_other'
                st[_k] = st.get(_k, 0) + 1
            # Gate on the TRUE insert (G-run excluded). The emitted read carries the
            # G-run, so testing len(ins) directly would let a read through on the
            # strength of its non-templated G and silently change which reads survive.
            if len(ins) - glen < a.min_insert:
                st['too_short_after_trim'] += 1
                continue

            # STRUCTURE LIVES IN THE READ NAME — this is the single source of truth.
            # Because the bobcode/UMI/spacer/G-run are TRIMMED OFF before alignment, the
            # BAM SEQ no longer contains them; any module that re-parsed structure from
            # SEQ would silently fail. Encoding them here means the FASTQ path and the
            # BAM path read the SAME values, and STAR preserves QNAME (it drops the
            # comment field, which is why none of this may live after a space).
            #   <origid>_<26N R1 UMI>_<6nt R2 UMI>_<bobcode>_<Grun len>_<T|F RT-end>
            # '_' is safe: Illumina read ids use ':' only.
            rt = 'T' if rt_end else 'F'
            name = f'{rid2}_{umi1}_{umi2}_{bc}_{glen}_{rt}'
            out.write(f'@{name}\n{ins}\n+\n{insq}\n')
            st['written'] += 1
            if PAIRED and r1out is not None:
                # the mate, untouched: the priming hexamer stays (STAR soft-clips what does
                # not match), no poly(A) trim (R1 reads away from the 3' end)
                r1out.write(f'@{name}\n{s1}\n+\n{q1}\n')
                st['mate_r1_written'] += 1
    if r1out is not None:
        r1out.close()

    n = st['n'] or 1
    st['pct_bobcode_called'] = round(100 * st['bobcode_called'] / n, 2)
    st['pct_written'] = round(100 * st['written'] / n, 2)

    # Measure -> derive -> require acknowledgement. See derive_capabilities().
    enforce_geometry(st, a.guide, offs, n)
    print(json.dumps(st, indent=1))
    if a.stats:
        json.dump(st, open(a.stats, 'w'), indent=1)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'merge-stats':
        merge_main(sys.argv[2:])
    else:
        main()

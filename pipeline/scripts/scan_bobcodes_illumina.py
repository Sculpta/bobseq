#!/usr/bin/env python3
"""Efficient, config-free auto-detection of the bobcodes used in an experiment,
by scanning read motifs. Answers three questions from the reads alone:

  1. Which chemistry/architecture?  (24plex vs tso-bob vs BOB-J) -- by which
     anchor signature dominates the sampled reads.
  2. Which bobcodes are present?    -- extract the barcode window downstream of
     the winning anchor and match it (exact or <=1 mismatch) to the reference
     registry (bobcode_reference), which carries the NAME / ID / species.
  3. In what proportion?            -- read fraction per detected bobcode.

Only a sample of reads is scanned (default 30k) and detection is str.find +
O(1) dict lookups over precomputed 1-mismatch variants, so it is cheap enough to
run inline during analyze_barcode. Returns a JSON-able dict; see scan().

CLI:  python3 scan_bobcodes_illumina.py <fastq[.gz]> [n_sample]
"""
import os, gzip
from collections import Counter
import bobcode_reference_illumina as ref

_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def rc(s):
    return s.translate(_COMP)[::-1]

# a detected anchor-signature maps to the set of registry architectures whose
# barcodes could follow it (24plex TruSeq tail is shared by A/B/C and D sets).
_ARCH_GROUP = {
    '24plex': ('24plex', '24plex-D'),
    'tso-bob (SMARTer)': ('tso-bob (SMARTer)',),
    'BOB-J (Nextera-R1)': ('BOB-J (Nextera-R1)',),
}


def _iter_fastq(path, n):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as fh:
        i = 0
        for k, line in enumerate(fh):
            if k & 3 == 1:              # every 4th line, offset 1 = sequence
                yield line.strip().upper()
                i += 1
                if i >= n:
                    return


def _lookup_for(arch_group, max_mm=1, extra=None):
    """{variant: (canonical, record)} and the sorted candidate lengths, built
    from the registry archs relevant to the detected chemistry PLUS any per-run
    `extra` bobcodes ({seq: record}). The extras are the run's own configured
    barcodes from its guide -- included so a bobcode that is genuinely in use but
    not yet in the reference sheets (e.g. GTTGCTT) is still detected and named."""
    ref.load()
    subset = {s: r for s, r in ref.BARCODES.items() if r['arch'] in arch_group}
    if extra:
        # a reference-sheet entry keeps precedence for its name; only add extras
        # that add a new sequence.
        for s, r in extra.items():
            subset.setdefault(s, r)
    exact = {s: (s, subset[s]) for s in subset}
    lut = dict(exact)
    if max_mm >= 1:
        for s in subset:
            for i in range(len(s)):
                for b in 'ACGT':
                    if b == s[i]:
                        continue
                    v = s[:i] + b + s[i + 1:]
                    if v not in exact:
                        lut.setdefault(v, (s, subset[s]))
    lengths = sorted({len(s) for s in subset}, reverse=True)   # longest first
    return lut, lengths


def _match_after(strand, start, lut, lengths):
    """Return (canonical_barcode, record) for the barcode window at `start`, or
    None. Longest candidate length first so a D-set 11-mer is not truncated to a
    spurious 7-mer hit."""
    for L in lengths:
        seg = strand[start:start + L]
        if len(seg) < L:
            continue
        hit = lut.get(seg)
        if hit:
            return hit
    return None


def scan(reads, n=30000, min_frac=0.002, min_abs=20, extra=None):
    """reads: a list[str] of sequences OR a path to a FASTQ[.gz].
    extra: optional {seq: record} of the run's own configured bobcodes (from its
    guide), merged into the registry so barcodes not yet in the reference sheets
    are still detected and named.
    Returns a dict describing the detected architecture + bobcodes, or None.

    Architecture is chosen by which anchor is most often *followed by a known
    bobcode* -- not by bare anchor presence -- so an incidental library adapter
    (e.g. a Nextera Read1 primer inside a 24plex library) cannot win, because
    nothing downstream of it matches the registry.

    Detection fires when the winning architecture has >= max(min_abs,
    min_frac*n_scanned) barcode-followed hits. The absolute floor keeps
    low-depth barcodes callable: a few dozen unambiguous <=1-mismatch hits
    downstream of the correct anchor are strong evidence, even at <1% of reads."""
    ref.load()
    sample = list(_iter_fastq(reads, n)) if isinstance(reads, str) else \
        [s.upper() for s in reads[:n]]
    ns = len(sample)
    if ns == 0:
        return None

    archs = [a for a, _ in ref.ARCH_SIGNATURES]
    motifs = {a: (m, rc(m)) for a, m in ref.ARCH_SIGNATURES}
    luts = {a: _lookup_for(_ARCH_GROUP[a], extra=extra) for a in archs}
    anchor_hits = {a: 0 for a in archs}          # anchor present (any downstream)
    matched_hits = {a: 0 for a in archs}         # anchor + known bobcode
    counts = {a: Counter() for a in archs}
    rec_by_seq = {a: {} for a in archs}

    for r in sample:
        rcr = None                               # lazily computed reverse strand
        for a in archs:
            motif, rcm = motifs[a]
            p = r.find(motif)
            strand = r
            if p < 0:
                if rcm in r:
                    if rcr is None:
                        rcr = rc(r)
                    strand = rcr
                    p = strand.find(motif)
                else:
                    continue
            if p < 0:
                continue
            anchor_hits[a] += 1
            lut, lengths = luts[a]
            hit = _match_after(strand, p + len(motif), lut, lengths)
            if hit:
                canon, rec = hit
                matched_hits[a] += 1
                counts[a][canon] += 1
                rec_by_seq[a][canon] = rec

    # winner = most anchor+bobcode matches; ties fall back to bare anchor count,
    # then to ARCH_SIGNATURES order (TruSeq before Nextera-only).
    best_arch = max(archs, key=lambda a: (matched_hits[a], anchor_hits[a],
                                          -archs.index(a)))
    matched = matched_hits[best_arch]
    if matched < max(min_abs, int(min_frac * ns)):
        return dict(architecture=None, n_scanned=ns, arch_counts=anchor_hits,
                    matched_by_arch=matched_hits, barcodes=[], n_matched=0,
                    frac_matched=0.0, anchor_signature=None)

    bars = []
    for seq, c in counts[best_arch].most_common():
        rec = rec_by_seq[best_arch][seq]
        bars.append(dict(seq=seq, name=rec['name'], id=rec['id'],
                         arch=rec['arch'], species=rec['species'],
                         count=c, frac=round(c / ns, 4)))
    return dict(architecture=best_arch, n_scanned=ns, arch_counts=anchor_hits,
                matched_by_arch=matched_hits,
                anchor_signature=motifs[best_arch][0], n_matched=matched,
                frac_matched=round(matched / ns, 4), barcodes=bars)


if __name__ == '__main__':
    import sys, json
    if len(sys.argv) < 2:
        sys.exit("usage: scan_bobcodes_illumina.py <fastq[.gz]> [n_sample]")
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 30000
    res = scan(sys.argv[1], n=n)
    if not res:
        sys.exit("no reads / empty")
    print(f"architecture : {res['architecture']}   "
          f"(anchor {res['anchor_signature']})")
    print(f"scanned      : {res['n_scanned']:,} reads   "
          f"matched {res['n_matched']:,} ({res['frac_matched']*100:.1f}%)")
    print(f"arch signal  : " + "  ".join(
        f"{a}: anchor={res['arch_counts'][a]}/matched={res['matched_by_arch'][a]}"
        for a in res['arch_counts']))
    print("detected bobcodes:")
    for b in res['barcodes'][:8]:
        sp = f" [{b['species']}]" if b['species'] else ""
        print(f"  {b['seq']:12s} {b['id'] or '':8s} {b['name']:16s}"
              f"{sp:9s}  {b['count']:>7,} ({b['frac']*100:.1f}%)")

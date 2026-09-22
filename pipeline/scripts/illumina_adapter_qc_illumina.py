#!/usr/bin/env python3
"""Illumina P5/P7 adapter QC for ONT reads — fast seed-and-verify matching.

Adapter matching is done with a pigeonhole seed + Hamming verify instead of a
whole-read fuzzy regex:

  * To match a pattern with <= m mismatches, split it into m+1 non-overlapping
    seeds — by the pigeonhole principle at least one seed must match exactly.
    `str.find` (C-speed) locates candidate seeds; the O(len) Hamming check runs
    ONLY at those few candidate offsets (not at every position like fuzzy regex).
  * `rc(read)` is computed once per read and reused.

Each candidate offset is verified by a Hamming fast path (substitutions only); if
that fails it falls back to a small banded edit-distance DP so ONT indels inside
the adapter are still caught — the DP therefore runs on only the handful of
indel-bearing / near-miss offsets, not on every read. The DP also returns the
indel-aware END of the match, so the i5/i7 index window is extracted at the right
offset even when the primer carries an insertion or deletion. The results are
equivalent to a whole-read fuzzy-regex search at a fraction of the cost.
"""
import gzip
import os
import collections

COMP = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(COMP)[::-1]

# --- primer / handle constants -----------------------------------------------
P5 = 'AATGATACGGCGACCACCGAGAT'
P7 = 'CAAGCAGAAGACGGCATACGAGAT'
P5_HALVES = ('AATGATACGGCG', 'CGACCACCGAGAT')
P7_HALVES = ('CAAGCAGAAGAC', 'GGCATACGAGAT')
R1 = 'CTACACGACGCTCTTCCGATCT'
R2 = 'AGACGTGTGCTCTTCCGATCT'
I5_OFFSET = 6
I7_OFFSET = 0
IDX_WIN = 12
# Index CORE reporting. The 12-nt window above starts at a fixed offset after
# the P5/P7 handle, so it carries linker bases and clips the index (67C primers: i5 window =
# index + 2 nt of the Read1 primer; i7 window = 5-nt linker GTAGA + 7 nt of the index),
# which is not what a sample sheet lists. The core is the 10 nt immediately upstream of
# the Read1 / Read2 sequencing primer, whatever the linker length. i5 is reported forward;
# i7 is reported as the reverse complement of the oligo (the "index 1" a sequencer sheet
# expects), matching the primer-table convention. The 12-nt windows stay in the JSON for continuity.
READ1_PRIMER = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'
READ2_PRIMER = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
CORE_LEN = 10
CORE_SCAN = 40            # nt after the handle end within which the read primer must start

# --- NEBNext dual-index identification --------------------------------------
# Match the detected i5/i7 index window against the saved NEBNext index list so
# the report can name the primer (e.g. i501/i703) or flag 'unknown' (non-NEB /
# TruSeq). Single source of truth = the dual-index primer sheet below, located via the
# BOBSEQ_NEB_INDEX_FILE env var with the config/ copy as fallback. The location matters
# more than it looks: _load_neb_indexes() swallows OSError and returns an EMPTY registry,
# after which match_neb_index() always returns None -- which the code documents as meaning
# "no NEB index present (i.e. a non-NEB/TruSeq primer)". A missing file would therefore be
# indistinguishable from a real biological finding, and every sample would be silently
# reported as non-NEB.
_NEB_INDEX_FILE = os.environ.get(
    'BOBSEQ_NEB_INDEX_FILE',
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'config', 'nebnext_e7600_dual_index_primers.tsv'))
_NEB_INDEXES = None

def _load_neb_indexes():
    """{'i5': {id: {seq,...}}, 'i7': {...}} — the in-primer, forward-read and
    revcomp-read forms plus their reverse complements (covers both ONT orientations
    and both sequencer workflows)."""
    global _NEB_INDEXES
    if _NEB_INDEXES is not None:
        return _NEB_INDEXES
    out = {'i5': {}, 'i7': {}}
    try:
        for ln in open(_NEB_INDEX_FILE):
            if ln.startswith('#') or ln.startswith('index_name') or not ln.strip():
                continue
            f = ln.rstrip('\n').split('\t')
            name, kind = f[0], f[2]
            seqs = {s for s in f[3:6] if s}
            seqs |= {rc(s) for s in list(seqs)}
            if kind in out:
                out[kind][name] = seqs
    except OSError:
        pass
    _NEB_INDEXES = out
    return out

def match_neb_index(window, kind):
    """Best NEBNext index id whose 8-mer appears in the detected `window` (i5/i7
    dominant sequence) within <=1 mismatch, slid across the window. Returns the id
    (e.g. 'i501') or None when no NEB index is present (i.e. a non-NEB/TruSeq primer)."""
    if not window:
        return None
    for name, seqs in _load_neb_indexes().get(kind, {}).items():
        for idx in seqs:
            L = len(idx)
            if L > len(window):
                continue
            for j in range(len(window) - L + 1):
                if sum(a != b for a, b in zip(window[j:j + L], idx)) <= 1:
                    return name
    return None
JUNC_WIN = 60


def _seeds(pat, mm):
    """(offset, seed) pigeonhole pieces of `pat`: for any Hamming<=mm occurrence,
    at least one seed matches exactly."""
    parts = mm + 1
    n = len(pat)
    sz = n // parts
    return tuple((i * sz, pat[i * sz: ((i + 1) * sz if i < parts - 1 else n)])
                 for i in range(parts))


def _ham_end(text, c, pat, mm):
    """Fast path: substitution-only. Returns c+len(pat) if Hamming<=mm else -1."""
    m = len(pat)
    if c < 0 or c + m > len(text):
        return -1
    d = 0
    tc = text
    for i in range(m):
        if tc[c + i] != pat[i]:
            d += 1
            if d > mm:
                return -1
    return c + m


def _edit_end(text, cand, pat, mm):
    """Indel-aware fallback: banded semiglobal DP over a small window around `cand`.
    Returns the END position (in `text`) of the best <=mm edit-distance match, else
    -1. Runs only when the Hamming fast path fails (indel-bearing or non-matches)."""
    m = len(pat)
    a = cand - mm
    if a < 0:
        a = 0
    seg = text[a: cand + m + mm]
    ns = len(seg)
    if ns < m - mm:
        return -1
    prev = [0] * (ns + 1)                       # free start gap
    for i in range(1, m + 1):
        pci = pat[i - 1]
        cur = [i]
        cjp = i
        rmin = i
        pj = prev
        for j in range(1, ns + 1):
            v = pj[j - 1] if seg[j - 1] == pci else pj[j - 1] + 1
            d = pj[j] + 1
            if d < v:
                v = d
            e2 = cjp + 1
            if e2 < v:
                v = e2
            cur.append(v)
            cjp = v
            if v < rmin:
                rmin = v
        if rmin > mm:
            return -1
        prev = cur
    bj = 0; bv = prev[0]                         # free end gap: best-scoring end column
    for j in range(1, ns + 1):
        if prev[j] < bv:
            bv = prev[j]; bj = j
    return (a + bj) if bv <= mm else -1


def _match_end(text, cand, pat, mm):
    """End position of a <=mm match near `cand` (Hamming fast path, DP fallback), or -1."""
    e = _ham_end(text, cand, pat, mm)
    return e if e >= 0 else _edit_end(text, cand, pat, mm)


def _first_end(text, pat, seeds, mm):
    """END position of the leftmost <=mm match of `pat` in `text`, else -1
    (indel-aware — used by index/contiguity extraction so the offset is correct)."""
    n = len(text); m = len(pat); best_c = -1; best_e = -1; find = text.find
    for off, sd in seeds:
        start = 0
        while True:
            p = find(sd, start)
            if p < 0:
                break
            c = p - off
            if -mm <= c <= n - m + mm and (best_c < 0 or c < best_c):
                cc = c if c > 0 else 0
                e = _match_end(text, cc, pat, mm)
                if e >= 0:
                    best_c = cc; best_e = e
            start = p + 1
    return best_e


def _first(text, pat, seeds, mm):
    """START position of the leftmost <=mm match of `pat` in `text`, else -1
    (used for the R1-vs-R2 nearest-handle comparison, which needs start offsets)."""
    n = len(text); m = len(pat); best = -1; find = text.find
    for off, sd in seeds:
        start = 0
        while True:
            p = find(sd, start)
            if p < 0:
                break
            c = p - off
            if -mm <= c <= n - m + mm and (best < 0 or c < best):
                cc = c if c > 0 else 0
                if _match_end(text, cc, pat, mm) >= 0:
                    best = cc
            start = p + 1
    return best


def _all(text, pat, seeds, mm):
    """Set of approx start positions of a <=mm match of `pat` in `text`."""
    n = len(text); m = len(pat); out = set(); find = text.find
    for off, sd in seeds:
        start = 0
        while True:
            p = find(sd, start)
            if p < 0:
                break
            c = p - off
            cc = c if c > 0 else 0
            if -mm <= c <= n - m + mm and cc not in out and _match_end(text, cc, pat, mm) >= 0:
                out.add(cc)
            start = p + 1
    return out


def _spec(full, halves, kfull=3, khalf=1):
    """[(pattern, seeds, mm, orient)] = full adapter (both strands) + each half."""
    out = [(full, _seeds(full, kfull), kfull, 'F'),
           (rc(full), _seeds(rc(full), kfull), kfull, 'R')]
    for h in halves:
        out.append((h, _seeds(h, khalf), khalf, 'F'))
        out.append((rc(h), _seeds(rc(h), khalf), khalf, 'R'))
    return out


_P5spec = _spec(P5, P5_HALVES)
_P7spec = _spec(P7, P7_HALVES)
_P5seeds = _seeds(P5, 3); _P7seeds = _seeds(P7, 3)
_R1seeds = _seeds(R1, 3); _R2seeds = _seeds(R2, 3)
_RD1seeds = _seeds(READ1_PRIMER, 3); _RD2seeds = _seeds(READ2_PRIMER, 3)
# TruSeq Read1 / Read2 primer handles as full-length (both-strand) specs, for the
# "R1/R2 on opposite ends" architecture rows. R1 and R2 share the 3' CTCTTCCGATCT
# motif but differ by >3 nt in their 5' halves, so mm=3 does not cross-match them.
_R1spec = _spec(R1, ())
_R2spec = _spec(R2, ())


def _hits(s, spec):
    out = []
    for pat, seeds, mm, ori in spec:
        for c in _all(s, pat, seeds, mm):
            out.append((c, ori))
    return out


def _index_at(s, rcs, primer, seeds, offset):
    """i5/i7 window after the first forward-primer match (read then rc). Uses the
    indel-aware primer END so the index offset stays correct on reads with an
    insertion/deletion inside the primer."""
    for strand in (s, rcs):
        e = _first_end(strand, primer, seeds, 3)
        if e >= 0:
            idx = strand[e + offset: e + offset + IDX_WIN]
            if len(idx) == IDX_WIN:
                return idx
    return None


def _index_core(s, rcs, primer, seeds, read_primer, rp_seeds):
    """10-nt index core = the bases immediately upstream of the Read1/Read2 sequencing
    primer, located after the first P5/P7 handle match (read, then rc). Independent of
    the linker length between handle and index. None when the read primer is not found
    within CORE_SCAN nt or fewer than CORE_LEN nt separate handle and primer."""
    for strand in (s, rcs):
        e = _first_end(strand, primer, seeds, 3)
        if e >= 0:
            w = strand[e: e + CORE_SCAN + len(read_primer)]
            p = _first(w, read_primer, rp_seeds, 3)
            if p >= CORE_LEN:
                return w[p - CORE_LEN: p]
            return None
    return None


def _contig(s, rcs, primer, seeds):
    """Nearer TruSeq handle just interior of the primer: 'R1' / 'R2' / None."""
    for strand in (s, rcs):
        e = _first_end(strand, primer, seeds, 3)
        if e >= 0:
            w = strand[e: e + JUNC_WIN]
            a = _first(w, R1, _R1seeds, 3)
            b = _first(w, R2, _R2seeds, 3)
            pa = a if a >= 0 else 10 ** 9
            pb = b if b >= 0 else 10 ** 9
            if pa == pb == 10 ** 9:
                return None
            return 'R1' if pa <= pb else 'R2'
    return None


def _dominant(idxs):
    if not idxs:
        return (None, None, 0)
    dom, _ = collections.Counter(idxs).most_common(1)[0]
    ham = lambda a, b: sum(x != y for x, y in zip(a, b))
    match = sum(1 for i in idxs if ham(i, dom) <= 1)
    return (dom, round(100 * match / len(idxs), 1), len(idxs))


def _both_ends_same(hits, L):
    return any(((x < 0.4 * L and y > 0.6 * L) or (y < 0.4 * L and x > 0.6 * L)) and o1 != o2
               for x, o1 in hits for y, o2 in hits)


def _opp_ends(hA, hB, L):
    """True if an A-feature and a B-feature sit at opposite read termini in opposite
    orientation (the correct dual-ended architecture: one motif near the 5' end, the
    other near the 3' end, pointing inward)."""
    return any(((x < 0.4 * L and y > 0.6 * L) or (y < 0.4 * L and x > 0.6 * L)) and o1 != o2
               for x, o1 in hA for y, o2 in hB)


def compute_illumina_qc(fastq_path):
    n = p5 = p7 = opp = symm = 0
    tr_opp = tr_symm = 0                               # TruSeq Read1/Read2 handles
    i5s, i7s = [], []
    i5cores, i7cores = [], []                         # 10-nt index cores (oligo orientation)
    p5r1 = p5r2 = p5j = p7r2 = p7r1 = p7j = 0
    reads = []                                        # cached for the index pass 2
    op = gzip.open if fastq_path.endswith('.gz') else open
    with op(fastq_path, 'rt') as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().strip(); fh.readline(); fh.readline()
            n += 1; L = len(s)
            reads.append(s)
            rcs = rc(s)                                   # once per read
            h5 = _hits(s, _P5spec); h7 = _hits(s, _P7spec)
            a, b = bool(h5), bool(h7)
            p5 += a; p7 += b
            if a and b and _opp_ends(h5, h7, L):
                opp += 1
            if _both_ends_same(h5, L) or _both_ends_same(h7, L):
                symm += 1
            # TruSeq Read1 / Read2 primer handles (inner ends of the P5 / P7 sides)
            r1 = _hits(s, _R1spec); r2 = _hits(s, _R2spec)
            if r1 and r2 and _opp_ends(r1, r2, L):
                tr_opp += 1
            if _both_ends_same(r1, L) or _both_ends_same(r2, L):
                tr_symm += 1
            i5 = _index_at(s, rcs, P5, _P5seeds, I5_OFFSET)
            i7 = _index_at(s, rcs, P7, _P7seeds, I7_OFFSET)
            if i5:
                i5s.append(i5)
            if i7:
                i7s.append(i7)
            c5 = _index_core(s, rcs, P5, _P5seeds, READ1_PRIMER, _RD1seeds)
            c7 = _index_core(s, rcs, P7, _P7seeds, READ2_PRIMER, _RD2seeds)
            if c5:
                i5cores.append(c5)
            if c7:
                i7cores.append(c7)
            c5 = _contig(s, rcs, P5, _P5seeds)
            if c5 == 'R1':
                p5r1 += 1; p5j += 1
            elif c5 == 'R2':
                p5r2 += 1; p5j += 1
            c7 = _contig(s, rcs, P7, _P7seeds)
            if c7 == 'R2':
                p7r2 += 1; p7j += 1
            elif c7 == 'R1':
                p7r1 += 1; p7j += 1
    i5win, i5wm, i5wn = _dominant(i5s)                # legacy 12-nt windows
    i7win, i7wm, i7wn = _dominant(i7s)
    i5dom, i5m, i5n = _dominant(i5cores)              # 10-nt cores, oligo orientation
    i7dom_oligo, i7m, i7n = _dominant(i7cores)
    i7dom = rc(i7dom_oligo) if i7dom_oligo else None  # reported as the sheet lists it
    if i5dom is None:                                 # read primer never found: fall back
        i5dom, i5m, i5n = i5win, i5wm, i5wn
    if i7dom is None:
        i7dom, i7m, i7n = i7win, i7wm, i7wn
    # ---- pass 2: the actual i5 / i7 INDEX sequences at opposite ends ------------
    # Stricter than the P5/P7 rows above (which only need the adapter backbone): here
    # the dominant i5 and i7 index windows themselves must be found at opposite read
    # termini (i5i7 = correct dual-index; i5i5/i7i7 = same index both ends, incorrect).
    idx_opp = idx_symm = 0
    idx_scored = 0
    if i5dom and i7dom and i5n >= 20 and i7n >= 20:
        i5spec = _spec(i5win or i5dom, (), kfull=1)   # 12-nt window, allow 1 mismatch
        i7spec = _spec(i7win or i7dom_oligo or i7dom, (), kfull=1)
        for s in reads:
            L = len(s)
            hi5 = _hits(s, i5spec); hi7 = _hits(s, i7spec)
            if not (hi5 or hi7):
                continue
            idx_scored += 1
            if hi5 and hi7 and _opp_ends(hi5, hi7, L):
                idx_opp += 1
            if _both_ends_same(hi5, L) or _both_ends_same(hi7, L):
                idx_symm += 1
        idx_avail = True
    else:
        idx_avail = False
    pct = lambda x, d: round(100 * x / d, 1) if d else None
    return {
        'n': n,
        'partial_p5_pct': pct(p5, n),
        'partial_p7_pct': pct(p7, n),
        'p5p7_opposite_ends_pct': pct(opp, n),
        'p5p5_or_p7p7_pct': pct(symm, n),
        'truseq_r1r2_opposite_ends_pct': pct(tr_opp, n),
        'truseq_r1r1_or_r2r2_pct': pct(tr_symm, n),
        'idx_i5i7_opposite_ends_pct': (pct(idx_opp, n) if idx_avail else None),
        'idx_i5i5_or_i7i7_pct': (pct(idx_symm, n) if idx_avail else None),
        'idx_scored_n': (idx_scored if idx_avail else None),
        'i5_dominant': i5dom, 'i5_match_pct': i5m, 'i5_n': i5n,
        'i7_dominant': i7dom, 'i7_match_pct': i7m, 'i7_n': i7n,
        'i7_dominant_oligo': i7dom_oligo, 'i5_window12': i5win, 'i7_window12': i7win,
        'index_convention': '10-nt core upstream of the Read1/Read2 primer; i5 forward, i7 reverse complement (sample-sheet form)',
        'p5_read1_pct': pct(p5r1, p5j), 'p5_read2_pct': pct(p5r2, p5j), 'p5_junc_n': p5j,
        'p7_read2_pct': pct(p7r2, p7j), 'p7_read1_pct': pct(p7r1, p7j), 'p7_junc_n': p7j,
    }


if __name__ == '__main__':
    import sys, json
    print(json.dumps(compute_illumina_qc(sys.argv[1]), indent=2))

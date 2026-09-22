#!/usr/bin/env python3
"""Gene + UMI duplication (ONT-tolerant, position-ignored) for UMI chemistries.

Collapses protein-coding (non-rRNA, non-MT) mRNA reads by (max-overlap gene, 6 nt
TSO UMI), with the UMI clustered at edit distance <=1 (ONT-error tolerant). Start/
end position is intentionally ignored: ONT reads of the SAME original molecule
truncate at variable points, so positional dedup wrongly splits true duplicates.

Returns {'n_mrna', 'unique', 'dup_reads', 'dup_pct'} or None if no callable UMIs
are found (i.e. the chemistry has no TSO-embedded UMI).
"""
import os
import subprocess
import re
import bisect
import collections
import math
import star_taxonomy_illumina as st
import override_24plex_illumina as ox

_CIG = re.compile(r'(\d+)([MIDN=X])')
def _positional():
    """True when the run is Illumina R2-positional (structure lives in the QNAME)."""
    return getattr(ox, 'READ_LAYOUT', 'ont_anchor') == 'r2_positional'


# --- WHICH UMI (Illumina) ---------------------------------------------------
# Options: '6N' (R2 bobcode UMI only) | '12N' (default) | '26N' (R1 only) | '32N' (both).
#
# WHY NOT 6N ALONE: only 4^6 = 4096 states, so it saturates by birthday collision as a
# dedup group (same gene + same 5' position) grows, INFLATING the duplicate rate.
# Measured against the independent 26N on a 24-plex library: of 12,001 6N-called duplicate
# groups, 86.8% had an identical 26N and 3.7% were within 1 mismatch, but 9.5% had a
# completely different 26N — i.e. ~9.5% of 6N duplicate calls were collisions, not PCR.
# Those false positives concentrated on low-complexity 6-mers (e.g. CGGGGG) and rDNA
# scaffolds, exactly where collisions are expected.
#
#   '12N' : 6 nt R2 bobcode UMI + the FIRST 6 nt of the R1 26N  <-- current choice
# Why the FIRST 6 of R1: the 26N runs from R1 pos 0-25 and is immediately followed by
# the polyT tract, so the LAST bases sit next to the homopolymer and pick up its bias
# (measured T content rises 37% -> 41% across the window). Positions 0-5 are the
# furthest from it and were measured at Q40 with ~1.94 bits/base entropy.
# Space: 4^12 = 16.7M vs 4^6 = 4096, so collision saturation stops being a practical
# concern at any group size this assay reaches.
UMI_SOURCE = '32N'
UMI1_LEN = 26                       # nt of R1 that are UMI; 0 when R1 carries none


def configure_umi(rt_primers):
    """Set UMI_SOURCE / UMI1_LEN / UMI_MAX_ED from the RT primers the run declares.

    Chemistry decides where the UMI is, and nothing else may. Two constructs exist:
      R1_26N_polydT   R1 = [26N UMI][polydT]...; TSO adds 6N   -> 32N composite, ed<=2
      R1long-6N / 9N  random priming. The 6N/9N on R1 is the PRIMING HEXAMER, i.e.
                      templated cDNA, NOT a UMI (the random hexamer is just the six
                      nucleotides and then a primer). The UMI is the TSO's whole degenerate
                      stretch between the bobcode and the G-run: 6N plus the spacer, which
                      is degenerate too (HH / HHWMW / HHMWMWMW; measured per position on
                      24-plex reads). Its length is
                      the code's G-run offset (8 / 11 / 14 / 7 nt) -> '6N+spacer', ed 0.
                      R1 is cDNA on this chemistry and is aligned as the mate of R2.
    Anything unrecognised stops the run: a hardcoded 26 would silently deduplicate a
    random-priming library on 20 nt of cDNA plus 6 nt of UMI.
    """
    global UMI_SOURCE, UMI1_LEN, UMI_MAX_ED
    names = [str(x) for x in (rt_primers if isinstance(rt_primers, (list, tuple))
                              else str(rt_primers or '').split(','))]
    joined = ' '.join(names).lower()
    if '26n' in joined:
        UMI_SOURCE, UMI1_LEN, UMI_MAX_ED = '32N', 26, 2
    elif any(t in joined for t in ('r1long-6n', '-6n', '_6n', '9n', 'random', 'hexamer')):
        # EXACT match for the 6-mer. Measured on a synthetic 24-plex (planted
        # truth): ed<=1 recovered 88% of planted molecules on the deep mouse lanes and
        # 97% on human ones, because a 6-mer has 18 neighbours in a 4,096-tag space and
        # same-position groups at rRNA-scale hotspots reach 500+ molecules; exact match
        # recovered 97% and 100%. The 2% of reads carrying a UMI sequencing error are the
        # price, and it is smaller than the collision cost.
        # 6N + HH = 8 nt is the UMI definition used for the benchmark. The read name still
        # carries the whole degenerate stretch (8 / 11 / 14 / 7 nt by code set); dedup uses its first 8 nt.
        # Measured on the 24-plex: -0.85% molecules vs the full stretch (0% on the 8-nt code set). The
        # full-stretch definition stays selectable for comparisons: BOBSEQ_UMI_SOURCE=6N+spacer.
        UMI_SOURCE, UMI1_LEN, UMI_MAX_ED = '6N+HH', 0, 0
        if os.environ.get('BOBSEQ_UMI_SOURCE') == '6N+spacer':
            UMI_SOURCE = '6N+spacer'
    else:
        raise SystemExit(f'umi_dedup: FATAL cannot place the UMI for rt_primers_used='
                         f'{names!r}. Known: R1_26N_polydT (32N), R1long-6N / 9N (TSO 6N). '
                         f'Refusing to guess: a wrong UMI length corrupts every duplicate '
                         f'number silently.')
    return {'umi_source': UMI_SOURCE, 'umi1_len': UMI1_LEN, 'umi_max_ed': UMI_MAX_ED}
_UMI_SPACE = {'6N': 4 ** 6, '6N+HH': 4 ** 6 * 9, '6N+spacer': 4 ** 6 * 9, '12N': 4 ** 12, '26N': 4 ** 26, '32N': 4 ** 32}
# '6N+spacer' space is the SMALLEST set's (6N + HH); the +11 and +14 sets have 8x and 64x more.
_SATURATION_WARN_GROUP = 100        # group size at which 6N loss passes ~1%

# Substitutions tolerated when deciding two UMIs are the same molecule. Sequencing
# error at Q40 is ~1e-4/base, so over 32 nt ~0.3% of UMIs carry >=1 error; allowing 2
# absorbs those without merging genuinely different UMIs (two random 32-mers differ in
# ~24 positions, so an accidental <=2 collision is vanishingly unlikely).
UMI_MAX_ED = 2                      # overridden per chemistry by configure_umi()


def cluster_umis_map(umis, max_ed=None):
    """{umi -> cluster id}, merging UMIs within `max_ed` substitutions.

    Split out from _cluster_umis so the SAME clustering can be used as a FILTER
    (keep one read per cluster), not merely counted. _cluster_umis is now a thin
    wrapper over this, so the reported unique-read count and any deduplicated
    read set are guaranteed to agree by construction rather than by coincidence.

    PIGEONHOLE BLOCKING, not all-vs-all. Split each UMI into max_ed+1 blocks: if two
    UMIs differ in <= max_ed positions then at least one block must match EXACTLY, so
    only UMIs sharing a block are ever compared. All-vs-all on 251k unique UMIs would be
    ~6e10 comparisons; blocking on ~10 nt keys spreads them over ~1e6 buckets, so the
    buckets are nearly empty and the work is close to linear."""
    if max_ed is None:
        max_ed = UMI_MAX_ED
    uniq = sorted(set(umis))
    n = len(uniq)
    parent = list(range(n))
    if n < 2 or max_ed <= 0:
        # Each UMI is its own cluster.
        return {u: i for i, u in enumerate(uniq)}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    L = len(uniq[0])
    nb = max_ed + 1
    for bi in range(nb):
        s, e = bi * L // nb, (bi + 1) * L // nb
        buckets = {}
        for i, u in enumerate(uniq):
            buckets.setdefault(u[s:e], []).append(i)
        for idxs in buckets.values():
            if len(idxs) < 2:
                continue
            for a in range(len(idxs)):
                for b in range(a + 1, len(idxs)):
                    ia, ib = idxs[a], idxs[b]
                    ra, rb = find(ia), find(ib)
                    if ra == rb:
                        continue
                    ua, ub = uniq[ia], uniq[ib]
                    d = 0
                    for x, y in zip(ua, ub):
                        if x != y:
                            d += 1
                            if d > max_ed:
                                break
                    if d <= max_ed:
                        parent[ra] = rb
        del buckets
    return {u: find(i) for i, u in enumerate(uniq)}


def _cluster_umis(umis, max_ed=None):
    """Number of distinct UMI clusters (the unique-read metric). Thin wrapper over
    cluster_umis_map so the count can never drift from the dedup filter's grouping."""
    return len(set(cluster_umis_map(umis, max_ed).values()))


def pick_umi(mq):
    """UMI string for a QNAME-parsed record, per UMI_SOURCE. '' if unavailable."""
    u2 = mq.get('umi2') or ''          # 6 nt, R2 (bobcode-side)
    u1 = mq.get('umi1') or ''          # 26 nt, R1
    if UMI_SOURCE == '6N':
        return u2 if len(u2) == 6 else ''
    if UMI_SOURCE == '6N+spacer':
        # prep wrote the whole stretch up to the G-run (8 / 11 / 14 / 7 nt by code set); the
        # length is the code's, so anything shorter than the 6N core is a failed call
        return u2 if len(u2) >= 6 else ''
    if UMI_SOURCE == '6N+HH':
        # the first 8 nt of that stretch (6N + HH; 7 for the 11-mer BOB25D): the benchmark's
        # UMI definition; prep still writes the full stretch into the QNAME
        return u2[:8] if len(u2) >= 6 else ''
    if UMI_SOURCE == '12N':
        return (u2 + u1[:6]) if (len(u2) == 6 and len(u1) >= 6) else ''
    if UMI_SOURCE == '26N':
        return u1 if len(u1) >= 26 else ''
    if UMI_SOURCE == '32N':
        return (u1 + u2) if (len(u1) >= 26 and len(u2) == 6) else ''
    return u2

# --- optical / ExAmp duplicate separation -----------------------------------
# Not all duplicates are PCR. On a PATTERNED flowcell (NovaSeq) two mechanisms put
# the SAME molecule into two clusters that are physical NEIGHBOURS:
#   * optical  - one cluster miscalled as two during imaging
#   * ExAmp    - exclusion-amplification seeds an adjacent well ("pad hopping")
# Both are sequencing artifacts, NOT library complexity, so they should not be
# reported as PCR duplication. They are separable because they are spatially
# clustered: same tile, small pixel distance. 2500 px is the Picard NovaSeq
# convention (OPTICAL_DUPLICATE_PIXEL_DISTANCE; 100 for older unpatterned flowcells).
OPTICAL_PIXEL_DIST = 2500


def _flowcell_xy(qname):
    """(tile, x, y) from an Illumina read id INSTR:RUN:FC:LANE:TILE:X:Y (prep may
    append _umi1_umi2_... which is stripped first). None if not parseable."""
    q = qname.split()[0]
    parts = q.rsplit('_', 5)                   # prep appends exactly 5 fields; the id itself may contain '_'
    base = parts[0] if len(parts) == 6 else q
    p = base.split(':')
    if len(p) < 7:
        return None
    try:
        return (p[4], int(p[5]), int(p[6]))
    except ValueError:
        return None


def _n_optical(coords):
    """Reads in a duplicate group that are a spatial neighbour of another member
    (same tile, <=OPTICAL_PIXEL_DIST px) — i.e. optical/ExAmp rather than PCR."""
    pts = [c for c in coords if c]
    if len(pts) < 2:
        return 0
    hit = set()
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            if pts[i][0] != pts[j][0]:
                continue
            if (pts[i][1] - pts[j][1]) ** 2 + (pts[i][2] - pts[j][2]) ** 2 <= OPTICAL_PIXEL_DIST ** 2:
                hit.add(i); hit.add(j)
    return max(0, len(hit) - 1) if hit else 0


def _refend(p, c):
    return p + sum(int(n) for n, o in _CIG.findall(c) if o in 'MDN=X') - 1

def _ed_le1(a, b):
    if a == b:
        return True
    la, lb = len(a), len(b)
    if abs(la - lb) > 1:
        return False
    if la == lb:
        return sum(x != y for x, y in zip(a, b)) <= 1
    if la > lb:
        a, b, la, lb = b, a, lb, la
    i = j = d = 0
    while i < la and j < lb:
        if a[i] == b[j]:
            i += 1; j += 1
        else:
            d += 1; j += 1
            if d > 1:
                return False
    return True

def _n_clusters(umis):
    """Distinct UMIs. EXACT matching, O(n).

    Merging UMIs within edit distance 1 (to absorb ONT basecalling error) would be an
    all-vs-all scan, i.e. O(n^2) in the number of distinct UMIs per gene. Not done here:
    Illumina UMI bases are Q40 (~1e-4 per base, ~0.1% chance of any error in 12 nt), so
    error-merging buys almost nothing and costs quadratic time at depth. It also
    actively over-merges: at 12 nt, ed<=1 has 36 neighbours per UMI, so distinct
    molecules get collapsed."""
    return len(set(umis))

def _mrna_gene(ch, sx, e, genes, gst, exons, est, ridna, rdst):
    """(species, gene_name) of the max-overlap protein-coding gene for an exonic,
    non-rRNA, non-MT read, else None."""
    if ch.endswith('_MT'):
        return None
    v = ridna.get(ch)
    if v:
        i = bisect.bisect_right(rdst[ch], e)
        for j in range(max(0, i - 1), -1, -1):
            a, b = v[j]
            if b < sx:
                break
            if a < e and b > sx:
                return None                              # rRNA
    gs = st._overlap(genes, gst, ch, sx, e, collect=True)
    if not gs or not st._overlap(exons, est, ch, sx, e):
        return None
    best = max(gs, key=lambda g: min(g[1], e) - max(g[0], sx))
    if best[2] != 'protein_coding' or not best[3]:
        return None
    return ('h' if ch.startswith('HUMAN_') else 'm', best[3])

def _n_mol(coords, W):
    """Distinct molecules among reads already sharing (gene, UMI[, strand]).

    W == 0 (ILLUMINA): coordinates are exact, so "same molecule" is just "same
    (start, end)" and the answer is len(set(coords)) — O(n). The union-find below is
    O(n^2) (every read compared to every other) which is invisible at subsample size
    but explodes with depth: a group of 70,000 reads would need ~2.5e9 comparisons.
    W > 0 (ONT): reads of one molecule truncate at variable points, so the +/-W
    union-find is still required to absorb that jitter."""
    n = len(coords)
    if n == 1:
        return 1
    if not W:
        return len(set(coords))
    par = list(range(n))
    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for i in range(n):
        for j in range(i + 1, n):
            if abs(coords[i][0] - coords[j][0]) <= W and abs(coords[i][1] - coords[j][1]) <= W:
                par[find(i)] = find(j)
    return len({find(i) for i in range(n)})


def estimate_library_size(n_reads, n_unique):
    """Estimated total distinct molecules in the LIBRARY, from the reads sequenced and
    the unique molecules seen (the Lander-Waterman / Picard EstimateLibraryComplexity
    model). Solves  U = C * (1 - exp(-N/C))  for C by bisection.

    Intuition: if you sequence N reads and only see U distinct molecules, the amount of
    re-sampling tells you how big the pool must be. No duplicates -> pool is effectively
    unbounded (returns None); heavy duplication -> pool is small.

    CAVEAT: the model assumes every molecule is equally likely to be sampled. RNA-seq is
    the opposite — expression spans orders of magnitude, so highly expressed transcripts
    duplicate long before rare ones do. That inflates the apparent duplicate rate and
    makes this a LOWER BOUND on true complexity, not a point estimate."""
    N, U = float(n_reads), float(n_unique)
    if N <= 0 or U <= 0 or U >= N:
        return None                      # no duplication -> cannot bound the pool
    def f(C):
        return C * (1.0 - math.exp(-N / C)) - U
    lo, hi = U, U * 2.0
    for _ in range(200):                 # expand until the root is bracketed
        if f(hi) > 0:
            break
        lo, hi = hi, hi * 2.0
    for _ in range(200):                 # bisect
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def saturation_curve(dup_hist, n_points=24):
    """Molecules recovered against reads sampled, computed EXACTLY from the duplicity
    histogram. Returns [(reads_sampled, expected_distinct_molecules), ...].

    WHY THIS AND NOT THE COMPLEXITY NUMBER. `estimate_library_size` assumes every
    molecule is equally likely to be sampled; RNA-seq expression spans orders of
    magnitude, so that estimate is a LOWER BOUND whose bias depends on how skewed the
    expression is. Comparing two such bounds ACROSS assays therefore compares their
    skew as much as their libraries. This curve assumes nothing about abundance: it is
    the exact expectation of how many distinct molecules d sampled reads would recover,
    given the observed distribution, so two libraries can be read off at a common depth.

    NO SUBSAMPLING, NO RNG. For a molecule seen k times out of N reads, the chance d
    reads all miss it is C(N-k, d) / C(N, d), evaluated in log-gamma space. Molecules
    with the same k share a term, so the whole curve costs O(len(hist)) per point and
    is bit-for-bit reproducible.

    dup_hist: {reads_per_molecule: how_many_molecules} (ints, or the string keys that
    survive a JSON round-trip).
    """
    hist = {}
    for k, m in (dup_hist or {}).items():
        try:
            k = int(k)
        except (TypeError, ValueError):
            continue
        if k > 0 and m:
            hist[k] = hist.get(k, 0) + int(m)
    if not hist:
        return []
    N = sum(k * m for k, m in hist.items())          # total reads
    total_mol = sum(hist.values())                   # molecules actually observed
    if N <= 0 or total_mol <= 0:
        return []

    lg = math.lgamma
    def expected(d):
        if d >= N:
            return float(total_mol)                  # sampling everything finds everything
        base = lg(N - d + 1) - lg(N + 1)
        out = 0.0
        for k, m in hist.items():
            rem = N - k                              # reads belonging to other molecules
            if d > rem:
                out += m                             # cannot avoid this molecule
            else:
                out += m * (1.0 - math.exp(lg(rem + 1) - lg(rem - d + 1) + base))
        return out

    # log-spaced depths so the early, informative part of the curve is not crushed
    lo = max(1, N // 10000)          # not a hard 1000: small inputs must still curve
    if lo >= N:
        return [(N, float(total_mol))]
    pts, seen = [], set()
    for i in range(n_points):
        d = int(round(lo * (N / lo) ** (i / (n_points - 1))))
        d = min(max(d, 1), N)
        if d not in seen:
            seen.add(d)
            pts.append((d, expected(d)))
    if pts[-1][0] != N:
        pts.append((N, float(total_mol)))
    return pts


def saturation_summary(dup_hist):
    """Coverage and marginal yield from the duplicity histogram. Returns {} if empty.

    These are GOOD-TURING quantities, and they are exact rather than modelled. Of the
    next read sequenced, the chance it lands on a molecule never seen before is f1/N,
    where f1 is the number of molecules observed exactly ONCE. So:

        coverage        = 1 - f1/N     fraction of the molecule POOL already sampled
        marginal yield  = f1/N         new molecules per additional read

    WHY THIS MATTERS MORE THAN THE COMPLEXITY NUMBER. `estimate_library_size`
    (Lander-Waterman) assumes uniform sampling. On real RNA-seq it can report a library
    as exhausted while a third of the molecules seen have been seen only once, which is
    the signature of a pool still far from sampled out. Measured on one library:
    LW estimated 1.17M against 1.17M observed (i.e. "done"), while f1 = 378,751 and
    2.29% of reads were still hitting new molecules. Good-Turing needs no abundance
    model, so it does not have that failure mode and is the number to compare ACROSS
    assays. LW is reported alongside for comparability with other tools.
    """
    hist = {}
    for k, m in (dup_hist or {}).items():
        try:
            k = int(k)
        except (TypeError, ValueError):
            continue
        if k > 0 and m:
            hist[k] = hist.get(k, 0) + int(m)
    if not hist:
        return {}
    N = sum(k * m for k, m in hist.items())
    U = sum(hist.values())
    if N <= 0 or U <= 0:
        return {}
    f1 = hist.get(1, 0)
    pts = saturation_curve(hist)
    # Reads to reach half the molecules: the EXACT crossing of the expectation curve,
    # found by bisection on d. Reading it off the 24 log-spaced grid depths (first depth at
    # or above half) would overstate the answer by up to one grid step (4% on a 100k
    # fixture, 12-23% on real libraries).
    def _expected(d):
        lg = math.lgamma; tot = 0.0
        for k, m in hist.items():
            p_miss = math.exp(lg(N - k + 1) - lg(N - k - d + 1) - lg(N + 1) + lg(N - d + 1)) if d <= N - k else 0.0
            tot += m * (1.0 - p_miss)
        return tot
    half = None
    if pts and _expected(N) >= 0.5 * U:
        lo, hi = 0, N
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if _expected(mid) < 0.5 * U: lo = mid
            else: hi = mid
        half = hi
    return {
        'n_reads': N,
        'molecules_observed': U,
        'singletons': f1,
        'singleton_pct': round(100.0 * f1 / U, 2),
        'coverage_pct': round(100.0 * (1.0 - f1 / N), 3),
        'new_molecules_per_1k_reads': round(1000.0 * f1 / N, 2),
        'reads_to_half_of_observed': half,
    }


def compute_dup_funnel(bam, gtf, rdna, W=10):
    """Duplicate rate, computed directly: of the aligned reads carrying a callable UMI,
    what fraction share their UMI with another read.

    Deliberately simple. A 4-step funnel (non-ribosomal mRNA -> same gene -> same
    start/end -> same UMI) is not used: measured on a 1M-read subset those steps all
    land within a few points of each other (68.3 / 70.1 / 72.7%), so the gene and
    position layers add complexity without changing the answer — and the gene lookup
    is the slowest part of the whole module.

    gtf/rdna are accepted but unused (kept so callers do not change)."""
    umis = collections.Counter()   # reads per UMI: gives the aligned-scope histogram
    fam = collections.defaultdict(list)
    n = 0
    p = subprocess.Popen(['samtools', 'view', '-F', '0x904', bam],
                         stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        f = line.split('\t')
        ch = f[2]
        if not (ch.startswith('HUMAN_') or ch.startswith('MOUSE_')):
            continue
        if _positional():
            mq = ox.measure_from_qname(f[0])
            if mq is None:
                continue
            u = pick_umi(mq)
        else:
            m = ox.measure_struct_24plex(f[9])
            if m is None or m.get('code') is None:
                continue
            sread = f[9] if m['orient'] == 'F' else ox.rc(f[9])
            u = sread[m['bc_pos'] + 7: m['bc_pos'] + 13]
        if not u:
            continue
        n += 1
        umis[u] += 1
        # (position, UMI) identifies a duplicate family; flowcell coords let us ask
        # whether the copies are physical NEIGHBOURS (optical/ExAmp) or unrelated (PCR).
        _xy = _flowcell_xy(f[0])
        if _xy:
            fam[(ch, int(f[3]), u)].append(_xy)
    p.wait()
    if not n:
        return None
    # Cluster once and reuse the mapping, so the unique count and the duplicity
    # histogram are the same clustering by construction rather than by coincidence.
    _cmap = cluster_umis_map(umis.keys())
    _per_cluster = collections.Counter()
    for _u, _c in umis.items():
        _per_cluster[_cmap[_u]] += _c
    uniq = len(_per_cluster)
    _hist_aligned = collections.Counter(_per_cluster.values())
    dup = n - uniq
    # FLOWCELL duplicates: copies sitting within OPTICAL_PIXEL_DIST on the SAME tile.
    # Optical (one cluster imaged as two) and NovaSeq ExAmp pad-hopping both put the same
    # molecule in a NEIGHBOURING well, so they are spatially clustered; genuine PCR copies
    # land at unrelated flowcell positions. NOTE this can only be seen in CONTIGUOUS read
    # data — a strided/random subsample almost never contains both members of a neighbour
    # pair and will report ~0 here regardless of the true rate.
    opt = sum(_n_optical(v) for v in fam.values() if len(v) > 1)
    lib = estimate_library_size(n, uniq)
    return {'n_umi': n, 'unique': uniq, 'dup_reads': dup,
            'dup_pct': round(100 * dup / n, 2),
            'lib_size': lib, 'lib_size_M': (round(lib / 1e6, 3) if lib else None),
            'dup_flowcell': opt, 'dup_flowcell_pct': round(100 * opt / n, 2),
            'optical_px': OPTICAL_PIXEL_DIST,
            'umi_source': UMI_SOURCE, 'umi_max_ed': UMI_MAX_ED,
            # ALIGNED-scope POSITIONAL histogram: clusters keyed on (contig, POS, UMI
            # family), the key this function's dup_pct / flowcell rows use. It is NOT the
            # histogram the report's complexity rows and saturation plot read: those prefer
            # 'duplicity_hist', which dedup_reads_illumina records over the same aligned
            # reads under definition D (bobcode + contig + strand + 5' anchor + UMI ed), so
            # that the duplicate rate, the complexity rows and the curve share one
            # definition. This one is kept as a diagnostic fallback.
            'duplicity_hist_aligned': {str(k): v for k, v in sorted(_hist_aligned.items())}}


def compute_gene_umi_dup(bam, gtf, rdna):
    genes, gst, exons, est = gtf
    ridna, rdst = rdna
    per_gene = collections.defaultdict(list)
    n = 0
    p = subprocess.Popen(['samtools', 'view', '-F', '0x904', bam],
                         stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        f = line.split('\t'); ch = f[2]
        if not (ch.startswith('HUMAN_') or ch.startswith('MOUSE_')):
            continue
        pos = int(f[3]); e = _refend(pos, f[5])
        g = _mrna_gene(ch, pos - 1, e, genes, gst, exons, est, ridna, rdst)
        if g is None:
            continue
        # Same QNAME-vs-SEQ split as compute_dup_funnel: on Illumina the UMI has been
        # trimmed out of the aligned SEQ, so it must come from the read name, and the
        # 26N R1 UMI must lead because the 6 nt UMI saturates (4^6) at this depth.
        if _positional():
            mq = ox.measure_from_qname(f[0])
            if mq is None:
                continue
            umi = pick_umi(mq)                   # composition set by UMI_SOURCE
            if umi:
                per_gene[g].append(umi)
                n += 1
            continue
        m = ox.measure_struct_24plex(f[9])
        if m is None or m.get('code') is None:
            continue
        strand = f[9] if m['orient'] == 'F' else ox.rc(f[9])
        umi = strand[m['bc_pos'] + 7: m['bc_pos'] + 13]     # 6 nt UMI after the 7-mer
        if len(umi) == 6:
            per_gene[g].append(umi)
            n += 1
    p.wait()
    if n == 0:
        return None
    def _summ(lists):
        nn = sum(len(us) for us in lists)
        uu = sum(_n_clusters(us) for us in lists)
        return {'n_mrna': nn, 'unique': uu, 'dup_reads': nn - uu,
                'dup_pct': round(100 * (nn - uu) / nn, 1) if nn else None}
    res = _summ(list(per_gene.values()))                          # combined (human + mouse)
    res['human'] = _summ([us for k, us in per_gene.items() if k[0] == 'h'])
    res['mouse'] = _summ([us for k, us in per_gene.items() if k[0] == 'm'])
    return res

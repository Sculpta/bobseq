#!/usr/bin/env python3
"""dedup_reads_illumina.py — deduplicate a barcode's reads BEFORE any analysis, on definition D.

The key is (bobcode, contig, strand, 5' position, UMI ed<=N); see build_keepset. With
paired-end alignment (random-priming chemistry) the key also carries the fragment length
|TLEN|, which with the R2 anchor fixes the priming site, and the UMI is the whole degenerate
stretch. The analysis BAM holds mate-1 (R2) records only; the full pair BAM
(<bc>_star_pairs.bam), when present, is filtered by the same read names and measured for
insert length into the sidecar (block `insert_size`).


WHY THIS EXISTS
    Without this step duplication would be a REPORTED METRIC only: every stage
    (composition, barcode accuracy, transcript diversity, top genes) would count PCR copies
    as independent observations. At typical depths that is a large effect — measured on a
    24-plex library, 95.7% of aligned reads at 10M are duplicates, so an abundant molecule
    that amplified well is counted thousands of times while a rare one is counted once.
    Deduplicating first makes every downstream number ONE VOTE PER MOLECULE.

WHAT IT DOES
    1. Computes the duplicate metrics on the FULL data and writes them to a sidecar,
       BEFORE filtering. This is not optional bookkeeping: run the existing metric on an
       already-deduplicated BAM and it reports ~0% duplicates by construction, which
       would look like a pristine library rather than a filtered one.
    2. Groups reads by molecule (definition D) and keeps ONE read per group.
    3. Rewrites the BAM and both FASTQs in place, preserving the originals as `.all.`.

TWO DECISIONS, both deliberate:

  * COLLAPSE IS BY POSITION + UMI (definition D), not globally by UMI. A global,
    position-independent UMI cluster (the grouping umi_dedup.compute_dup_funnel also
    reports) is only safe for a long UMI: with a 32 nt UMI and ~0.4M molecules chance
    collisions are negligible, but with a 6 nt UMI it would collapse most of a library
    into 4,096 tags. Position is what tells two molecules with the same short tag apart.

  * THE SURVIVING READ IS THE LEXICOGRAPHICALLY SMALLEST READ ID in its cluster, NOT the
    first encountered. "First" would depend on BAM record order, which depends on STAR's
    --runThreadN. A min-id rule is reproducible across thread counts and reruns.

    NB: the representative is chosen WITHOUT regard to alignment quality. Within a UMI
    cluster the copies are the same molecule, so this is a coin flip between near-identical
    records; it is not a quality-based pick and should not be described as one.

THE KEEP-SET IS BUILT FROM THE FASTQ, NOT THE BAM. compute_dup_funnel only sees reads
aligned to HUMAN_/MOUSE_ contigs; keying the filter off that would silently delete every
unmapped read, changing the mapping rate into an artifact of deduplication. Every prepped
read carries its UMI in the QNAME, so the FASTQ is the honest universe.

USAGE:
  python3 dedup_reads_illumina.py <run_dir> [barcode01 ...] [--force]
"""
import json
import re
import os
import subprocess
import sys

SELF = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SELF)


def _qname_umi(rid, ox, ud):
    mq = ox.measure_from_qname(rid)
    if mq is None:
        return None
    return ud.pick_umi(mq) or None


def _ref_span(cigar):
    """Reference bases consumed by a CIGAR string (M/D/N/=/X)."""
    n = 0
    for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        if op in 'MDN=X':
            n += int(length)
    return n


def build_keepset(bam, prepped_fq, ox, ud):
    """-> (keepset, n_reads, n_clusters, duplicity_hist, meta). One read per MOLECULE,
    where a molecule is definition D of the benchmark:

        (bobcode, contig, strand, 5'-end position, [|TLEN|,] UMI within UMI_MAX_ED substitutions)

    |TLEN| enters the key on paired alignments (flag 0x1): together with the R2 anchor it
    fixes the mate's 5' end, the random-priming site, so one molecule = one template-switch
    point + one priming site + one UMI. It is 0 when the mate did not align, which reduces
    to the single-end key for that fragment (a copy whose R1 failed to align is therefore
    not merged with its mates-mapped siblings; the same convention as Picard's).

    A global position-independent UMI cluster is only defensible for a 32 nt UMI; with the
    TSO 6N as the only UMI (R1long-6N) it would collapse most of a library into 4,096 tags.
    Position is what tells two molecules with the same short tag apart, and these are 3'
    end methods, so same position = one capture event.

    Mapped reads: grouped on the BAM's primary alignment. The anchor is the read's 5' end
    (POS for forward reads, alignment end for reverse reads), so a clipped 3' end cannot
    move it. Unmapped reads have no position and fall back to (bobcode, UMI) clustering,
    counted separately in meta. Two passes over the BAM so only integers are held per read.
    """
    max_ed = ud.UMI_MAX_ED
    groups = {}          # (bobcode, contig, strand, anchor, |tlen|) -> [(umi, recno)]
    unmapped = {}        # (bobcode, umi) -> min recno   (position unknown)
    n = n_mapped = n_unmapped = n_nocall = 0
    paired = None
    both = r2_only = r1_only = pair_unmapped = 0
    rd = subprocess.Popen(['samtools', 'view', '-F', '0x900', bam],
                          stdout=subprocess.PIPE, text=True)
    for recno, line in enumerate(rd.stdout):
        n += 1
        f = line.split('\t', 9)
        mq = ox.measure_from_qname(f[0])
        u = ud.pick_umi(mq) if mq else None
        if not u:
            n_nocall += 1
            continue
        bc = mq['seven']
        flag = int(f[1])
        if paired is None:
            paired = bool(flag & 1)
        if paired:
            # the analysis BAM holds mate-1 (R2) records only; the mate's fate is in the flags
            if flag & 4:
                if flag & 8: pair_unmapped += 1
                else:        r1_only += 1
            elif flag & 8:   r2_only += 1
            else:            both += 1
        if flag & 4:
            n_unmapped += 1
            key = (bc, u)
            if key not in unmapped:
                unmapped[key] = recno
            continue
        n_mapped += 1
        pos = int(f[3])
        anchor = pos + _ref_span(f[5]) - 1 if flag & 16 else pos
        tlen = abs(int(f[8])) if paired else 0
        groups.setdefault((bc, f[2], flag & 16, anchor, tlen), []).append((u, recno))
    rd.wait()

    keep_rec = set()
    hist = {}
    for members in groups.values():
        if len(members) == 1:
            keep_rec.add(members[0][1]); hist[1] = hist.get(1, 0) + 1
            continue
        cmap = ud.cluster_umis_map([u for u, _ in members], max_ed)
        best, size = {}, {}
        for u, r in members:
            c = cmap[u]
            size[c] = size.get(c, 0) + 1
            if c not in best or r < best[c]:
                best[c] = r
        keep_rec.update(best.values())
        for k in size.values():
            hist[k] = hist.get(k, 0) + 1
    n_clusters = sum(hist.values())
    # Definition-D molecules among ALIGNED reads = the groups above (every group is a
    # mapped anchor). This is the one duplicate rate the report shows.
    kept_aligned = len(keep_rec)
    # unmapped: one per (bobcode, UMI), exact -- there is no position to do better with
    keep_rec.update(unmapped.values())

    keep = set()
    rd = subprocess.Popen(['samtools', 'view', '-F', '0x900', bam],
                          stdout=subprocess.PIPE, text=True)
    for recno, line in enumerate(rd.stdout):
        if recno in keep_rec:
            keep.add(line[:line.index('\t')])
    rd.wait()
    meta = {'dedup_key': ('bobcode+contig+strand+5prime_position+fragment_length+UMI' if paired
                          else 'bobcode+contig+strand+5prime_position+UMI'),
            'paired': bool(paired),
            'umi_source': ud.UMI_SOURCE, 'umi_max_ed': max_ed,
            'bam_primary': n, 'mapped_with_umi': n_mapped, 'unmapped_with_umi': n_unmapped,
            'no_umi_call': n_nocall, 'unmapped_kept': len(unmapped),
            # definition D over aligned reads: molecules / reads -> the reported duplicate rate
            'dedup_kept_aligned': kept_aligned,
            'dup_pct_defD_aligned': round(100.0 * (1 - kept_aligned / n_mapped), 3) if n_mapped else None}
    if paired:
        # fragments by mate fate (records with a callable UMI): "mapped" everywhere in the
        # report means R2 aligned, i.e. both + r2_only
        meta.update({'mate_pairs_both_mapped': both, 'r2_only_mapped': r2_only,
                     'r1_only_mapped': r1_only, 'pairs_unmapped': pair_unmapped})
    if not keep:
        return None, n, 0, {}, meta
    return keep, n, n_clusters, hist, meta


_CIG_ALL = re.compile(r'(\d+)([MIDNSHP=X])')


def _cigar_walk(pos, cigar):
    """-> (ref_end, leading_S, trailing_S, query_aligned, [(intron_start, intron_end)]) with
    1-based inclusive reference coordinates; introns are the N operations."""
    ops = _CIG_ALL.findall(cigar)
    lead = int(ops[0][0]) if ops and ops[0][1] == 'S' else 0
    trail = int(ops[-1][0]) if len(ops) > 1 and ops[-1][1] == 'S' else 0
    ref = pos; qa = 0; introns = []
    for L, op in ops:
        L = int(L)
        if op in 'M=X':
            ref += L; qa += L
        elif op == 'I':
            qa += L
        elif op == 'D':
            ref += L
        elif op == 'N':
            introns.append((ref, ref + L - 1)); ref += L
    return ref - 1, lead, trail, qa, introns


INSERT_CAP = 1500     # histogram bin ceiling; longer inserts are counted in the top bin


def insert_summary(block):
    """Median / quartiles / mean and the shares from the histogram, so a merged sidecar and a
    single one are summarised by the same code. Mutates and returns `block`."""
    hist = {int(k): int(v) for k, v in (block.get('hist') or {}).items()}
    n = sum(hist.values())
    block['n'] = n
    if not n:
        block.update({'median': None, 'q1': None, 'q3': None, 'mean': None})
        return block
    acc = 0; qs = {}
    for k in sorted(hist):
        acc += hist[k]
        for name, q in (('q1', 0.25), ('median', 0.5), ('q3', 0.75)):
            if name not in qs and acc >= q * n:
                qs[name] = k
    block.update(qs)
    block['mean'] = round(sum(k * v for k, v in hist.items()) / n, 1)
    for k in ('n_mates_overlap', 'n_shorter_than_r1', 'n_shorter_than_r2'):
        block[f'pct_{k[2:]}'] = round(100.0 * block.get(k, 0) / n, 2)
    return block


def insert_size_stats(pairs_bam):
    """Insert length per fragment from the FULL pair BAM, pre-dedup (one value per sequenced
    fragment, as asked: 'insert length per read'). Both mates primary and mapped to the same
    contig with TLEN != 0. In cDNA space:

        insert = |TLEN| - intron bases spanned by either mate (union of N operations)
                 + R2 5' soft clip beyond the G-run + R1 5' soft clip

    The outer clips are insert (untemplated G-run excluded via the QNAME's G-run length; the
    priming hexamer's mismatches are clipped by STAR but are cDNA); inner clips are adapter
    read-through, not insert. Mates that do not overlap can hide an intron between them, so
    their share is reported (n_no_overlap) next to the median rather than corrected."""
    hist = {}
    n_overlap = n_short_r1 = n_r2_short = 0
    r1_lens = {}
    rd = subprocess.Popen(['samtools', 'view', '-F', '0x900', pairs_bam],
                          stdout=subprocess.PIPE, text=True)
    pending = None
    for line in rd.stdout:
        f = line.split('\t', 10)
        if pending is None or pending[0] != f[0]:
            pending = f
            continue
        m1, m2 = (pending, f) if int(pending[1]) & 0x40 else (f, pending)
        pending = None
        fl1, fl2 = int(m1[1]), int(m2[1])
        if (fl1 & 4) or (fl2 & 4) or m1[2] != m2[2]:
            continue
        tlen = abs(int(m1[8]))
        if tlen == 0:
            continue
        p1, p2 = int(m1[3]), int(m2[3])
        e1, lead1, trail1, qa1, in1 = _cigar_walk(p1, m1[5])
        e2, lead2, trail2, qa2, in2 = _cigar_walk(p2, m2[5])
        # union of intron intervals within the fragment span
        nlen = 0; cur = None
        for a, b in sorted(in1 + in2):
            if cur is None or a > cur[1] + 1:
                if cur: nlen += cur[1] - cur[0] + 1
                cur = [a, b]
            else:
                cur[1] = max(cur[1], b)
        if cur: nlen += cur[1] - cur[0] + 1
        clip1 = trail1 if fl1 & 16 else lead1          # R2 5' end
        clip2 = trail2 if fl2 & 16 else lead2          # R1 5' end
        try:
            glen = int(m1[0].rsplit('_', 5)[4])
        except (IndexError, ValueError):
            glen = 0
        ins = tlen - nlen + max(0, clip1 - glen) + clip2
        ins = max(1, min(ins, INSERT_CAP))
        hist[ins] = hist.get(ins, 0) + 1
        if tlen - nlen < qa1 + qa2:
            n_overlap += 1
        r1len = len(m2[9]); r1_lens[r1len] = r1_lens.get(r1len, 0) + 1
        if ins < r1len:
            n_short_r1 += 1
        if ins < len(m1[9]):
            n_r2_short += 1
    rd.wait()
    r1_mode = max(r1_lens, key=r1_lens.get) if r1_lens else None
    return insert_summary({'scope': 'pre-deduplication, both mates aligned to one contig, per fragment',
                           'definition': '|TLEN| - introns spanned by either mate + outer soft clips (R2 beyond the G-run, R1)',
                           'cap': INSERT_CAP, 'r1_read_len_mode': r1_mode,
                           'hist': {str(k): v for k, v in sorted(hist.items())},
                           'n_mates_overlap': n_overlap, 'n_shorter_than_r1': n_short_r1,
                           'n_shorter_than_r2': n_r2_short})


def filter_fastq(src, dst, keep, key=lambda r: r):
    kept = 0
    with open(src) as fi, open(dst, 'w') as fo:
        while True:
            h = fi.readline()
            if not h:
                break
            s = fi.readline(); p = fi.readline(); q = fi.readline()
            if key(h[1:].split()[0]) in keep:
                fo.write(h); fo.write(s); fo.write(p); fo.write(q)
                kept += 1
    return kept


def filter_bam(src, dst, keep):
    """samtools 1.20 here has no -N/--qname-file, so filter the stream ourselves."""
    rd = subprocess.Popen(['samtools', 'view', '-h', src], stdout=subprocess.PIPE, text=True)
    wr = subprocess.Popen(['samtools', 'view', '-b', '-o', dst, '-'],
                          stdin=subprocess.PIPE, text=True)
    kept = 0
    for line in rd.stdout:
        if line.startswith('@'):
            wr.stdin.write(line)
            continue
        if line[:line.index('\t')] in keep:
            wr.stdin.write(line)
            kept += 1
    wr.stdin.close(); rd.wait(); wr.wait()
    return kept


def dedup_barcode(run_dir, bc, force=False):
    import star_taxonomy_illumina as st
    import umi_dedup_illumina as ud
    import override_24plex_illumina as ox
    st._get_adapter()          # MUST run first: configures ox, else measure_from_qname
    _guides = [os.path.join(run_dir, 'machine_readable_guide', g)
               for g in os.listdir(os.path.join(run_dir, 'machine_readable_guide'))
               if g.endswith('_analysisguide.txt')]
    _rt = ''
    for _g in _guides:
        infl = False
        for ln in open(_g):
            t = ln.strip()
            if t == 'BEGIN_CONFIG': infl = True; continue
            if t == 'END_CONFIG': break
            if infl and ln.startswith('rt_primers_used\t'):
                _rt = ln.rstrip('\n').split('\t', 1)[1]
    if not _rt:
        sys.exit('dedup_reads: FATAL no rt_primers_used in the guide; UMI length undefined')
    _umi = ud.configure_umi(_rt)
    print(f"  {bc}: UMI from chemistry {_rt} -> {_umi}")
                               # returns None for EVERY read and the keep-set is empty.

    fqd = os.path.join(run_dir, 'fastq', bc)
    prepped = os.path.join(fqd, f'prepped_{bc}.fastq')
    noad = os.path.join(fqd, f'noadapter_combined_{bc}.fastq')
    bam = os.path.join(run_dir, 'individual-analyses', bc, f'{bc}_star_combined.bam')
    side = os.path.join(run_dir, 'individual-analyses', bc, f'{bc}_dup_prefilter.json')
    # paired-end (random priming): the full pair BAM and the R1 mate FASTQ, filtered by the
    # same read names; both optional so the poly-dT path and old run dirs are untouched
    pairs = os.path.join(run_dir, 'individual-analyses', bc, f'{bc}_star_pairs.bam')
    mate = os.path.join(fqd, f'prepped_R1_{bc}.fastq')
    has_pairs = os.path.exists(pairs) and subprocess.run(
        ['samtools', 'view', '-c', pairs], capture_output=True, text=True).stdout.strip() not in ('', '0')
    has_mate = os.path.exists(mate) and os.path.getsize(mate) > 0

    if os.path.exists(prepped + '.all') and not force:
        print(f'  {bc}: already deduplicated (found {os.path.basename(prepped)}.all) — skipping')
        return None
    for p in (prepped, bam):
        if not os.path.exists(p):
            print(f'  {bc}: missing {p} — skipping'); return None

    # 1. PRE-FILTER metrics, on the untouched BAM. Must happen before any filtering.
    print(f'  {bc}: computing pre-filter duplicate metrics …')
    gtf = rdna = None
    df = ud.compute_dup_funnel(bam, gtf, rdna)
    if df:
        df['scope'] = 'pre-deduplication (aligned reads, HUMAN_/MOUSE_ only)'

    # 2. keep-set from the FASTQ (all prepped reads, mapped or not)
    print(f'  {bc}: clustering UMIs …')
    keep, n_reads, n_clusters, dup_hist, meta = build_keepset(bam, prepped, ox, ud)
    if not keep:
        print(f'  {bc}: no callable UMIs — NOT filtering'); return None
    if df is not None:
        # NB `fastq_reads` is the BAM primary-record count (== prepped reads, since STAR
        # keeps unmapped reads Within). Key name kept for readers; see 'bam_primary'.
        df['fastq_reads'] = n_reads
        df['fastq_unique'] = n_clusters
        df['dedup_kept'] = len(keep)
        # keys are strings in JSON; the reader casts back to int.
        df['duplicity_hist'] = {str(k): v for k, v in sorted(dup_hist.items())}
        df.update(meta)
        df['scope'] = 'pre-deduplication, aligned reads, definition D (same 5-prime position)'
        if has_pairs:
            print(f'  {bc}: insert length from the pair BAM …')
            df['insert_size'] = insert_size_stats(pairs)
            print(f"  {bc}: insert median {df['insert_size']['median']} nt over "
                  f"{df['insert_size']['n']:,} both-mapped fragments")
    with open(side, 'w') as f:
        json.dump(df or {}, f, indent=1)
    print(f'  {bc}: {n_reads:,} prepped reads -> {len(keep):,} unique '
          f'({100*len(keep)/n_reads:.1f}%)')

    # 3. rewrite, keeping originals
    os.rename(prepped, prepped + '.all')
    k1 = filter_fastq(prepped + '.all', prepped, keep)
    os.rename(bam, bam + '.all')
    k2 = filter_bam(bam + '.all', bam, keep)
    if has_pairs:
        os.rename(pairs, pairs + '.all')
        kp = filter_bam(pairs + '.all', pairs, keep)
        print(f'  {bc}: pair BAM {kp:,} records kept (both mates of every kept fragment)')
    if has_mate:
        os.rename(mate, mate + '.all')
        km = filter_fastq(mate + '.all', mate, keep)
        print(f'  {bc}: R1 mate FASTQ {km:,} reads kept')
        if km != k1:
            sys.exit(f'dedup_reads: FATAL {bc} kept {k1} R2 reads but {km} R1 mates; the mate '
                     f'FASTQ is out of step with the prepped reads')
    k3 = None
    if os.path.exists(noad):
        os.rename(noad, noad + '.all')
        # raw R2 names have no structure suffix; the prepped id is <origid>_<umi…>,
        # so match on the original id prefix.
        # Strip exactly the five fields prep appends (_umi1_umi2_bobcode_grun_rt), from the
        # RIGHT: instrument read ids can themselves contain '_' (some providers emit ids
        # like 'XXX-1830_XXX-2409'), so split('_')[0] would match nothing, leaving the
        # no-adapter FASTQ empty and analyze_barcode dividing by zero.
        pref = {r.rsplit('_', 5)[0] for r in keep}
        k3 = filter_fastq(noad + '.all', noad, pref)
    print(f'  {bc}: prepped {k1:,} reads · BAM {k2:,} records'
          + (f' · noadapter {k3:,} reads' if k3 is not None else ''))
    return {'n_reads': n_reads, 'kept': len(keep), 'bam_records': k2}



# ---- per-sample mode: one library sidecar from per-sample sidecars -----------------------
# Definition D keys on the bobcode, so the per-sample deduplications are exactly the pooled
# one split by sample: every count and every duplicity histogram adds, and the rates the
# report shows (definition-D duplicate rate, flowcell rate) are recomputed from the sums with
# the formulas of build_keepset / compute_dup_funnel. The two fields of the legacy UMI-only
# definition that are NOT additive (unique, dup_reads) are summed anyway and flagged in
# `merged_from`: they are logged, never shown, and a global UMI clustering across samples
# has no meaning once the UMI is the 6N.
_SIDE_COUNTS = ('n_umi', 'unique', 'dup_reads', 'dup_flowcell', 'fastq_reads', 'fastq_unique',
                'dedup_kept', 'bam_primary', 'mapped_with_umi', 'unmapped_with_umi',
                'no_umi_call', 'unmapped_kept', 'dedup_kept_aligned',
                'mate_pairs_both_mapped', 'r2_only_mapped', 'r1_only_mapped', 'pairs_unmapped')
_SIDE_HISTS = ('duplicity_hist_aligned', 'duplicity_hist')
_SIDE_CONSTANTS = ('optical_px', 'umi_source', 'umi_max_ed', 'scope', 'dedup_key', 'paired')
_INSERT_COUNTS = ('n_mates_overlap', 'n_shorter_than_r1', 'n_shorter_than_r2')


def merge_sidecars(parts):
    import umi_dedup_illumina as ud
    parts = [q for q in parts if q and not q.get('empty_sample')]
    if not parts:
        return {'n_reads': 0, 'dedup_kept': 0, 'empty_sample': True}
    out = {}
    for k in _SIDE_COUNTS:
        if any(k in q for q in parts):
            out[k] = sum(int(q.get(k, 0)) for q in parts)
    n = out.get('n_umi', 0)
    out['dup_pct'] = round(100 * out.get('dup_reads', 0) / n, 2) if n else None
    lib = ud.estimate_library_size(n, out.get('unique', 0)) if n else None
    out['lib_size'] = lib
    out['lib_size_M'] = (round(lib / 1e6, 3) if lib else None)
    out['dup_flowcell_pct'] = round(100 * out.get('dup_flowcell', 0) / n, 2) if n else None
    for k in _SIDE_CONSTANTS:
        vals = [q.get(k) for q in parts if k in q]
        if vals and any(v != vals[0] for v in vals):
            sys.exit(f'dedup merge-sidecars: samples disagree on {k}: {sorted(set(map(str, vals)))}')
        if vals:
            out[k] = vals[0]
    for k in _SIDE_HISTS:
        h = {}
        for q in parts:
            for b, v in (q.get(k) or {}).items():
                h[int(b)] = h.get(int(b), 0) + int(v)
        out[k] = {str(b): v for b, v in sorted(h.items())}
    m = out.get('mapped_with_umi', 0)
    out['dup_pct_defD_aligned'] = (round(100.0 * (1 - out.get('dedup_kept_aligned', 0) / m), 3)
                                   if m else None)
    ins = [q['insert_size'] for q in parts if q.get('insert_size')]
    if ins:
        h = {}
        for b in ins:
            for k, v in (b.get('hist') or {}).items():
                h[int(k)] = h.get(int(k), 0) + int(v)
        blk = {k: ins[0].get(k) for k in ('scope', 'definition', 'cap', 'r1_read_len_mode')}
        blk['hist'] = {str(k): v for k, v in sorted(h.items())}
        for k in _INSERT_COUNTS:
            blk[k] = sum(int(b.get(k, 0)) for b in ins)
        out['insert_size'] = insert_summary(blk)
    out['merged_from'] = {'samples': len(parts),
                          'note': 'counts and duplicity histograms summed over the per-sample '
                                  'sidecars (definition D keys on the bobcode, so this equals the '
                                  'pooled deduplication); unique/dup_reads/lib_size are sums of '
                                  'per-sample UMI-only values, not a global clustering'}
    return out


def merge_main(argv):
    """dedup_reads_illumina.py merge-sidecars --out library.json s1.json s2.json ..."""
    import argparse
    ap = argparse.ArgumentParser(prog='dedup_reads_illumina.py merge-sidecars')
    ap.add_argument('parts', nargs='+'); ap.add_argument('--out', required=True)
    a = ap.parse_args(argv)
    out = merge_sidecars([json.load(open(p)) for p in a.parts])
    with open(a.out, 'w') as f:
        json.dump(out, f, indent=1)
    print(f'dedup merge-sidecars: {len(a.parts)} sidecar(s) -> {a.out}: '
          f'reads {out.get("bam_primary", 0):,}, kept {out.get("dedup_kept", 0):,}')


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    force = '--force' in sys.argv
    if not args:
        sys.exit(__doc__)
    run_dir = args[0]
    bcs = args[1:]
    if not bcs:
        fq = os.path.join(run_dir, 'fastq')
        bcs = sorted(d for d in os.listdir(fq) if os.path.isdir(os.path.join(fq, d)))
    os.chdir(run_dir)
    # Import the pipeline modules THROUGH the run's own scripts/ link. star_taxonomy
    # locates the analysisguide relative to its own __file__ (its parent directory), not
    # the cwd — importing from a shared scripts folder directly makes it search there,
    # find no guide, and silently fall back to the WRONG chemistry.
    _run_scripts = os.path.join(os.path.abspath(run_dir), 'scripts')
    if os.path.isdir(_run_scripts):
        sys.path.insert(0, _run_scripts)
    for bc in bcs:
        dedup_barcode(run_dir, bc, force=force)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'merge-sidecars':
        merge_main(sys.argv[2:])
    else:
        main()

#!/usr/bin/env python3
"""Pure-STAR per-read taxonomy + species-barcode swap rate (no minimap2).

For each barcode this reads ONLY STAR-derived / experimental inputs:
  * {bc}_star_combined.bam  -- STARlong alignment to the combined human+mouse
                               genome index (species = chrom prefix; cross-species
                               multimapping = ambiguous).
  * combined_genome.gtf      -- genome annotation (gene_biotype + exons) for the
                               molecular category at each read's locus.
  * rdna_loci.bed            -- data-driven rDNA regions; reads overlapping these
                               are rRNA (the asymmetric-rDNA reads STAR can't
                               species-resolve).
  * {bc} FASTQ               -- the species barcode physically in the read
                               (measure_structure); the swap ground-truth label.

It does NOT use per_read.json species/transcript/genomic_class, read_taxonomy,
the combined_transcriptome, or any minimap2 index -- this is a pure-STAR path.

Writes results['star_taxonomy'] = {category: {human,mouse,ambiguous,na}} (the
composition is over MAPPED reads only -- STAR's unmapped/flag-4 records are not
counted; results['star_taxonomy_scope'] records this) and
results['star_swap'] (per-category + rRNA-excluded specificity) into each
{bc}_speciesmix_results.json.
"""
import os, sys, re, json, bisect, subprocess, math, pickle, hashlib
import multiprocessing as _mp


def _disk_cache(srcpath, tag, builder):
    """Return builder() result, memoized to a pickle on disk keyed on the source
    file's path+mtime. Skips the (slow) GTF parse on reruns when the GTF is
    unchanged. Falls back to builder() on any cache error."""
    try:
        mt = int(os.path.getmtime(srcpath))
        cdir = os.path.join(os.path.dirname(os.path.abspath(srcpath)), '.gtf_cache')
        os.makedirs(cdir, exist_ok=True)
        key = hashlib.md5(f'{srcpath}|{mt}|{tag}'.encode()).hexdigest()[:16]
        cpath = os.path.join(cdir, f'{tag}_{key}.pkl')
        if os.path.exists(cpath):
            with open(cpath, 'rb') as f:
                return pickle.load(f)
        obj = builder()
        tmp = cpath + '.tmp'
        with open(tmp, 'wb') as f:
            pickle.dump(obj, f, protocol=4)
        os.replace(tmp, cpath)
        return obj
    except Exception:
        return builder()


# --- parallel per-barcode enrichment ---------------------------------------
# Workers are fork()ed AFTER the GTF/rdna/terminal/adapter caches are built in
# the parent, so the big read-only structures are shared copy-on-write and never
# pickled. _process_bc_worker reads them from module globals (inherited on fork).
_SHARED_GTF = None
_SHARED_RDNA = None


_CLASSIFY_ONE = None          # set just before forking; inherited by fork()
PARALLEL_MIN_READS = 150_000  # below this, fork + IPC costs more than it saves


def _ox_mod():
    """override_24plex, imported lazily (only the 24plex path needs it)."""
    import override_24plex_illumina as _ox
    return _ox


def _classify_chunk(rids):
    """Worker: category + gene for a slice of read ids. Returns compact tuples
    (rid, category_index, geneinfo) — the index keeps the pickled payload small."""
    fn = _CLASSIFY_ONE
    idx = {c: i for i, c in enumerate(CATEGORY_ORDER)}
    out = []
    for r in rids:
        k, gi = fn(r)
        out.append((r, idx.get(k, -1), gi))
    return out


def _process_bc_worker(bc):
    try:
        process_barcode(bc, _SHARED_GTF, _SHARED_RDNA)
        return (bc, None)
    except Exception as e:
        return (bc, repr(e))
from collections import defaultdict, Counter

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RUN_DIR = os.path.dirname(SCRIPT_DIR)
INDIV = os.path.join(RUN_DIR, 'individual-analyses')
FASTQ = os.path.join(RUN_DIR, 'fastq')
# Reference locations are read from env vars with a literal fallback, so an integrator
# can repoint them for their environment without editing this file.
REF = os.environ.get('BOBSEQ_REF_DIR', 'reference')   # directory holding combined_genome.gtf and rdna_loci.bed; set BOBSEQ_REF_DIR
GTF = os.environ.get('BOBSEQ_GTF', f'{REF}/annotations/combined_genome.gtf')
RDNA_BED = os.environ.get('BOBSEQ_RDNA_BED', f'{REF}/annotations/rdna_loci.bed')

CATEGORY_ORDER = ['mRNA (protein-coding)', 'Ribosomal-protein mRNA',
                  'Other ncRNA (lncRNA/NMD/ret-intron/pseudo)', 'rRNA',
                  'Mitochondrial', 'pre-mRNA / intronic', 'Intergenic / genomic',
                  'Unmapped / contaminant']
SPK = ['human', 'mouse', 'ambiguous', 'na']
EXP = {'p1': 'human', 'p2': 'mouse'}   # polydT species barcode -> species (default orient.)
CODE_SEQ = {'p1': 'ATCGCT', 'p2': 'TCGATA'}  # species barcode -> 6mer sequence (display)
# Molecular categories EXCLUDED from every barcode-correctness metric:
#   rRNA                 -- rDNA is human-only on the combined genome, so rRNA
#                           reads are force-typed human and are unresolvable.
#   Intergenic / genomic -- dominated by insert-less / poly-dT-primed molecules
#                           whose polyT tail maps to genomic poly-A/T tracts in
#                           EITHER species, making the species call arbitrary
#                           (verified by per-read spot-check).
#   Unmapped / contaminant -- no species (already excluded by the species gate).
# Molecular categories excluded from the barcode-correctness denominator (both
# chemistries): rRNA (rDNA asymmetric across the combined genome -> force-typed
# human), Mitochondrial (kept out of correctness by convention), and
# Intergenic/genomic + Unmapped (no confident transcript -> arbitrary species).
CORRECTNESS_EXCLUDE = {'rRNA', 'Mitochondrial', 'Intergenic / genomic', 'Unmapped / contaminant'}
# Low-complexity filter: a read is also dropped from barcode-correctness metrics
# if its STAR-aligned bases (the sequence that determined the species call) have
# a normalized dinucleotide entropy below this cutoff — i.e. the call rests on a
# homopolymer / simple-repeat (e.g. polyT mapping to a genomic poly-A/T tract),
# regardless of read length. Calibrated on real data: entropy<0.65 removes ~70%
# of suspect short-swap reads while dropping ~0.1% of real (long-correct) reads.
LOWCX_ENTROPY_MIN = 0.65
# [^K]*$ : no ribosomal-protein gene name carries a K, but RPS6KA1-6 / RPS6KB1-2 / RPS6KC1 /
# RPS6KL1 (kinases) do and must not be counted as ribosomal-protein mRNA.
_RIBO = re.compile(r'^(RP[LS]|RPLP|MRP[LS]|Rp[ls]|Rplp|Mrp[ls])[^K]*$', re.I)
_RRNA_BT = {'rRNA', 'Mt_rRNA', 'rRNA_pseudogene'}
# PCA gene-exclusion regexes — the correct-vs-swap PCA universe:
# cytoplasmic+mito ribosomal-protein subunits (the [^Kk] guard spares RPS6K*
# kinases) and mitochondrially-encoded genes (mt-* / MT-*).
_PCA_RIBO = re.compile(r'^[Mm]?[Rr][Pp][LlSs][^Kk]*$')
_MITO_RE = re.compile(r'^[Mm][Tt]-')


def _refend(pos, cig):
    rc = 0; n = ''
    for x in cig:
        if x.isdigit(): n += x
        elif x in 'MDN=X': rc += int(n); n = ''
        else: n = ''
    return pos + rc


def _qaln(cig):
    """Query bases aligned (M/I/=/X) -- the STAR 'aligned insert length'."""
    q = 0; n = ''
    for x in cig:
        if x.isdigit(): n += x
        elif x in 'MI=X': q += int(n); n = ''
        else: n = ''
    return q


_NH_RE = re.compile(r'\bNH:i:(\d+)')   # STAR multimapper count
_NM_RE = re.compile(r'nM:i:(\d+)')   # STAR mismatch tag (substitutions, not indels)


def _cigar_mid(cig):
    """(M, I, D) base counts from a CIGAR. M = aligned columns (match+mismatch,
    incl =/X); I = insertion; D = deletion. N (splice)/S/H/P are not errors and
    are ignored. Combined with the nM tag this gives the ONT per-base error."""
    M = I = D = 0; n = ''
    for x in cig:
        if x.isdigit(): n += x
        elif x in 'M=X': M += int(n); n = ''
        elif x == 'I': I += int(n); n = ''
        elif x == 'D': D += int(n); n = ''
        else: n = ''
    return M, I, D


_CIG_OPS = re.compile(r'(\d+)([MIDNSHP=X])')


def _aln_query_bases(seq, cig):
    """The query bases under M/=/X CIGAR ops — i.e. the aligned read sequence
    that determined the alignment (soft-clips and insertions skipped)."""
    i = 0; out = []
    for num, op in _CIG_OPS.findall(cig):
        num = int(num)
        if op in 'M=X':
            out.append(seq[i:i + num]); i += num
        elif op in 'IS':
            i += num
    return ''.join(out)


def _dinuc_entropy(s):
    """Normalized dinucleotide Shannon entropy of `s` in [0,1]: ~0 for a
    homopolymer / simple-repeat, ~1 for complex sequence. Sequences < 3 nt
    return 1.0 (too short to judge on complexity grounds)."""
    n = len(s)
    if n < 3:
        return 1.0
    # Counter(zip(s, s[1:])) counts OVERLAPPING dinucleotides in C (Counter's
    # _count_elements fast path) instead of a Python loop that slices n-1 new
    # 2-char strings. Called ~920k times per 1M reads. Keys are ('A','C') tuples
    # rather than 'AC' strings, which is invisible here: only .values() is used.
    # NB: str.count() cannot replace this -- it skips overlaps, so 'AAA' would give
    # 1 instead of 2, exactly wrong for the homopolymers this flags.
    di = Counter(zip(s, s[1:]))
    tot = n - 1
    h = -sum((c / tot) * math.log2(c / tot) for c in di.values())
    return h / 4.0   # log2(16 possible dinucleotides)


def _median(xs):
    xs = sorted(xs); n = len(xs)
    if not n: return None
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2


def load_rdna(path):
    ivl = defaultdict(list)
    for line in open(path):
        if not line.strip() or line[0] == '#' or line.startswith(('track', 'browser')):
            continue  # skip provenance header / BED meta lines
        c, s, e = line.split()[:3]
        ivl[c].append((int(s), int(e)))
    for c in ivl: ivl[c].sort()
    return ivl, {c: [a for a, _ in v] for c, v in ivl.items()}


def parse_gtf(path):
    return _disk_cache(path, 'gtf', lambda: _parse_gtf_impl(path))


def _parse_gtf_impl(path):
    genes = defaultdict(list); exons = defaultdict(list)
    bt_re = re.compile(r'gene_biotype "([^"]+)"'); nm_re = re.compile(r'gene_name "([^"]+)"')
    with open(path) as f:
        for line in f:
            if '\tgene\t' in line:
                p = line.split('\t')
                genes[p[0]].append((int(p[3]) - 1, int(p[4]),
                                    (bt_re.search(p[8]) or [None, ''])[1],
                                    (nm_re.search(p[8]) or [None, ''])[1]))
            elif '\texon\t' in line:
                p = line.split('\t'); exons[p[0]].append((int(p[3]) - 1, int(p[4])))
    for c in genes: genes[c].sort()
    for c in exons: exons[c].sort()
    return (genes, {c: [g[0] for g in v] for c, v in genes.items()},
            exons, {c: [e[0] for e in v] for c, v in exons.items()})


_TID_RE = re.compile(r'transcript_id "([^"]+)"')
_TERM_CACHE = None


def load_terminal_exons(path):
    return _disk_cache(path, 'termexons', lambda: _load_terminal_exons_impl(path))


def _load_terminal_exons_impl(path):
    """Per-chrom sorted intervals of 3'-terminal exons (the 3'-most exon of each
    transcript: max-end on + strand, min-start on - strand), unioned across all
    transcripts. Returns (intervals, starts) like load_rdna. polydT/3'-capture
    cDNA should align here, so it marks the high-confidence 'real 3' mRNA' set."""
    term = {}   # transcript_id -> (chrom, start, end, strand); keep only 3'-most
    with open(path) as f:
        for line in f:
            if '\texon\t' not in line:
                continue
            p = line.split('\t')
            m = _TID_RE.search(p[8])
            if not m:
                continue
            tid = m.group(1); s = int(p[3]) - 1; e = int(p[4]); strand = p[6]
            cur = term.get(tid)
            if cur is None:
                term[tid] = (p[0], s, e, strand)
            elif strand == '+':
                if e > cur[2]:
                    term[tid] = (p[0], s, e, strand)
            else:
                if s < cur[1]:
                    term[tid] = (p[0], s, e, strand)
    ivl = defaultdict(set)
    for chrom, s, e, _ in term.values():
        ivl[chrom].add((s, e))
    ivl = {c: sorted(v) for c, v in ivl.items()}
    return ivl, {c: [a for a, _ in v] for c, v in ivl.items()}


def _get_terminal():
    """Cached terminal-exon intervals (parsed once per process)."""
    global _TERM_CACHE
    if _TERM_CACHE is None:
        _TERM_CACHE = load_terminal_exons(GTF)
    return _TERM_CACHE


_MAXLEN = {}   # id(ivl) -> {contig: longest interval length}; the exact scan-back bound


def _maxlen(ivl, c):
    d = _MAXLEN.get(id(ivl))
    if d is None:
        d = _MAXLEN[id(ivl)] = {}
    m = d.get(c)
    if m is None:
        m = d[c] = max((seg[1] - seg[0] for seg in ivl[c]), default=0)
    return m


def _overlap(ivl, st, c, s, e, collect=False):
    """Intervals on contig c overlapping [s, e). Intervals are sorted by start; walk back
    from the first start >= e and stop as soon as a start is more than the contig's
    longest interval before s, because nothing earlier can reach s. This is EXACT.
    Stopping after a fixed number of intervals would not be: in dense annotation
    (many short features starting just before a long gene) a genuinely overlapping gene
    could go unexamined and the read would fall to intergenic / intronic."""
    v = ivl.get(c)
    if not v: return [] if collect else False
    i = bisect.bisect_right(st[c], e); out = []
    lo = s - _maxlen(ivl, c)
    for j in range(i - 1, -1, -1):
        seg = v[j]
        if seg[0] < lo: break
        if seg[1] < s: continue
        if seg[0] < e and seg[1] > s:
            if not collect: return True
            out.append(seg)
    return out if collect else False


def classify_bc(bc, gtf, rdna, term):
    genes, gst, exons, est = gtf
    ridna, rdst = rdna
    term_ivl, term_st = term
    bam = os.path.join(INDIV, bc, f'{bc}_star_combined.bam')
    allaln = defaultdict(list); sp = defaultdict(set)
    pv = subprocess.Popen(['samtools', 'view', bam], stdout=subprocess.PIPE, text=True)
    prim = {}; qaln = {}; lowcx = set(); term_reads = set(); multi_nh = set()
    posfp = {}                                 # rid -> (chrom, strand, start, end) positional dup key
    aq_by_read = {}                            # rid -> aligned insert seq (kept for low-complexity gating)
    err = {'M': 0, 'I': 0, 'D': 0, 'sub': 0}   # ONT error tally (primary alns only)
    for line in pv.stdout:
        p = line.split('\t', 6); rid = p[0]; flag = int(p[1])
        if flag & 4: continue
        ch = p[2]; pos = int(p[3]); end = _refend(pos, p[5])
        allaln[rid].append((ch, pos, end))
        if ch.startswith('HUMAN_'): sp[rid].add('human')
        elif ch.startswith('MOUSE_'): sp[rid].add('mouse')
        if not (flag & 0x900):
            prim[rid] = (ch, pos, end); qaln[rid] = _qaln(p[5])
            posfp[rid] = (ch, '-' if (flag & 0x10) else '+', pos, end)   # positional dup key
            _nh = _NH_RE.search(p[6]) if len(p) > 6 else None
            if _nh and int(_nh.group(1)) > 1: multi_nh.add(rid)
            M, I, D = _cigar_mid(p[5])
            mm = _NM_RE.search(p[6]) if len(p) > 6 else None
            err['M'] += M; err['I'] += I; err['D'] += D
            err['sub'] += min(int(mm.group(1)), M) if mm else 0
            # low-complexity flag on the aligned bases (the species-call sequence).
            # SEQ is field 4 of the rest-string left by split('\t', 6).
            rest = p[6].split('\t', 5) if len(p) > 6 else []
            aq = _aln_query_bases(rest[3], p[5]) if len(rest) >= 4 else ''
            if aq and _dinuc_entropy(aq) < LOWCX_ENTROPY_MIN:
                lowcx.add(rid)
            if _overlap(term_ivl, term_st, ch, pos, end):
                term_reads.add(rid)
            if aq:                       # store; the dup fingerprint is built (gated) after categorization
                aq_by_read[rid] = aq

    def species(r):
        s = sp.get(r)
        return 'unmapped' if not s else ('human' if s == {'human'} else 'mouse' if s == {'mouse'} else 'ambiguous')

    def in_rdna(r):
        for c, s, e in allaln[r]:
            v = ridna.get(c)
            if not v: continue
            i = bisect.bisect_right(rdst[c], e)
            for j in range(max(0, i - 1), -1, -1):
                a, b = v[j]
                if b < s: break
                if a < e and b > s: return True
        return False

    def classify_one(r):
        """(category, geneinfo|None). geneinfo = (gene_name, biotype) of the
        max-overlap gene at the primary locus, captured for ALL gene-overlapping
        reads (incl. MT/rRNA) to feed transcript-diversity / PCA gene counts."""
        if r not in prim:
            return 'Unmapped / contaminant', None
        c, s, e = prim[r]
        gs = _overlap(genes, gst, c, s, e, collect=True)
        gi = None
        if gs:
            best = max(gs, key=lambda g: min(g[1], e) - max(g[0], s))
            if best[3]:
                gsp = ('human' if c.startswith('HUMAN_')
                       else 'mouse' if c.startswith('MOUSE_') else None)
                gi = (best[3], best[2], gsp)   # (gene_name, biotype, species)
        if in_rdna(r):
            return 'rRNA', gi
        if c.endswith('_MT'):
            return 'Mitochondrial', gi
        if not gs:
            return 'Intergenic / genomic', None
        if not _overlap(exons, est, c, s, e):
            return 'pre-mRNA / intronic', gi
        bts = [g[2] for g in gs]
        if any(b in _RRNA_BT for b in bts):
            return 'rRNA', gi
        pc = [g for g in gs if g[2].startswith('protein_coding')]
        if pc:
            return ('Ribosomal-protein mRNA' if any(_RIBO.match(g[3] or '') for g in pc)
                    else 'mRNA (protein-coding)'), gi
        return 'Other ncRNA (lncRNA/NMD/ret-intron/pseudo)', gi

    # ---- gene/category assignment: the dominant cost, parallelised ----------
    # classify_one() is a pure per-read lookup into the annotation intervals, so reads
    # are independent. Workers are FORKED, which matters: the GTF/rDNA interval tables
    # are large and fork shares them copy-on-write, so nothing is copied or pickled.
    # Only the read ids go out and (id, category index, gene) come back.
    # Skipped for small inputs (fork + IPC would cost more than it saves) and inside a
    # daemon worker (a multiprocessing daemon may not spawn children).
    cat = {}; geneinfo = {}
    _rids = list(allaln)
    try:
        _daemon = _mp.current_process().daemon
    except Exception:
        _daemon = False
    _nw = min(os.cpu_count() or 1, 8)
    # Module constant (not a literal) so a test can force the serial path and check
    # that parallel classification returns byte-identical categories.
    if len(_rids) >= PARALLEL_MIN_READS and _nw > 1 and not _daemon:
        global _CLASSIFY_ONE
        _CLASSIFY_ONE = classify_one
        _chunks = [_rids[i::_nw] for i in range(_nw)]
        try:
            with _mp.get_context('fork').Pool(_nw) as _pool:
                for _part in _pool.imap_unordered(_classify_chunk, _chunks):
                    for r, ki, gi in _part:
                        cat[r] = CATEGORY_ORDER[ki] if ki >= 0 else CATEGORY_ORDER[-1]
                        if gi:
                            geneinfo[r] = gi
        except Exception as _e:                       # any failure -> serial fallback
            print(f'  (parallel classify failed, falling back to serial: {_e})')
            cat = {}; geneinfo = {}
            for r in _rids:
                k, gi = classify_one(r)
                cat[r] = k
                if gi:
                    geneinfo[r] = gi
    else:
        for r in _rids:
            k, gi = classify_one(r)
            cat[r] = k
            if gi:
                geneinfo[r] = gi
    isrr = {r: (cat[r] == 'rRNA') for r in allaln}
    # PCR duplicates = reads with an IDENTICAL STAR-aligned insert sequence.
    # Reported at TWO scopes: total (all aligned reads) and protein-coding mRNA only.
    # NB: sequence-based (not position dedup) and, without a UMI, an exact-sequence
    # UPPER BOUND -- abundant/structured RNA (esp. the GC-rich rDNA contig) basecalls
    # near-identically across INDEPENDENT molecules and inflates the TOTAL, so the
    # protein-coding-mRNA rate is the informative one.
    # PCR (positional) duplicates: reads sharing the same primary-alignment footprint
    # (chrom, strand, start, end). Position-based, tolerates ONT basecall errors.
    def _dupstats(rids):
        c = Counter(posfp[r] for r in rids if r in posfp)
        return {'exact_dup_reads': sum(v for v in c.values() if v > 1),
                'exact_dup_groups': sum(1 for v in c.values() if v > 1),
                'unique_inserts': len(c), 'n': sum(c.values())}
    star_dup = _dupstats(posfp.keys())                            # total (all aligned)
    star_dup['mrna'] = _dupstats(                                 # protein-coding mRNA only
        [r for r in posfp if cat.get(r) == 'mRNA (protein-coding)'])
    return allaln, sp, species, cat, isrr, qaln, err, geneinfo, lowcx, term_reads, star_dup, multi_nh


def _load_chem_config():
    """(barcode_location, species_map, tso_arch) from the experiment config; defaults
    to the BOB-polydT chemistry if the config can't be read."""
    try:
        sys.path.insert(0, SCRIPT_DIR)
        from analyze_speciesmix_illumina import EXPERIMENT_CONFIG
        return (EXPERIMENT_CONFIG.get('barcode_location', 'BOBcode'),
                EXPERIMENT_CONFIG.get('primary_species_map', {}) or {},
                EXPERIMENT_CONFIG.get('tso_arch'))
    except Exception as e:
        print(f'  (star_taxonomy: config load failed ({e}); defaulting to BOBcode)')
        return 'BOBcode', {}, None


_ADAPTER = None


def _get_adapter():
    """Cached chemistry adapter (BOB-polydT 6mer or TSO 7mer), selected by
    EXPERIMENT_CONFIG['barcode_location']."""
    global _ADAPTER
    if _ADAPTER is None:
        sys.path.insert(0, SCRIPT_DIR)
        import barcode_adapter_illumina as barcode_adapter
        bl, smap, arch = _load_chem_config()
        _ADAPTER = barcode_adapter.get_adapter(bl, smap, arch)
    return _ADAPTER


def barcodes(bc, adapter):
    # ILLUMINA r2_positional: scan the PREPPED (trimmed) FASTQ, because its read names
    # are the ones STAR wrote into the BAM — that is the only key space the BAM can be
    # joined on. Structure comes from the name, not the sequence, so trimming is fine.
    # ONT keeps the untrimmed structure FASTQ.
    import override_24plex_illumina as _ox        # local: only the 24plex path needs it
    _POSITIONAL = getattr(_ox, 'READ_LAYOUT', 'ont_anchor') == 'r2_positional'
    fq = (os.path.join(FASTQ, bc, f'noadapter_combined_{bc}.fastq'))
    if _POSITIONAL:
        for _cand in (os.path.join(FASTQ, bc, f'prepped_{bc}.fastq'),
                      os.path.join(FASTQ, bc, 'prepped_R2.fastq')):
            if os.path.exists(_cand):
                fq = _cand
                break
    code = {}; struct = {}
    el_count = {k: 0 for k, _lbl, _d in adapter.ELEMENTS}   # structural-element tally
    # TSO backbone-completeness breakdown. A 5'-truncated backbone fails the strict
    # full-backbone 'backbone' element, but the barcode anchors on the backbone TAIL
    # next to the 7-mer, so it is still recovered. Report three numbers:
    #   bb_perfect = el_count['backbone'] (full backbone present)
    #   bb_partial = barcode found but WITHOUT a full backbone (truncated)
    #   barcode    = a species barcode was called (any backbone completeness)
    # Only meaningful for the TSO chemistry (has a multiplicity detector).
    _tso_bb = getattr(adapter, 'multiplicity', None) is not None
    if _tso_bb:
        el_count['bb_partial'] = 0
        el_count['barcode'] = 0
    n_reads = 0
    with open(fq) as f:
        while True:
            h = f.readline()
            if not h: break
            s = f.readline().rstrip('\n'); f.readline(); f.readline()
            n_reads += 1
            _full_bb = False
            for k, _lbl, det in adapter.ELEMENTS:
                if det(s):
                    el_count[k] += 1
                    if k == 'backbone':
                        _full_bb = True
            # ILLUMINA r2_positional: the JOIN KEY must match the BAM. prep_reads
            # rewrites the read name to <origid>_<umi1>_<umi2>_<bobcode>_<grun>_<T|F>
            # and that is what STAR aligned, so keying this table on the UNTRIMMED
            # FASTQ's original id would produce a key space that shares NOTHING with the
            # BAM — star_crosstab would come out all zeros and the species-direction
            # sign-off (detect_species_map) would be silently dead. Structure is taken
            # from the prepped read name, which also carries the R1-derived RT-end
            # that adapter.measure(s) cannot see (it gets no mate here).
            if _POSITIONAL:
                m = _ox.measure_from_qname(h[1:])
                if m and m.get('code') in adapter.EXP:
                    rid = h[1:].split()[0]
                    code[rid] = m['code']
                    struct[rid] = m
                    if _tso_bb:
                        el_count['barcode'] += 1
                continue
            m = adapter.measure(s)
            if m and m.get('code') in adapter.EXP:
                rid = h[1:].split()[0]
                code[rid] = m['code']
                struct[rid] = m       # full structure for funnel levels / retention
                if _tso_bb:
                    el_count['barcode'] += 1
                    if not _full_bb:
                        el_count['bb_partial'] += 1
    return code, struct, el_count, n_reads


def _star_crosstabs(code, struct, species, cat, isrr, lowcx, qaln, adapter, multi_tso=frozenset()):
    """Barcode x STAR-species crosstabs (all / MT+rRNA-removed / retained-clean)
    and funnel-level barcode specificity (rRNA+intergenic+low-cx+>=2-TSO chimeras
    excluded). The funnel levels come from the chemistry adapter; the deepest level
    == the 'Full structure' retained set. Pure-STAR."""
    SPK3 = ('human', 'mouse', 'ambiguous')
    # CANONICAL funnel — BARCODE-FIRST so every level is a scorable (barcoded +
    # hu/ms-aligned) population; the reads-passing counts then EQUAL the accuracy
    # denominators. tso-bob fills has7/exact/mrna/insert/polyt/c28 (UMI/TruSeq blank).
    # No C28 step: this construct has no C28 primer, and on the 24plex path the 'c28'
    # flag is never set, so such a step would retain 0 reads and read as a total
    # structural failure rather than an inapplicable filter.
    # The polydT step is DROPPED when the run's R1 is too short to carry a
    # >=10 nt tract (declared r1_polyt=absent, cross-checked against the real R1 length
    # by prep_reads). Keeping it would retain 0 reads and read as a total structural
    # failure rather than an inapplicable filter -- the same reason there is no C28 step.
    import override_24plex_illumina as _ox24
    _POLYT_STEP = bool(getattr(_ox24, 'POLYT_MEASURABLE', True))
    CANON = [('has7', 'has 7-mer barcode'), ('exact', '+ exact 1× barcode'),
             ('mrna', '+ mRNA only'), ('insert', '+ insert >30 bp (STAR)')]
    if _POLYT_STEP:
        CANON.append(('polyt', '+ R1 poly-dT RT end (flag from prep)'))
    _LAST = len(CANON)                            # full-structure keep-set = passes all steps
    def _blank():
        return {seq: {k: 0 for k in SPK3} for seq in adapter.CODE_SEQ.values()}
    def _canon_level(v, k, cd, m, aln, is_lowcx, is_multi):
        """Barcode-first canonical level 0..5: has7 → exact → mRNA → insert → polyT.
        Level 0 = NOT scorable (no barcode / not human-mouse / low-complexity /
        ≥2-TSO chimera). 'mRNA only' = not rRNA/MT/intergenic (this is the accuracy lift)."""
        if cd is None or v not in ('human', 'mouse') or is_lowcx or is_multi:
            return 0                              # not scorable
        if m is None or (m.get('dist') if m.get('dist') is not None else 9) != 0:
            return 1                              # has 7-mer, but not exact
        if not k or k in CORRECTNESS_EXCLUDE:
            return 2                              # exact, but not mRNA
        if aln <= 30:
            return 3                              # mRNA, but insert <30 bp
        # RT end is a FLAG on the positional (Illumina) path: the tract is on R1 and
        # prep records >=10 nt of T after the 26N as rt_end. Legacy anchor-mode
        # measurements still carry a polyt length; accept either form.
        if _POLYT_STEP and not (m.get('rt_end') if m.get('polyt') is None else (m.get('polyt') or 0) >= 10):
            return 4                              # insert, but no RT end / polyT
        return _LAST                              # full structure (polyT step may be absent)
    allx, mtrep, clean = _blank(), _blank(), _blank()
    funnel = {k: {'n': 0, 'co': 0} for k, _ in CANON}   # n == scored (barcode-first)
    _mrna_i = [c[0] for c in CANON].index('mrna')       # reached mRNA <=> lvl > _mrna_i
    sp_mrna = {'human': [0, 0], 'mouse': [0, 0]}        # [scored, correct] at the mRNA level, by ALIGNED species
    retained = set()
    for r in cat:
        v = species(r); k = cat.get(r); aln = qaln.get(r, 0)
        cd = code.get(r); m = struct.get(r)
        lvl = _canon_level(v, k, cd, m, aln, r in lowcx, r in multi_tso)
        ret = (lvl == _LAST)
        if ret:
            retained.add(r)
        # species-mix crosstabs — over ALL coded reads (independent of the funnel)
        if cd is not None:
            seq = adapter.CODE_SEQ[cd]
            if v in SPK3:
                allx[seq][v] += 1
                if k not in ('Mitochondrial', 'rRNA'):
                    mtrep[seq][v] += 1
                    if ret:
                        clean[seq][v] += 1
        # funnel — scorable reads only (level >=1); n and correct move together, so
        # the counts table and the accuracy table share the same denominators.
        if lvl >= 1:
            ok = 1 if v == adapter.EXP[cd] else 0
            for i in range(lvl):
                d = funnel[CANON[i][0]]
                d['n'] += 1; d['co'] += ok
            # Mouse/Human % correct are reported at the mRNA level, split by ALIGNED species.
            if lvl > _mrna_i and v in ('human', 'mouse'):
                sp_mrna[v][0] += 1; sp_mrna[v][1] += ok
    funnel_spec = {k: {'correct': d['co'], 'n': d['n'], 'n_reach': d['n'], 'label': lbl,
                       'specificity': round(d['co'] / d['n'] * 100, 2) if d['n'] else None}
                   for (k, lbl), d in ((c, funnel[c[0]]) for c in CANON)}
    for _asp, _pfx in (('mouse', 'ms'), ('human', 'hu')):   # mRNA-level ms/hu for the Master Summary
        _sc, _co = sp_mrna[_asp]
        funnel_spec[f'_{_pfx}_correct'] = round(100 * _co / _sc, 2) if _sc else None
        funnel_spec[f'_{_pfx}_n'] = _sc
    return allx, mtrep, clean, funnel_spec, retained


def _star_summary(code, retained, species, cat, isrr, qaln, tax, total_reads, lowcx, term_reads, adapter, multi_tso=frozenset()):
    """Per-BC STAR metrics for the cross-barcode master summary. Correctness is
    over confident-species reads, EXCLUDING rRNA / intergenic categories,
    low-complexity alignments (lowcx), and >=2-TSO chimeras (multi_tso: two
    barcodes on one read). `term_reads` = reads whose primary alignment overlaps
    a 3'-terminal exon (high-confidence 3' cDNA)."""
    _codes = list(adapter.EXP.keys())
    hu_code = next((c for c, s in adapter.EXP.items() if s == 'human'), None)
    ms_code = next((c for c, s in adapter.EXP.items() if s == 'mouse'), None)
    by_sp = {'human': {c: 0 for c in _codes}, 'mouse': {c: 0 for c in _codes}}
    hu_len, ms_len, le = [], [], []
    ret_c = ret_t = 0
    term_c = term_t = 0
    for r, cd in code.items():
        v = species(r)
        if (v not in ('human', 'mouse') or cat.get(r) in CORRECTNESS_EXCLUDE
                or r in lowcx or r in multi_tso):
            continue
        ok = 1 if v == adapter.EXP[cd] else 0
        by_sp[v][cd] += 1
        ln = qaln.get(r, 0)
        (hu_len if v == 'human' else ms_len).append(ln)
        if ln > 0:
            le.append((ln, ok))
        if r in retained:
            ret_t += 1; ret_c += ok
        if r in term_reads:
            term_t += 1; term_c += ok
    ms_den = sum(by_sp['mouse'].values())
    hu_den = sum(by_sp['human'].values())
    mrna = tax['mRNA (protein-coding)']
    tot = total_reads or 0
    # per-50bp-bin barcode correctness on the aligned length; drives both the
    # length-effect Spearman and the "% correct vs length" line plot.
    clen = defaultdict(lambda: [0, 0])   # bin index (ln//50) -> [correct, total]
    for ln, ok in le:
        a = clen[ln // 50]; a[1] += 1; a[0] += ok
    length_rho = length_p = None
    if len(le) >= 30:
        bins = sorted(b for b, a in clen.items() if a[1] >= 5)
        if len(bins) >= 3:
            ys = [clen[b][0] / clen[b][1] for b in bins]
            try:
                from scipy.stats import spearmanr
                rho, p = spearmanr(bins, ys)
                if rho == rho:  # not NaN
                    length_rho, length_p = float(rho), float(p)
            except Exception:
                pass
    all_len = hu_len + ms_len
    aln_hist = defaultdict(int)
    for ln in all_len:
        aln_hist[(ln // 50) * 50] += 1          # 50-bp bins for the per-BC line plot
    return {
        'pct_correct_retained': round(ret_c / ret_t * 100, 2) if ret_t else None,
        'terminal_correct': round(term_c / term_t * 100, 2) if term_t else None,
        'ms_spec': (by_sp['mouse'][ms_code] / ms_den) if (ms_den and ms_code) else None,
        'hu_spec': (by_sp['human'][hu_code] / hu_den) if (hu_den and hu_code) else None,
        'hu_mrna_pct': (mrna['human'] / tot * 100) if tot else None,
        'ms_mrna_pct': (mrna['mouse'] / tot * 100) if tot else None,
        'hu_aln': _median(hu_len), 'ms_aln': _median(ms_len),
        'aln_median': _median(all_len),
        # means over the SAME length lists as the medians above (so the master
        # summary can show "mean / median" from a consistent read set).
        'hu_mean': (round(sum(hu_len) / len(hu_len)) if hu_len else None),
        'ms_mean': (round(sum(ms_len) / len(ms_len)) if ms_len else None),
        'aln_mean': (round(sum(all_len) / len(all_len)) if all_len else None),
        'aln_hist': {str(k): v for k, v in sorted(aln_hist.items())},
        'length_rho': length_rho, 'length_p': length_p,
        'correct_vs_len': {str(b * 50): v for b, v in sorted(clen.items())},
    }


def _star_gene_stats(allaln, species, code, geneinfo, adapter, isrr=None, multi_nh=None):
    """Returns (diversity, top_genes, pca_counts) from STAR gene assignments.

    rRNA MASK: reads flagged rRNA (isrr — they overlap a `rdna_loci.bed` locus) are
    EXCLUDED from gene-level counting. Genomes carry degenerate rDNA fragments inside
    protein-coding-gene introns (Filip1l/HFM1…); abundant rRNA maps there and, since
    geneinfo is the max-overlap gene regardless of the rRNA call, would otherwise be
    miscounted as that gene's mRNA (dominating top-gene lists, ~2x-inflating mRNA
    reads in rRNA-heavy libraries). Matches the composition table's rRNA call.

      diversity  -- {total_mapped, unique_genes, simpson, top1/5/10_pct} over
                    PROTEIN-CODING mRNAs only (ribosomal-protein + mito genes
                    excluded). 'total_mapped' = mRNA reads. Keyed by gene name.
      top_genes  -- {'human': {'total': N, 'top': [[gene, count], ...]},
                     'mouse': {...}}  top protein-coding mRNAs split by species.
      pca_counts -- {species: {'correct'/'swap': {gene: count}}} over the same
                    protein-coding, non-ribo, non-mito universe (PCA features).
    """
    # Restrict diversity to protein-coding mRNAs (non-ribosomal-protein,
    # non-mitochondrial) — the biologically informative universe; rRNA/mt-rRNA
    # would otherwise dominate every top-gene list and the concentration stats.
    _rr = isrr or {}
    n_rrna_masked = 0
    mrna_total = Counter()
    mrna_sp = {'human': Counter(), 'mouse': Counter()}   # split by gene species
    multi_sp = {'human': Counter(), 'mouse': Counter()}  # of those, STAR multimappers (NH>1)
    other_sp = Counter(); other_multi = Counter()        # genes on a contig of neither species
    _mn = multi_nh or set()
    for r in allaln:
        gi = geneinfo.get(r)
        if not gi or not gi[0] or gi[1] != 'protein_coding':
            continue
        if _rr.get(r):                       # rRNA on an rDNA-in-gene locus — not real mRNA
            n_rrna_masked += 1
            continue
        name = gi[0]
        if _PCA_RIBO.match(name) or _MITO_RE.match(name):
            continue
        mrna_total[name] += 1
        gsp = gi[2] if len(gi) > 2 else None
        if gsp in mrna_sp:
            mrna_sp[gsp][name] += 1
            if r in _mn: multi_sp[gsp][name] += 1
        else:
            other_sp[name] += 1
            if r in _mn: other_multi[name] += 1
    total = sum(mrna_total.values())
    top_genes = {sp: {'total': sum(c.values()),
                      # third field = fraction of the gene's reads that are STAR multimappers
                      # (NH>1): >= 0.9 marks a multi-copy / repeat-like locus, not expression
                      'top': [[g, n, round(multi_sp[sp][g] / n, 2)] for g, n in c.most_common(7)]}
                 for sp, c in mrna_sp.items()}
    # Complete per-gene table, so a library-level unit can be merged from per-sample units
    # with exact unique-gene counts, top-N shares and multimapper fractions (per-sample
    # mode). {gene species: {gene: [reads, multimapper reads]}}; 'other' = neither prefix.
    top_genes['by_gene'] = {sp: {g: [n, multi_sp[sp][g]] for g, n in c.items()} for sp, c in mrna_sp.items()}
    top_genes['by_gene']['other'] = {g: [n, other_multi[g]] for g, n in other_sp.items()}
    diversity = None
    if total:
        counts = sorted(mrna_total.values(), reverse=True)
        simpson = 1.0 - sum((c / total) ** 2 for c in counts)   # Gini–Simpson
        def topk(k):
            return round(sum(counts[:k]) / total * 100, 1)
        # Per-species unique gene counts. NOT derivable from unique_genes: that key is
        # over gene NAMES pooled across species, and human/mouse orthologues often share
        # a name up to case (ACTB / Actb) — mrna_total keys on the name alone, so pooling
        # under-counts. mrna_sp is keyed by the gene's own species (chrom prefix), so
        # these two can legitimately sum to MORE than unique_genes.
        diversity = {'total_mapped': total, 'unique_genes': len(mrna_total),
                     'unique_genes_human': len(mrna_sp['human']),
                     'unique_genes_mouse': len(mrna_sp['mouse']),
                     'rrna_masked': n_rrna_masked,   # reads dropped as rDNA-in-gene rRNA
                     'simpson': round(simpson, 4),
                     'top1_pct': topk(1), 'top5_pct': topk(5), 'top10_pct': topk(10)}
    pca = {'mouse': {'correct': Counter(), 'swap': Counter()},
           'human': {'correct': Counter(), 'swap': Counter()}}
    for r, cd in code.items():
        v = species(r)
        if v not in ('mouse', 'human'):
            continue
        gi = geneinfo.get(r)
        if not gi or gi[1] != 'protein_coding':
            continue
        if _rr.get(r):                       # same rRNA mask for the PCA feature counts
            continue
        name = gi[0]
        if not name or _PCA_RIBO.match(name) or _MITO_RE.match(name):
            continue
        kind = 'correct' if v == adapter.EXP[cd] else 'swap'
        pca[v][kind][name] += 1
    pca_counts = {sp_: {k: dict(c) for k, c in d.items()} for sp_, d in pca.items()}
    return diversity, top_genes, pca_counts


def _tso_multiplicity_stats(bc, adapter, species, qaln, struct):
    """TSO-unit multiplicity + barcode concordance for the TSO-structure page
    (pure-STAR; runs only when adapter.multiplicity is set, i.e. TSO chemistry).
    Re-scans the FASTQ, counts TSO units per read (1 = canonical, >=2 = TSO-TSO
    chimera), and for 2-TSO reads records the barcode pair (same- vs different-
    species), the orientation, the random-pairing expectation, and the
    productivity (% aligned / retained / species-called) of >=2-TSO vs 1-TSO
    reads."""
    fq = os.path.join(FASTQ, bc, f'noadapter_combined_{bc}.fastq')
    n_dist = Counter(); unit_codes = Counter()
    unit_seqs = Counter(); pair_seq = Counter()    # by CODE (any two samples)
    pair = Counter(); orient = Counter()
    one_rids = []; multi_rids = []
    n_reads = 0
    with open(fq) as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline().rstrip('\n'); f.readline(); f.readline()
            n_reads += 1
            m = adapter.multiplicity(s)
            n = m['n_tso']
            n_dist[min(n, 3)] += 1
            for c in m['codes']:
                unit_codes[c] += 1
            for cs in (m.get('code_seqs') or []):
                unit_seqs[cs] += 1
            rid = h[1:].split()[0]
            if n == 1:
                one_rids.append(rid)
            elif n >= 2:
                multi_rids.append(rid)
                if n == 2:
                    pair[tuple(sorted(m['codes']))] += 1
                    pair_seq['diff' if len(set(m.get('code_seqs') or [])) > 1 else 'same'] += 1
                    orient[''.join(m['orientations'])] += 1

    def _prod(rids):
        n = len(rids)
        if not n:
            return {'n': 0, 'aligned_pct': None, 'retained_pct': None, 'spcall_pct': None}
        na = sum(1 for r in rids if qaln.get(r, 0) > 0)
        nr = sum(1 for r in rids
                 if r in struct and adapter.retained_reason(struct[r], qaln.get(r, 0))[0])
        ns = sum(1 for r in rids if species(r) in ('human', 'mouse'))
        return {'n': n, 'aligned_pct': round(100 * na / n, 1),
                'retained_pct': round(100 * nr / n, 1), 'spcall_pct': round(100 * ns / n, 1)}

    n2 = sum(pair.values())
    diff = pair.get(('h', 'm'), 0)
    same = pair.get(('h', 'h'), 0) + pair.get(('m', 'm'), 0)
    th = unit_codes.get('h', 0); tm = unit_codes.get('m', 0); tu = th + tm or 1
    exp_diff = 2.0 * (th / tu) * (tm / tu)   # P(two independent draws differ)
    ge2 = len(multi_rids)
    stats = {
        'n_total': n_reads,
        'count_dist': {'0': n_dist[0], '1': n_dist[1], '2': n_dist[2], '3+': n_dist[3]},
        'n_ge2': ge2,
        'ge2_pct': round(100 * ge2 / n_reads, 1) if n_reads else None,
        'pair2': {'hh': pair.get(('h', 'h'), 0), 'mm': pair.get(('m', 'm'), 0), 'hm': diff},
        'same2': same, 'diff2': diff,
        'same2_pct': round(100 * same / n2, 1) if n2 else None,
        'diff2_pct': round(100 * diff / n2, 1) if n2 else None,
        'exp_diff2_pct': round(100 * exp_diff, 1),
        # Pairing by CODE: works on any run with >= 2 codes (a 24-plex screen has 21 human
        # codes, so 'same species' hides most chimera pairs). Random expectation for two
        # independent draws being different codes = 1 - sum(p_i^2).
        'same_code2': pair_seq.get('same', 0), 'diff_code2': pair_seq.get('diff', 0),
        'same_code2_pct': round(100 * pair_seq.get('same', 0) / n2, 1) if n2 else None,
        'diff_code2_pct': round(100 * pair_seq.get('diff', 0) / n2, 1) if n2 else None,
        'exp_diff_code2_pct': round(100 * (1.0 - sum((v / max(sum(unit_seqs.values()), 1)) ** 2 for v in unit_seqs.values())), 1) if unit_seqs else None,
        'species_split_informative': (min(th, tm) / tu) >= 0.10,
        'unit_code_seqs': dict(unit_seqs),
        'orient2': dict(orient),
        'prod_ge2': _prod(multi_rids),
        'prod_one': _prod(one_rids),
        'unit_codes': {'h': th, 'm': tm},
    }
    # (stats, exclusion-set) — the >=2-TSO rids are dropped from every barcode-
    # correctness metric (ambiguous: two barcodes on one read).
    return stats, set(multi_rids)


def process_barcode(bc, gtf, rdna):
    adapter = _get_adapter()
    term = _get_terminal()
    allaln, sp, species, cat, isrr, qaln, err, geneinfo, lowcx, term_reads, star_dup, multi_nh = classify_bc(bc, gtf, rdna, term)
    code, struct, el_count, n_reads = barcodes(bc, adapter)
    # TSO-TSO chimera exclusion: reads carrying >=2 TSO units have two barcodes
    # (96% cross-species) and are ambiguous -> dropped from EVERY barcode-
    # correctness metric below (analogous to lowcx). Empty set for chemistries
    # without a multiplicity detector (e.g. polydT), so this is a no-op there.
    tso_mult = None
    multi_tso = set()
    # ILLUMINA r2_positional: DISABLED, on both correctness and cost grounds.
    #
    # CORRECTNESS. A TSO-TSO chimera is defined by seeing TWO TSO units inside ONE
    # read. A 150 bp mate sees one end of the molecule, so the second unit is not
    # observable -- the PDF already drops the '% reads with >=2 TSO (chimera)' row as
    # "not measurable" for exactly this reason. What the positional detector actually
    # finds is chance: it slides a <=1-mismatch 7-mer match across ~140 windows, and
    # 22/16384 patterns match per code per window. Measured on 100k reads of
    # a 2x150 library: real reads 1.55% >=2 units vs 1.89% on a control where
    # everything after the true bobcode is SHUFFLED (chimeras impossible). The
    # control rate is HIGHER, i.e. zero signal -- yet `multi_tso` was removing those
    # ~16k reads from the barcode-accuracy denominator.
    #
    # COST. That window scan was 193.6s of the 239.5s barcode stage (81%) at 1M reads.
    import override_24plex_illumina as _ox        # local: only the 24plex path needs it
    _positional = getattr(_ox, 'READ_LAYOUT', 'ont_anchor') == 'r2_positional'
    if getattr(adapter, 'multiplicity', None) and not _positional:
        tso_mult, multi_tso = _tso_multiplicity_stats(bc, adapter, species, qaln, struct)
    # composition
    tax = {k: {s: 0 for s in SPK} for k in CATEGORY_ORDER}
    for r in allaln:
        tax[cat[r]][species(r)] += 1
    # swap: per category + overall rRNA-excluded (chimeras excluded)
    percat = {}
    ov_c = ov_t = 0
    for r in code:
        v = species(r)
        if v not in ('human', 'mouse'): continue
        ok = (v == adapter.EXP[code[r]])
        k = cat[r]; d = percat.setdefault(k, [0, 0]); d[1] += 1; d[0] += ok
        # Overall (the headline barcode-accuracy number) drops rRNA/intergenic
        # categories, low-complexity alignments, and >=2-TSO chimeras; the raw
        # per-category rows above are left unfiltered (diagnostic), matching lowcx.
        if k not in CORRECTNESS_EXCLUDE and r not in lowcx and r not in multi_tso:
            ov_t += 1; ov_c += ok
    swap = {'per_category': {k: {'correct': v[0], 'n': v[1],
                                 'specificity': round(v[0] / v[1] * 100, 2) if v[1] else None}
                             for k, v in percat.items()},
            'overall_rrna_excluded': {'correct': ov_c, 'n': ov_t,
                                      'specificity': round(ov_c / ov_t * 100, 2) if ov_t else None}}
    # write into results JSON
    rp = os.path.join(INDIV, bc, f'{bc}_speciesmix_results.json')
    res = json.load(open(rp))
    res['star_taxonomy'] = tax
    res['star_duplicates'] = star_dup
    res['star_taxonomy_order'] = CATEGORY_ORDER
    res['star_taxonomy_scope'] = 'mapped reads'  # composition excludes unmapped/flag-4 records
    res['star_swap'] = swap
    crossall, crossmtrep, crossclean, funnel_spec, retained = _star_crosstabs(
        code, struct, species, cat, isrr, lowcx, qaln, adapter, multi_tso)
    # Barcode accuracy vs non-templated G-run length, mRNA (protein-coding) reads only.
    # Reads reaching here are already UMI-deduplicated (dedup_reads_illumina runs right
    # after STAR), so this is one vote per MOLECULE — which matters: PCR copies of a
    # mis-trimmed read all carry the same G-run and the same barcode error, so counting
    # reads dilutes the effect. Measured on the 1M subset, deduplicating widened the
    # G-run 0 deficit (98.60% -> 98.13%) while G-run >=3 rose to ~99.94%.
    # mRNA-only on purpose: on unfiltered reads the high-G-run bins are rRNA/intergenic
    # enriched, which manufactures a fake accuracy decline that is pure composition.
    _ga = defaultdict(lambda: [0, 0])
    for _r, _bcd in code.items():
        _v = species(_r)
        if _v not in ('human', 'mouse') or cat.get(_r) != 'mRNA (protein-coding)':
            continue
        _mq = _ox_mod().measure_from_qname(_r)
        if _mq is None or _mq.get('grun') is None:
            continue
        _d = _ga[_mq['grun']]
        _d[1] += 1
        _d[0] += (_v == adapter.EXP[_bcd])
    res['star_grun_accuracy'] = {str(g): v for g, v in sorted(_ga.items())}
    res['star_crosstab'] = crossall
    res['star_crosstab_mtrep'] = crossmtrep
    res['star_crosstab_clean'] = crossclean
    res['star_funnel_spec'] = funnel_spec
    # detected per-element length distributions (the chemistry adapter declares
    # which numeric structure fields to histogram: TSO/C28 for polydT,
    # backbone/polyT for TSO). Keyed by field name.
    _hists = {field: defaultdict(int) for field, _lbl in adapter.LEN_HIST}
    for m in struct.values():
        for field, _lbl in adapter.LEN_HIST:
            _v = m.get(field)
            if _v is None:
                continue                  # not observable on this read layout
            _hists[field][int(_v)] += 1
    res['star_len_hist'] = {field: {str(k): v for k, v in sorted(h.items())}
                            for field, h in _hists.items()}
    res['star_len_hist_labels'] = {field: lbl for field, lbl in adapter.LEN_HIST}
    summ = _star_summary(code, retained, species, cat, isrr, qaln, tax,
                         res.get('total_reads', 0), lowcx, term_reads, adapter, multi_tso)
    summ['n_lowcx_excluded'] = len(lowcx)
    summ['terminal_pct'] = round(len(term_reads) / len(allaln) * 100, 1) if allaln else None
    summ['terminal_n'] = len(term_reads); summ['aln_n'] = len(allaln)   # counts behind the % (library merge)
    summ['pct_correct'] = swap['overall_rrna_excluded']['specificity']
    # ONT per-base error over the aligned read cores (substitutions from STAR's
    # nM tag, indels from CIGAR I/D; N/soft-clips excluded). This is a lower
    # bound on the raw error since STAR soft-clips error-dense read ends.
    _aln = err['M'] + err['I'] + err['D']
    summ['err_pct'] = round((err['sub'] + err['I'] + err['D']) / _aln * 100, 2) if _aln else None
    summ['err_sub_pct'] = round(err['sub'] / err['M'] * 100, 2) if err['M'] else None
    summ['err_indel_pct'] = round((err['I'] + err['D']) / _aln * 100, 2) if _aln else None
    summ['err_aln_bases'] = _aln
    # chemistry-aware structural-element rates (% of all reads), for the Master
    # Summary "% reads with …" rows; empty for chemistries with no ELEMENTS.
    summ['struct_pct'] = {k: (round(el_count[k] / n_reads * 100, 1) if n_reads else None)
                          for k in el_count}
    summ['struct_labels'] = {k: lbl for k, lbl, _d in adapter.ELEMENTS}
    res['star_summary'] = summ
    # TSO-structure multiplicity (computed above to build the chimera-exclusion
    # set). 1 TSO = canonical, >=2 = TSO-TSO chimera; barcode concordance among
    # 2-TSO reads. (summ is the same object as res['star_summary'].)
    if tso_mult is not None:
        res['tso_multiplicity'] = tso_mult
        summ['tso_ge2_pct'] = tso_mult.get('ge2_pct')
    # transcript diversity + top genes + correct-vs-swap PCA gene counts (STAR)
    diversity, top_genes, pca_counts = _star_gene_stats(allaln, species, code, geneinfo, adapter, isrr, multi_nh)
    res['star_transcript_diversity'] = diversity
    res['star_gene_reads'] = top_genes.pop('by_gene', {})   # complete per-gene table (library merge)
    res['star_top_genes'] = top_genes
    res['star_pca_gene_counts'] = pca_counts
    # read-structure retention for §1 (figure + funnel-counts table), chemistry-
    # aware and pure-STAR: FASTQ structure (struct) + STAR aligned length (qaln)
    # as the insert; drop reasons / funnel from the adapter. No per_read.json.
    _drop = Counter(); _fcounts = {k: 0 for k, _l, _g in adapter.FUNNEL}; _pthist = Counter()
    # parallel aggregation with rRNA reads excluded — for the §1 "rRNA-excluded"
    # diagnostic figure (TSO mis-priming on rRNA dominates the short-polyT peak;
    # subtracting it shows the real cDNA-capture quality).
    _drop_nr = Counter(); _fcounts_nr = {k: 0 for k, _l, _g in adapter.FUNNEL}; _pthist_nr = Counter()
    n_struct_rrna   = sum(1 for r in struct if cat.get(r) == 'rRNA')
    n_aligned_rrna  = sum(1 for r in allaln if cat.get(r) == 'rRNA')
    n_unstruct_rrna = n_aligned_rrna - n_struct_rrna
    n_no_bc         = max(n_reads - len(struct), 0)
    _drop[adapter.NO_CODE_REASON]    += n_no_bc
    _drop_nr[adapter.NO_CODE_REASON] += max(n_no_bc - n_unstruct_rrna, 0)
    _last_lvl = len(adapter.FUNNEL) - 1
    for rid, m in struct.items():
        aln = qaln.get(rid, 0)
        _kept, _reason = adapter.retained_reason(m, aln)
        _lvl = adapter.funnel_level(m, aln)
        # a read with >=2 TSO units (two 7-mer barcodes) is a chimera and is
        # disqualified regardless of where else it fails. ALL such reads go in
        # the 'chimera' drop bucket (so the §1 figure agrees with the §-Multiplicity
        # chimera %); the funnel cumulative count is capped below 'retained'.
        if rid in multi_tso:
            _reason = 'chimera'
            if _lvl == _last_lvl:
                _lvl = _last_lvl - 1
        _drop[_reason] += 1
        if m.get('polyt') is not None:
            _pthist[int(m.get('polyt'))] += 1
        for _i in range(_lvl + 1):
            _fcounts[adapter.FUNNEL[_i][0]] += 1
        if cat.get(rid) != 'rRNA':
            _drop_nr[_reason] += 1
            if m.get('polyt') is not None:
                _pthist_nr[int(m.get('polyt'))] += 1
            for _i in range(_lvl + 1):
                _fcounts_nr[adapter.FUNNEL[_i][0]] += 1
    res[adapter.RET_KEY] = {
        'n_total': n_reads,
        'drop_reason': dict(_drop),
        'funnel': _fcounts,
        'funnel_labels': {k: lbl for k, lbl, _g in adapter.FUNNEL},
        'polyt_hist': {str(k): v for k, v in sorted(_pthist.items())},
    }
    res[adapter.RET_KEY + '_no_rrna'] = {
        'n_total': n_reads - n_aligned_rrna,
        'drop_reason': dict(_drop_nr),
        'funnel': _fcounts_nr,
        'funnel_labels': {k: lbl for k, lbl, _g in adapter.FUNNEL},
        'polyt_hist': {str(k): v for k, v in sorted(_pthist_nr.items())},
    }
    json.dump(res, open(rp, 'w'), indent=2)
    return tax, swap


def main():
    print('parsing combined GTF (genome annotation)...')
    gtf = parse_gtf(GTF)
    rdna = load_rdna(RDNA_BED)
    bcs = sorted(d for d in os.listdir(INDIV)
                 if d.startswith('barcode')
                 and os.path.exists(os.path.join(INDIV, d, f'{d}_star_combined.bam')))
    print(f'{len(bcs)} barcodes with STAR BAMs\n')
    print(f"{'barcode':<11}{'mapped':>8}{'rRNA%':>7}{'mRNA%':>7}{'MT%':>6}"
          f"{'spec(rRNA-excl)':>17}")
    for bc in bcs:
        tax, swap = process_barcode(bc, gtf, rdna)
        tot = sum(sum(v.values()) for v in tax.values())
        rr = sum(tax['rRNA'].values()) / tot * 100
        mr = sum(tax['mRNA (protein-coding)'].values()) / tot * 100
        mt = sum(tax['Mitochondrial'].values()) / tot * 100
        sp = swap['overall_rrna_excluded']['specificity']
        print(f"{bc:<11}{tot:>8,}{rr:>6.1f}%{mr:>6.1f}%{mt:>5.1f}%{str(sp)+'%':>17}")
    print('\nWrote star_taxonomy + star_swap into each *_speciesmix_results.json')


if __name__ == '__main__':
    main()

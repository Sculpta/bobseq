#!/usr/bin/env python3
"""Full ordered-molecule funnel for the 24plex chemistry (§1 read-structure figure).

Counts reads carrying the canonical molecule elements, added in this order
(insert-out — the payload first, structural ends last):

  1. insert >30 bp     STAR aligned query length > 30 bp. House style (structure
                       from the FASTQ, insert length from the STAR alignment); as
                       the first filter it removes primer-dimers, which do not
                       align (their "insert" is ~2 bp).
  2. + polydT (>10)    a poly-T/A run >=10 (the RT end).
  3. + 26N UMI         an 18-34 nt spacer sits 5' of the polyT (room for the 26N
                       RT-primer UMI), i.e. the polyT does not butt onto an adapter.
  4. + TruSeq1         the Read1 handle (CTACACGACG..CTCTTCCGATCT) ends 18-34 nt 5'
                       of the polyT, i.e. that 26N gap is bounded by TruSeq1.
  5. + TruSeq2 = FULL  measure_struct_24plex finds the TSO backbone + species 7-mer.

Each level is cumulative (a subset of the one above); the deepest level a read
reaches is its class, and step n = reads reaching AT LEAST level n. Reads come
from the STAR BAM (read IDs + CIGAR for the insert + SEQ for the structure), so
no separate FASTQ pass is needed. 24plex-specific (TruSeq1-26N-polydT RT end)."""
import os
import re
import bisect
import subprocess

import override_24plex_illumina as ox
import star_taxonomy_illumina as st

# Opt-in analysis hook (off unless BOBSEQ_DUMP_DISCORDANT names a file): appends the
# read ids this funnel scores as discordant. This is the funnel behind the published
# 'Mouse/Human % correct' rows, so these ids ARE that population -- not a
# reconstruction of it. No effect on any published number.
_DUMP_DISCORDANT = os.environ.get('BOBSEQ_DUMP_DISCORDANT', '')
_DUMP_FH = open(_DUMP_DISCORDANT, 'a') if _DUMP_DISCORDANT else None


def _rc(s):
    return s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


ANCHOR = 'CTCTTCCGATCT'       # shared R1 3' handle (3' end of both TruSeq1 and TruSeq2)
TS1_D = 'CTACACGACG'          # distinctive 5' of TruSeq1 (Read1 / RT primer)
TS2_D = 'AGACGTGTG'           # distinctive 5' of TruSeq2 (TSO backbone)
# "partial" handle = the shared anchor + at least half of the distinctive 5' bases.
# ONT clips the outer ends, so requiring the full handle wrongly fails clipped reads;
# half the distinctive part (+ the 12 bp anchor) still identifies WHICH handle it is.
_TS1_MINMATCH = (len(TS1_D) + 1) // 2          # >= 5 of 10
_TS2_MINMATCH = (len(TS2_D) + 1) // 2          # >= 5 of 9
_POLY = re.compile(r'T{10,}')
_CIG = re.compile(r'(\d+)([MIDNSHP=X])')
_UMI_LO, _UMI_HI = 18, 40     # 26N tolerance window (widened slightly for ONT indels)

# barcode->insert distance step (inserted right after 'mRNA only'). The STAR-aligned
# insert must sit ~(UMI+spacer + GGG) nt 3' of the barcode; reads whose insert is on
# the 5' side of the barcode (wrong orientation) fail, reads whose alignment overlaps
# the barcode ('straddle', an alignment artifact) are kept. GRUN_OFFSET is the per-
# barcode UMI+spacer offset, set from the run config by analyze_speciesmix (empty ->
# a default of 8 for every barcode).
GRUN_OFFSET = {}
GGG_LEN = 3
DIST_TOL = 25                 # |observed gap - expected| tolerance (nt)

# Canonical funnel — BARCODE-FIRST so every level is a scorable (barcoded +
# species-aligned) population and the funnel counts == the accuracy denominators.
# Shared with tso-bob (which fills has7/exact/mrna/insert/polyt/c28, leaves
# umi/ts1/ts2 blank; 24plex fills all but c28).
STEPS = [
    ('has7',    'has one barcode'),
    ('mrna',    '+ mRNA only'),
    ('dist',    '+ G-run present (distance enforced at prep)'),
    ('grun3',   '+ G-run ≥3'),
    ('insert',  '+ insert >30 bp (STAR) = full structure'),
]
# polyT / 26N UMI / TruSeq1 / TruSeq2 are NOT funnel rows. On Illumina those four rows
# would be ONE observation (an R1 polyT tract) reported four times, and three of them could
# not fail: the UMI row would restate the same rt_end test (which already requires the polyT
# to start at 18-40, i.e. the 26N spacer), and TruSeq1/TruSeq2 are the SEQUENCING
# PRIMERS -- absent them the read would never have been sequenced, or would have been
# dropped at prep for lacking a bobcode. All three would read 100% by construction, which is
# worse than a row that reads 0%: a failing row is visibly broken, a row pinned at
# 100% looks like verified data.


def _aligned_qlen(cig):
    if cig == '*':
        return 0
    return sum(int(n) for n, o in _CIG.findall(cig) if o in 'MI=X')


def _refend(pos, cig):
    """1-based reference end of the alignment."""
    return pos + sum(int(n) for n, o in _CIG.findall(cig) if o in 'MDN=X') - 1


def _in_rdna(ridna, rdst, c, s, e):
    v = ridna.get(c)
    if not v:
        return False
    i = bisect.bisect_right(rdst[c], e)
    for j in range(max(0, i - 1), -1, -1):
        a, b = v[j]
        if b < s:
            break
        if a < e and b > s:
            return True
    return False


def _rrna_mt_cat(ch, s, e, gtf, rdna):
    """'rrna' / 'mt' / None for the primary alignment, matching
    star_taxonomy.classify_one's decision order: rDNA overlap -> rRNA; then _MT
    chromosome -> Mitochondrial; then rRNA gene-biotype on an exonic overlap ->
    rRNA (rRNA and MT are mutually exclusive categories)."""
    genes, gst, exons, est = gtf
    ridna, rdst = rdna
    if _in_rdna(ridna, rdst, ch, s, e):
        return 'rrna'                                   # rRNA (rDNA locus)
    if ch.endswith('_MT'):
        return 'mt'                                     # Mitochondrial
    gs = st._overlap(genes, gst, ch, s, e, collect=True)
    if gs and st._overlap(exons, est, ch, s, e) and any(g[2] in st._RRNA_BT for g in gs):
        return 'rrna'                                   # rRNA by gene biotype
    return None


def _dist_match(pre, distinctive):
    """# of matching bases of `distinctive` against the 3'-aligned tail of `pre`."""
    t = pre[-len(distinctive):]
    if len(t) < len(distinctive):
        return -1
    return sum(1 for a, b in zip(t, distinctive) if a == b)


def _ts2_partial(seq):
    """A TruSeq2 handle present anywhere (either strand): the shared anchor with >=
    half of the AGACGTGTG distinctive 5' bases immediately before it."""
    for s in (seq, _rc(seq)):
        st = 0
        while True:
            j = s.find(ANCHOR, st)
            if j < 0:
                break
            st = j + 1
            if _dist_match(s[:j], TS2_D) >= _TS2_MINMATCH:
                return True
    return False


def _ts1_partial_before(s, ps):
    """A (partial) TruSeq1 handle ends 18-40 nt 5' of the polyT start ps: the shared
    anchor with >= half of the CTACACGACG distinctive 5' bases before it."""
    lo = max(0, ps - (_UMI_HI + len(ANCHOR) + 4))
    win = s[lo:ps]
    for m in re.finditer(ANCHOR, win):
        e = lo + m.end()
        if _UMI_LO <= ps - e <= _UMI_HI and _dist_match(s[:lo + m.start()], TS1_D) >= _TS1_MINMATCH:
            return True
    return False


def _rt_reach(seq):
    """Deepest RT-side element reached: 0 none, 1 polyT, 2 +26N UMI, 3 +TruSeq1."""
    best = 0
    for s in (seq, _rc(seq)):
        for m in _POLY.finditer(s):
            ps = m.start()
            best = max(best, 1)                 # polyT
            if ps >= _UMI_LO:                   # room for a 26N spacer 5' of the polyT
                best = max(best, 2)             # + 26N UMI
                if _ts1_partial_before(s, ps):
                    best = max(best, 3)         # + TruSeq1 (partial) bounds the 26N
    return best


def _cat_full(ch, s, e, gtf, rdna):
    """Primary-alignment category: 'rrna' / 'mt' / 'mrna' (protein-coding) /
    'other' (intergenic / non-coding), matching star_taxonomy's decision order
    (rDNA → MT chromosome → rRNA biotype → protein-coding biotype)."""
    genes, gst, exons, est = gtf
    ridna, rdst = rdna
    if _in_rdna(ridna, rdst, ch, s, e):
        return 'rrna'
    if ch.endswith('_MT'):
        return 'mt'
    gs = st._overlap(genes, gst, ch, s, e, collect=True)
    ex = st._overlap(exons, est, ch, s, e)
    if gs and ex and any(g[2] in st._RRNA_BT for g in gs):
        return 'rrna'
    if gs and ex and any(str(g[2]).startswith('protein_coding') for g in gs):
        return 'mrna'
    return 'other'


_OFFSET_DEFAULTED = set()


def _default_offset(code):
    """Offset for a code the run config did not declare: 8 (the +8 human spacer),
    reported ONCE per code on stderr instead of silently."""
    if code not in _OFFSET_DEFAULTED:
        _OFFSET_DEFAULTED.add(code)
        import sys as _sys
        _sys.stderr.write(f'full_mol_funnel: WARNING no grun_offset declared for {code}; '
                          f'using 8. Declare it in the run config (grun_offsets).\n')
    return 8


def _dist_ok(seq, cig, m):
    """barcode->insert distance test (24plex): the STAR-aligned insert must sit
    ~(UMI+spacer+GGG) nt 3' of the barcode, within +/-DIST_TOL. Reads whose insert is
    entirely 5' of the barcode (wrong orientation) fail; reads whose alignment overlaps
    the barcode ('straddle', an alignment artifact) are kept. Expected gap = the per-
    barcode UMI+spacer offset (GRUN_OFFSET) + GGG_LEN."""
    if not cig or cig == '*' or m is None:
        return False
    qas = 0
    for nn, oo in _CIG.findall(cig):
        if oo in 'SH':
            qas += int(nn)
        else:
            break
    aq = sum(int(nn) for nn, oo in _CIG.findall(cig) if oo in 'MI=X')
    qae = qas + aq
    L = len(seq)
    aS, bS = (qas, qae) if m.get('orient') == 'F' else (L - qae, L - qas)
    # the matched code's own length (set D 11-mers, set C 14-mers), not a fixed 7
    e = m['bc_pos'] + int(m.get('bc_len') or len(m.get('code_seq') or m.get('seven') or '') or 7)
    if aS >= e - 3:                              # insert on the barcode's 3' side
        # Prefer the 7-mer the detector actually matched. Using the module constants
        # ox.BC_HUMAN/ox.BC_MOUSE (module defaults) would look up the WRONG key for any
        # other bobcode pair and silently fall back to offset 8, so the expected gap
        # would be wrong for every +11 barcode.
        seven = m.get('seven')
        if not seven:
            seven = ox.BC_HUMAN if m.get('code') == 'h' else ox.BC_MOUSE
        exp = (GRUN_OFFSET[seven] if seven in GRUN_OFFSET else _default_offset(seven)) + GGG_LEN
        return abs((aS - e) - exp) <= DIST_TOL
    if bS <= e + 3:                              # insert entirely 5' of the barcode
        return False                            # wrong side (chimera)
    return True                                 # straddle: alignment overlaps barcode -> keep


def _level(seq, aqlen, m, cat, asp, dist_ok):
    """Deepest CANONICAL funnel level 0..6, BARCODE-FIRST so levels 1..6 are all
    scorable (one barcode, one unit, hu/ms-aligned): has-barcode → mRNA-only →
    barcode-insert distance → G-run ≥3 → insert. Level 0 =
    not scorable. `m` = measure_struct_24plex(seq); `cat` = 'mrna'/'rrna'/'mt'/'other';
    `asp` = 'human'/'mouse'/None; `dist_ok` = _dist_ok(seq, cigar, m)."""
    if m is None or asp is None:
        return 0                                # not scorable (no 7-mer or no hu/ms alignment)
    _pos = getattr(ox, 'READ_LAYOUT', 'ont_anchor') == 'r2_positional'
    # A read with two or more TSO units is a TSO-TSO chimera: not scorable. On Illumina
    # (r2_positional) the aligned SEQ holds no barcode and the chimera screen already
    # happened at prep (a QNAME bobcode exists only if exactly one was called), so the
    # unit test is skipped. There is no separate 'exact 1x barcode' level: prep records
    # no mismatch count, so on Illumina it would duplicate the level above.
    if not _pos and len(ox.find_24plex_units(seq)) != 1:
        return 0
    if cat != 'mrna':
        return 1                                # has one barcode, but not protein-coding mRNA
    if not dist_ok:
        return 2                                # mRNA, but barcode->insert distance wrong
    # G-run >=3. On the Illumina path `dist_ok` is already `grun > 0`, so this only
    # removes runs of 1-2 (~2.4% of reads, ~10x the error rate of G-run >=3).
    if (m.get('grun') or 0) < 3:
        return 3                                # + distance, but G-run <3
    if aqlen <= 30:
        return 4                                # + G-run, but insert <30 bp
    return 5                                    # + insert = full structure
    # NOTE: _rt_reach() / _ts2_partial() are unreferenced from this ladder. They are
    # left in the module (they are the ONT full-molecule probes) but nothing in the
    # Illumina path calls them — do not read their presence as evidence those handles
    # are being checked.


def compute_full_mol_funnel(bam, gtf=None, rdna=None):
    """Return {'n_total', 'steps':[{key,label,n,pct,acc,acc_n,acc_correct}]} for a
    STAR BAM, or None. Per level, `acc` = barcode accuracy = of reads reaching that
    level that carry a species 7-mer AND a primary alignment on a human/mouse
    chromosome, the fraction whose 7-mer species matches the aligned species. gtf/
    rdna (star_taxonomy caches) drive the rRNA/MT filter; loaded if not supplied."""
    if not bam or not os.path.exists(bam):
        return None
    if gtf is None:
        gtf = st.parse_gtf(st.GTF)
    if rdna is None:
        rdna = st.load_rdna(st.RDNA_BED)
    try:
        pv = subprocess.Popen(['samtools', 'view', '-F', '0x900', bam],
                              stdout=subprocess.PIPE, text=True)
    except OSError:
        return None
    nlev = len(STEPS) + 1
    counts = [0] * nlev
    scored = [0] * nlev                          # barcoded + species-aligned, by level
    correct = [0] * nlev                         # of those, 7-mer species == aligned
    # per-aligned-species scored/correct, by level (for the mouse/human split)
    sp_scored = {'human': [0] * nlev, 'mouse': [0] * nlev}
    sp_correct = {'human': [0] * nlev, 'mouse': [0] * nlev}
    # 50-bp aligned-length bin -> [correct, total] over reads that reached the
    # barcode-insert distance step (level >= dist), for the length-accuracy plot.
    cvl = {}
    _dist_lvl = [k for k, _ in STEPS].index('dist') + 1
    _DUMP_TAG = os.path.basename(bam)
    # Reads whose alignment records hit BOTH genomes cannot be given a species; scored by
    # their primary contig they read as barcode swaps (0.15 points on one library). Exclude
    # them from the scorable population and report how many, as star_taxonomy's own
    # correctness metric already does.
    _seen = {'human': set(), 'mouse': set()}
    _pa = subprocess.Popen(['samtools', 'view', '-F', '0x4', bam], stdout=subprocess.PIPE, text=True)
    for _l in _pa.stdout:
        _q, _, _c = _l.split('\t', 3)[:3]
        if _c.startswith('HUMAN_'): _seen['human'].add(_q)
        elif _c.startswith('MOUSE_'): _seen['mouse'].add(_q)
    _pa.wait()
    ambiguous = _seen['human'] & _seen['mouse']; n_amb = 0
    del _seen
    n = 0
    for line in pv.stdout:
        f = line.split('\t')
        seq = f[9]
        if seq == '*':
            continue
        n += 1
        mapped = not (int(f[1]) & 0x4)
        ch = f[2]
        asp = None
        cat = None
        if mapped:
            aq = _aligned_qlen(f[5])
            pos = int(f[3])
            cat = _cat_full(ch, pos - 1, _refend(pos, f[5]), gtf, rdna)
            asp = 'human' if ch.startswith('HUMAN_') else 'mouse' if ch.startswith('MOUSE_') else None
            if f[0] in ambiguous:
                asp = None; n_amb += 1
        else:
            aq = 0
        # ILLUMINA r2_positional: structure comes from the QNAME (prep_reads writes it
        # there because the bobcode/UMI/G-run are TRIMMED OFF before alignment and so are
        # NOT in the BAM SEQ). Falling back to the sequence path would return None for
        # every read and silently zero the whole funnel.
        if getattr(ox, 'READ_LAYOUT', 'ont_anchor') == 'r2_positional':
            m = ox.measure_from_qname(f[0])
            # The barcode->insert distance is enforced at TRIM time (the insert begins
            # exactly after the confirmed G-run), so the equivalent evidence here is that
            # a G-run was actually found at the expected per-bobcode offset.
            do = bool(m and m.get('grun', 0) > 0)
        else:
            m = ox.measure_struct_24plex(seq)
            do = _dist_ok(seq, f[5], m) if (m is not None and asp is not None) else False
        lvl = _level(seq, aq, m, cat, asp, do)
        counts[lvl] += 1
        if lvl >= 1:                             # scorable: barcoded + hu/ms-aligned (levels 1..10)
            ok = int(m['species'] == asp)
            scored[lvl] += 1; correct[lvl] += ok
            sp_scored[asp][lvl] += 1; sp_correct[asp][lvl] += ok
            if _DUMP_DISCORDANT and lvl >= _dist_lvl and not ok:
                # same gate as the ms/hu split below (bi == _dist_lvl), so the dumped
                # set has exactly the reported hu_n / ms_n size.
                _DUMP_FH.write(f'{f[0]}\t{asp}\t{_DUMP_TAG}\n')
            if lvl >= _dist_lvl:                 # reads through the barcode-insert distance step
                _e = cvl.setdefault((aq // 50) * 50, [0, 0]); _e[0] += ok; _e[1] += 1
    pv.wait()
    if not n:
        return None
    steps = []
    for i, (k, lbl) in enumerate(STEPS, start=1):
        cum = sum(counts[i:])                   # reads reaching AT LEAST this level
        sc = sum(scored[i:]); co = sum(correct[i:])
        steps.append({'key': k, 'label': lbl, 'n': cum,
                      'pct': round(100 * cum / n, 1),
                      'acc': round(100 * co / sc, 2) if sc else None,
                      'acc_n': sc, 'acc_correct': co})
    # Mouse/Human % correct split, by ALIGNED species: of reads aligning to that
    # species, the fraction correctly barcoded. Split at the barcode-insert distance
    # level (mRNA + distance-filtered), so it tracks the accuracy after the distance QC.
    bi = [k for k, _ in STEPS].index('dist') + 1    # mouse/human split at the distance-filter level
    out = {'n_total': n, 'steps': steps, 'correct_vs_len': cvl, 'n_ambiguous_excluded': n_amb}
    for asp, pfx in (('mouse', 'ms'), ('human', 'hu')):
        sn = sum(sp_scored[asp][bi:]); sc = sum(sp_correct[asp][bi:])
        out[f'{pfx}_correct'] = round(100 * sc / sn, 2) if sn else None
        out[f'{pfx}_n'] = sn
        out[f'{pfx}_correct_n'] = sc          # count behind the % (library merge adds these)
    return out


def render_full_mol_funnel_figure(all_results, bc_list, bc_short, output_dir):
    """Horizontal cumulative funnel (one colored bar-pair per BC). Returns PNG path."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    data = {}
    for bc in bc_list:
        ff = all_results[bc].get('fullmol_funnel')
        if ff and ff.get('n_total'):
            data[bc] = ff
    if not data:
        return None
    labels = [s[1].replace(' = FULL molecule', '\n= FULL molecule') for s in STEPS]
    cmap = plt.get_cmap('tab10')
    fig, ax = plt.subplots(figsize=(9, 0.55 * len(labels) + 1.4 + 0.25 * len(data)))
    y = np.arange(len(labels))[::-1]
    nbc = len(data)
    h = 0.8 / max(nbc, 1)
    for i, (bc, ff) in enumerate(data.items()):
        off = (nbc - 1) / 2.0 * h - i * h
        pct = [s['pct'] for s in ff['steps']]
        cnt = [s['n'] for s in ff['steps']]
        col = cmap(i % 10)
        ax.barh(y + off, pct, height=h * 0.92, color=col, alpha=0.85,
                label=f"{bc_short.get(bc, bc)} (n={ff['n_total']:,})")
        for yy, c, p in zip(y + off, cnt, pct):
            ax.text(p + 0.5, yy, f'{c:,} ({p:.1f}%)', va='center', fontsize=6.5, color=col)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    xmax = max(s['pct'] for ff in data.values() for s in ff['steps'])
    ax.set_xlim(0, min(100, xmax * 1.35 + 6))
    ax.set_xlabel('% of reads (cumulative — each row is a subset of the one above)', fontsize=8)
    ax.set_title('Read structure — full ordered-molecule funnel\n'
                 'insert>30 · 7mer · polydT · 26N UMI · TruSeq1 · TruSeq2',
                 fontsize=10, fontweight='bold')
    ax.legend(frameon=False, fontsize=7.5, loc='lower right')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    out = os.path.join(output_dir, 'read_structure_full_molecule_funnel.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return out

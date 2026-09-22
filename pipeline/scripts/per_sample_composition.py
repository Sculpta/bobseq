#!/usr/bin/env python3
"""Per-sample read composition and rRNA-excluded species accuracy.

WHY THIS EXISTS
    Both library-level accuracy paths collapse a multiplex run into two
    species buckets:

      star_crosstab   keyed by adapter.CODE_SEQ, which analyze_speciesmix fills with ONE
                      human and ONE mouse code taken by dict order. Since the 24plex set
                      has minimum pairwise Hamming distance 4, no other code is ever
                      within the 1-mismatch call radius, so the other 22 samples score as
                      "no barcode".
      filter_funnel   DOES see all 24 barcodes (DETECT_BARCODES) but maps each to 'h'/'m'
                      and aggregates, so there is no per-code breakdown.

    Neither is wrong for a two-barcode species-mixing run, which is what they were built
    for. Neither can answer "how much does SAMPLE 17 swap", which is what a 24-plex
    experiment needs.

    This computes it per sample, as a separate step, reusing star_taxonomy's own
    classify/parse primitives so the category definitions cannot drift from the report.

THE rRNA CAVEAT THIS EXISTS TO HANDLE
    rDNA is represented in GRCh38 but largely absent from GRCm39, so conserved rRNA
    force-maps to human. A mouse sample therefore looks far less mouse-specific over ALL
    reads than over mRNA. Both views are emitted; the mRNA-only view is the one that
    reflects barcode swapping rather than reference asymmetry.

CLI: per_sample_composition.py --bams <bam>... --run-json <cfg> --out <tsv>
"""
import argparse, collections, json, os, subprocess, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'run'))
import star_taxonomy_illumina as st                                                # noqa: E402

CATS = ['mRNA (protein-coding)', 'Ribosomal-protein mRNA',
        'Other ncRNA (lncRNA/NMD/ret-intron/pseudo)', 'rRNA', 'Mitochondrial',
        'pre-mRNA / intronic', 'Intergenic / genomic', 'Unmapped']
# 'Unmapped' = reads assigned to the sample that produced no species-called alignment.
# It belongs in the composition: with it every class is a share of the sample's OWN reads
# and the shares sum to 100. Without it the denominator is the mapped subset, so a sample
# where most reads fail to align looks compositionally normal. Measured on a 24-plex
# library: the unmapped share runs 21% to 63% per sample, so leaving it out would be
# materially misleading.
# The categories star_taxonomy excludes from every barcode-correctness metric, for the
# reasons documented there: rRNA (reference asymmetry), MT (excluded by convention),
# intergenic (poly-dT primed, arbitrary species), unmapped (no species).
EXCLUDE = set(st.CORRECTNESS_EXCLUDE)


def overlap(ivl, starts, c, s, e, collect=False):
    return st._overlap(ivl, starts, c, s, e, collect=collect)


def library_reads(fastqs, library):
    """Read names belonging to one ONT barcode, taken from the FASTQ STAR was given.

    Dorado writes the barcode into every read header (SM:Z:barcodeNN), and that tag
    survives demux and adapter trimming, so the library a read came from is carried by
    the read itself. That makes the per-library split exact and alignment-free: no
    re-running STAR, no read-name bookkeeping that could drift out of sync.
    """
    tag = 'SM:Z:' + library
    keep = set()
    for fq in fastqs:
        if not os.path.exists(fq):
            continue
        with open(fq) as fh:
            for i, line in enumerate(fh):
                if i % 4:
                    continue
                if tag in line:
                    keep.add(line[1:].split(chr(9))[0].split()[0])
    return keep


def classify_bam(path, genes, gst, exons, est, ridna, rdst, extra=None, keep=None):
    """category -> {human, mouse} for one sample BAM, using star_taxonomy's own rules.

    `extra` also accumulates the metrics the library report gives per ONT barcode but
    which are genuinely per-SAMPLE properties, using star_taxonomy's own definitions:
      insert   _qaln(CIGAR): query bases aligned (M/I/=/X), the "aligned insert length"
      error    (sub + I + D) / (M + I + D), with substitutions from STAR's nM tag
      dup      reads sharing a primary-alignment footprint (chrom, strand, start, end).
               Position-based, so it tolerates ONT basecall error. Reported over ALL
               aligned reads AND over protein-coding mRNA only -- star_taxonomy notes the
               total is inflated because GC-rich rDNA basecalls near-identically across
               INDEPENDENT molecules, so the mRNA rate is the informative one.
    """
    out = collections.defaultdict(lambda: {'human': 0, 'mouse': 0})
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return out
    p = subprocess.Popen(['samtools', 'view', '-F', '0x904', path],
                         stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        f = line.split('\t', 12)
        if keep is not None and f[0] not in keep:
            continue                       # a read from the other ONT barcode
        c = f[2]
        sp = 'human' if c.startswith('HUMAN_') else 'mouse' if c.startswith('MOUSE_') else None
        if sp is None:
            continue
        s = int(f[3]) - 1
        e = st._refend(int(f[3]), f[5])
        # same order of tests as star_taxonomy.classify_one
        if overlap(ridna, rdst, c, s, e):
            k = 'rRNA'
        elif c.endswith('_MT'):
            k = 'Mitochondrial'
        else:
            gs = overlap(genes, gst, c, s, e, collect=True)
            if not gs:
                k = 'Intergenic / genomic'
            elif not overlap(exons, est, c, s, e):
                k = 'pre-mRNA / intronic'
            else:
                bts = [g[2] for g in gs]
                if any(b in st._RRNA_BT for b in bts):
                    k = 'rRNA'
                else:
                    pc = [g for g in gs if g[2].startswith('protein_coding')]
                    if pc:
                        k = ('Ribosomal-protein mRNA'
                             if any(st._RIBO.match(g[3] or '') for g in pc)
                             else 'mRNA (protein-coding)')
                    else:
                        k = 'Other ncRNA (lncRNA/NMD/ret-intron/pseudo)'
        out[k][sp] += 1
        if extra is not None:
            cig = f[5]
            extra['insert'].append(st._qaln(cig))
            M, I, D = st._cigar_mid(cig)
            nm = st._NM_RE.search(line)
            extra['err_M'] += M; extra['err_I'] += I; extra['err_D'] += D
            extra['err_sub'] += min(int(nm.group(1)), M) if nm else 0
            strand = '-' if (int(f[1]) & 0x10) else '+'
            fp = (c, strand, int(f[3]), e)
            extra['fp_all'][fp] += 1
            if k == 'mRNA (protein-coding)':
                extra['fp_mrna'][fp] += 1
    p.wait()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bams', nargs='+', required=True)
    ap.add_argument('--run-json', required=True)
    ap.add_argument('--per-sample', default=None,
                    help='per_sample.tsv, for reads_assigned. Without it the composition '
                         'denominator is species-called reads only, which HIDES the '
                         'unmapped fraction -- 21 to 63 percent on 260826.')
    ap.add_argument('--library', default=None,
                    help='restrict to one ONT barcode, e.g. barcode12. The two barcodes '
                         'are two PCR primer sets on the SAME pool, so this is the axis '
                         'that comparison runs along.')
    ap.add_argument('--library-fastq', nargs='*', default=None,
                    help='the trimmed FASTQs STAR was given; their read headers carry '
                         'the ONT barcode. Required with --library.')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    if args.library and not args.library_fastq:
        ap.error('--library needs --library-fastq')

    # A MISSING --per-sample file must not degrade silently: `assigned` would stay
    # empty, the denominator would fall back to mapped-only, Unmapped would come out 0,
    # and the report would draw a funnel whose "mapped" row is actually the assigned
    # count. Fail instead.
    if args.per_sample and not os.path.exists(args.per_sample):
        ap.error(f'--per-sample file not found: {args.per_sample}. Without it the '
                 f'composition denominator silently becomes mapped-only, which hides '
                 f'the unmapped fraction and inflates every other class.')
    assigned = {}
    if args.per_sample:
        rows = [l.rstrip(chr(10)).split(chr(9)) for l in open(args.per_sample)
                if not l.startswith('#')]
        hh = rows[0]
        for r in rows[1:]:
            d = dict(zip(hh, r))
            assigned[d['sample']] = int(d['reads_assigned'])

    cfg = json.load(open(args.run_json))
    illumina = str(cfg.get('platform', '')).lower().startswith('illumina')
    # unlabeled run: the sample IS the code (demux / per_sample_table name it the same way)
    labels = cfg.get('bobcode_labels') or {c: c for c in (cfg.get('tso_species_map') or {})}
    species = cfg.get('tso_species_map') or {}
    by_sample = {}
    for seq, name in labels.items():
        by_sample[name] = species.get(seq, '')

    genes, gst, exons, est = st.parse_gtf(os.environ['BOBSEQ_GTF'])
    ridna, rdst = st.load_rdna(os.environ['BOBSEQ_RDNA_BED'])

    # the same sample appears under both ONT barcodes; pool them
    grouped = collections.defaultdict(list)
    for b in args.bams:
        stem = os.path.basename(b).replace('_star_combined.bam', '')
        grouped[stem].append(b)

    def fname(s):
        k = ''.join(c if (c.isalnum() or c in '-_.') else '_' for c in s)
        while '__' in k:
            k = k.replace('__', '_')
        return k.strip('_')

    # In library mode the denominator is this library's own reads for this sample -- the
    # exact set STAR was handed -- so unmapped stays visible per library instead of being
    # diluted by the other primer set.
    keep_by = {}
    if args.library:
        fq_by = collections.defaultdict(list)
        for f in (args.library_fastq or []):
            stem = os.path.basename(f)
            for pre in ('noadapter_combined_', 'combined_'):
                if stem.startswith(pre):
                    stem = stem[len(pre):]
                    break
            fq_by[stem.replace('.fastq', '')].append(f)
        for stem, fqs in fq_by.items():
            keep_by[stem] = library_reads(fqs, args.library)

    rows = []
    for sample, exp in by_sample.items():
        tal = collections.defaultdict(lambda: {'human': 0, 'mouse': 0})
        ex = {'insert': [], 'err_M': 0, 'err_I': 0, 'err_D': 0, 'err_sub': 0,
              'fp_all': collections.Counter(), 'fp_mrna': collections.Counter()}
        keep = keep_by.get(fname(sample)) if args.library else None
        if args.library and keep is None:
            keep = set()                    # sample absent from this library: zero reads
        for b in grouped.get(fname(sample), []):
            for k, v in classify_bam(b, genes, gst, exons, est, ridna, rdst,
                                     ex, keep).items():
                tal[k]['human'] += v['human']; tal[k]['mouse'] += v['mouse']
        mapped = sum(v['human'] + v['mouse'] for v in tal.values())
        # denominator = the sample's own assigned reads when known, so the unmapped
        # share is visible; fall back to mapped-only when it is not.
        tot = len(keep) if args.library else assigned.get(sample, mapped)
        tal['Unmapped'] = {'human': max(0, tot - mapped), 'mouse': 0}
        # accuracy is over SPECIES-CALLED reads only: an unmapped read has no species,
        # so counting it would not be a swap rate.
        allh = sum(v['human'] for k, v in tal.items() if k != 'Unmapped')
        allm = sum(v['mouse'] for k, v in tal.items() if k != 'Unmapped')
        # 'Unmapped' must be excluded here too. EXCLUDE carries star_taxonomy's
        # vocabulary ('Unmapped / contaminant'); the bucket added above is named
        # 'Unmapped', so `k not in EXCLUDE` alone would let every unmapped read through
        # as a HUMAN read -- which barely moves human samples but collapses mouse ones
        # (a mouse sample at ~28% instead of ~96%).
        mh = sum(v['human'] for k, v in tal.items()
                 if k not in EXCLUDE and k != 'Unmapped')
        mm = sum(v['mouse'] for k, v in tal.items()
                 if k not in EXCLUDE and k != 'Unmapped')
        def pct(h, m):
            n = h + m
            if not n:
                return 0.0
            return round(100.0 * (h if exp == 'human' else m) / n, 2)
        ins = sorted(ex['insert'])
        aln = ex['err_M'] + ex['err_I'] + ex['err_D']
        def duppct(counter):
            n = sum(counter.values())
            d = sum(v for v in counter.values() if v > 1)
            return round(100.0 * d / n, 2) if n else 0.0
        r = {'sample': sample, 'expected_species': exp, 'reads_classified': mapped,
             'reads_assigned': tot,
             'pct_correct_all_reads': pct(allh, allm),
             'pct_correct_mrna_only': pct(mh, mm),
             'n_scorable_mrna': mh + mm,
             'insert_median': ins[len(ins) // 2] if ins else 0,
             'insert_mean': round(sum(ins) / len(ins), 1) if ins else 0.0,
             'aln_error_pct': round((ex['err_sub'] + ex['err_I'] + ex['err_D'])
                                    / aln * 100, 2) if aln else 0.0,
             # Illumina BAMs are post-dedup (definition D), so a positional duplicate
             # rate over them is not a duplicate rate: the per-sample rate lives in the
             # dedup sidecar (dup_pct_defD_aligned). NA keeps the column, not the number.
             'dup_pct_all_INFLATED': 'NA' if illumina else duppct(ex['fp_all']),
             'dup_pct_mrna': 'NA' if illumina else duppct(ex['fp_mrna'])}
        for c in CATS:
            r[c] = round(100.0 * (tal[c]['human'] + tal[c]['mouse']) / tot, 2) if tot else 0.0
        rows.append(r)

    cols = (['sample', 'expected_species', 'reads_assigned', 'reads_classified',
             'pct_correct_all_reads',
             'pct_correct_mrna_only', 'n_scorable_mrna', 'insert_median', 'insert_mean',
             'aln_error_pct', 'dup_pct_all_INFLATED', 'dup_pct_mrna'] + CATS)
    rows.sort(key=lambda r: r['sample'])
    with open(args.out, 'w') as fh:
        if illumina:
            fh.write('# Illumina: the per-sample BAMs are already deduplicated (definition D), so\n'
                     '# dup_pct_all_INFLATED and dup_pct_mrna are NA here; the per-sample duplicate\n'
                     '# rate is dup_pct_defD_aligned in samples/<sample>/<sample>_dup_prefilter.json.\n'
                     '# insert_* is the aligned length per read (capped by the read length);\n'
                     '# aln_error_pct = substitutions + indels over aligned bases (STAR nM).\n')
        else:
            fh.write('# dup_pct_all_INFLATED is an UPPER BOUND, not a duplicate rate. GC-rich rDNA\n'
                 '# basecalls near-identically across INDEPENDENT molecules, so it counts them\n'
                 '# as duplicates. Measured on 260826: median 31.18% vs 2.84% for dup_pct_mrna,\n'
                 '# an 11x inflation. QUOTE dup_pct_mrna. Never quote dup_pct_all_INFLATED.\n'
                 '# pct_correct_all_reads is depressed for MOUSE samples by reference\n'
                 '# asymmetry: rDNA is in GRCh38 but largely absent from GRCm39, so rRNA\n'
                 '# force-maps human. pct_correct_mrna_only excludes rRNA/MT/intergenic\n'
                 '# and is the view that reflects barcode swapping.\n')
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r[c]) for c in cols) + '\n')
    sys.stderr.write(f'per_sample_composition: {len(rows)} samples -> {args.out}\n')


if __name__ == '__main__':
    main()

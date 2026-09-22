#!/usr/bin/env python3
"""One row per sample: reads, assignable %, alignment, condition and replicate.

The unit here is the BOBCODE, pooled across ONT barcodes. The two ONT barcodes in a
pooled run are PCR primer sets applied to one library, so they are a technical dimension,
not a sample dimension; splitting every sample into two half-depth rows would hide the
sample and show the primer set.

WHAT IS DELIBERATELY ABSENT
    Barcode accuracy and swapping. Both are statements about cross-sample confusion, so
    they have no meaning inside a single sample. They stay in the library-level report.

ASSIGNABLE %
    Reported per sample and in total, because a per-sample number rests on the subset of
    reads carrying a readable bobcode (~73% of primary-mapped reads on a 24-plex library). If that
    fraction differs BY SAMPLE it is a bias rather than a loss: a sample whose bobcode is
    harder to read looks under-represented for a technical reason. The columns needed to
    test that are here.
"""
import argparse, bisect, collections, json, os, re, subprocess, sys

REP = re.compile(r'[-\s](\d+)$')


def condition_of(sample):
    """'Ris dose 25mM-2' -> ('Ris dose 25mM', 2). Sample sheets use either separator --
    'Hek control 1' a space, 'Ris dose 25mM-1' a dash -- so both are handled.
    A name with no trailing index is its own condition at replicate 1."""
    m = REP.search(sample)
    if not m:
        return sample, 1
    return sample[:m.start()].strip(), int(m.group(1))


def library_reads(fastqs, library):
    """Read names belonging to one ONT barcode (see per_sample_composition for why the
    read header is the authoritative source)."""
    tag = 'SM:Z:' + library
    keep = set()
    for fq in fastqs:
        if not os.path.exists(fq):
            continue
        with open(fq) as fh:
            for i, line in enumerate(fh):
                if i % 4 == 0 and tag in line:
                    keep.add(line[1:].split(chr(9))[0].split()[0])
    return keep


def load_rdna(bed):
    """{contig: [(start, end), ...]} from the rDNA loci BED. Empty dict if no BED given."""
    loci = collections.defaultdict(list)
    if not bed or not os.path.exists(bed):
        return {}
    with open(bed) as fh:
        for line in fh:
            if line.startswith(('#', 'track', 'browser')):
                continue
            f = line.split('\t')
            if len(f) < 3:
                continue
            loci[f[0]].append((int(f[1]), int(f[2])))
    return {c: sorted(v) for c, v in loci.items()}


def _in_rdna(loci, chrom, pos):
    v = loci.get(chrom)
    if not v:
        return False
    i = bisect.bisect_right(v, (pos, float('inf'))) - 1
    return i >= 0 and v[i][0] <= pos <= v[i][1]


def bam_counts(path, keep=None, rdna=None):
    """(mapped, human, mouse, human_rdna, mouse_rdna) primary records. An empty
    placeholder BAM is a real result: a sample with no reads, not a failure.

    Species comes from the contig prefix on the combined reference, which is the same
    signal star_taxonomy uses.

    WHY THE rDNA COUNTS. rDNA is represented in GRCh38 but largely absent from GRCm39, so
    conserved rRNA from a MOUSE sample force-maps to human and reads as a barcode swap.
    Measured on a 24-plex library: 81-87% of the human-mapped reads in the three RAW
    (mouse) samples fall in rDNA loci, and excluding them moves apparent specificity from
    82.5/80.4/91.2% to 97.1/95.2/98.5%. A 14.6 point artifact, entirely one-directional
    and entirely against the minority species. The size scales with rRNA content, so it is
    worst exactly when purification underperforms."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return 0, 0, 0, 0, 0
    try:
        p = subprocess.Popen(['samtools', 'view', '-F', '0x904', path],
                             stdout=subprocess.PIPE, text=True)
        loci = rdna or {}
        n = hu = mo = hur = mor = 0
        for line in p.stdout:
            f = line.split('\t', 4)
            if keep is not None and f[0] not in keep:
                continue                   # a read from the other ONT barcode
            n += 1
            if f[2].startswith('HUMAN_'):
                hu += 1
                if loci and _in_rdna(loci, f[2], int(f[3])):
                    hur += 1
            elif f[2].startswith('MOUSE_'):
                mo += 1
                if loci and _in_rdna(loci, f[2], int(f[3])):
                    mor += 1
        p.wait()
        return n, hu, mo, hur, mor
    except Exception as e:
        sys.stderr.write(f'per_sample_table: cannot read {path}: {e}\n')
        return 0, 0, 0, 0, 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run-json', required=True)
    ap.add_argument('--demux-stats', nargs='+', required=True)
    ap.add_argument('--bams', nargs='*', default=[])
    ap.add_argument('--library', default=None,
                    help='restrict to one ONT barcode, e.g. barcode12.')
    ap.add_argument('--library-fastq', nargs='*', default=None,
                    help='trimmed FASTQs whose headers carry the ONT barcode.')
    ap.add_argument('--rdna-bed', default=None,
                    help='rDNA loci BED. Without it the rDNA-excluded column cannot be '
                         'computed and the table falls back to the inflated view only, '
                         'which is exactly the confusion this flag exists to prevent.')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    rdna_loci = load_rdna(args.rdna_bed)
    if not rdna_loci:
        sys.stderr.write(
            'per_sample_table: WARNING no rDNA BED, so pct_expected_species_no_rdna will\n'
            '  equal the inflated column. On a two-species run this UNDERSTATES the minority\n'
            '  species by up to 15 points, because conserved rRNA force-maps to the genome\n'
            '  that actually carries rDNA (GRCh38 does, GRCm39 largely does not).\n')
    else:
        sys.stderr.write(f'per_sample_table: rDNA loci loaded for '
                         f'{len(rdna_loci)} contigs\n')

    cfg = json.load(open(args.run_json))
    species = cfg.get('tso_species_map') or {}
    # unlabeled run: the sample IS the code (demux names the lanes the same way)
    labels = cfg.get('bobcode_labels') or {c: c for c in species}

    # Pool the per-ONT-barcode demux stats: the sample is the bobcode, not the primer set.
    reads = collections.Counter()
    tot_all = tot_assigned = 0
    per_ont = {}
    for f in args.demux_stats:
        s = json.load(open(f))
        per_ont[s.get('ont_barcode') or os.path.basename(f)] = s
        tot_all += s['reads_total']
        tot_assigned += s['reads_assigned']
        for code, rec in (s.get('per_bobcode') or {}).items():
            reads[code] += rec['reads']

    bam_by_sample = {}
    for b in args.bams:
        # One BAM per sample: the primer sets are merged before alignment, so the
        # filename is the sample name and needs no de-prefixing.
        stem = os.path.basename(b).replace('_star_combined.bam', '')
        bam_by_sample.setdefault(stem, []).append(b)

    def fname(sample):
        keep = ''.join(c if (c.isalnum() or c in '-_.') else '_' for c in sample)
        while '__' in keep:
            keep = keep.replace('__', '_')
        return keep.strip('_')

    # Per-library mode: 'mapped' must be counted over THIS library's reads only,
    # otherwise a per-library denominator gets the pooled numerator.
    lib_keep = {}
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
            lib_keep[stem] = library_reads(fqs, args.library)

    rows = []
    for code, sample in labels.items():
        cond, rep = condition_of(sample)
        n = reads.get(code, 0)
        mapped = hu = mo = hur = mor = 0
        for x in bam_by_sample.get(fname(sample), []):
            # NOT a, b, c -- 'a' is the argparse namespace in this scope and shadowing it
            # would make args.out an int several lines later.
            n_aln, n_hu, n_mo, n_hur, n_mor = bam_counts(
                x, lib_keep.get(fname(sample), set()) if args.library else None,
                rdna=rdna_loci)
            mapped += n_aln; hu += n_hu; mo += n_mo; hur += n_hur; mor += n_mor
        exp = species.get(code, '')
        called = hu + mo
        # % of species-called reads that landed on the species this sample IS.
        # Ground truth is the experiment, not the pipeline: HEK is human, RAW is mouse.
        correct = (hu if exp == 'human' else mo) if exp else 0
        pct_correct = round(100.0 * correct / called, 2) if called else 0.0
        # The SAME calculation with rDNA removed. This is the one to quote: conserved rRNA
        # from the minority species force-maps to the other genome and reads as a swap.
        hu2, mo2 = hu - hur, mo - mor
        called2 = hu2 + mo2
        correct2 = (hu2 if exp == 'human' else mo2) if exp else 0
        pct_correct2 = round(100.0 * correct2 / called2, 2) if called2 else 0.0
        rows.append({
            'sample': sample, 'condition': cond, 'replicate': rep,
            'bobcode': code, 'species': species.get(code, ''),
            'reads_assigned': n,
            'pct_of_assigned': round(100.0 * n / tot_assigned, 4) if tot_assigned else 0.0,
            'mapped': mapped,
            'pct_mapped': round(100.0 * mapped / n, 2) if n else 0.0,
            'human': hu, 'mouse': mo,
            'human_rdna': hur, 'mouse_rdna': mor,
            'pct_expected_species_no_rdna': pct_correct2,
            'pct_expected_species_ALLREADS_rRNA_INFLATED': pct_correct,
        })
    rows.sort(key=lambda r: (r['condition'], r['replicate']))

    cols = ['sample', 'condition', 'replicate', 'bobcode', 'species',
            'reads_assigned', 'pct_of_assigned', 'mapped', 'pct_mapped',
            'human', 'mouse', 'human_rdna', 'mouse_rdna',
            # The corrected column comes FIRST so the eye lands on it. The uncorrected one
            # keeps its warning in the name, the same convention as dup_pct_all_INFLATED.
            'pct_expected_species_no_rdna',
            'pct_expected_species_ALLREADS_rRNA_INFLATED']
    with open(args.out, 'w') as fh:
        fh.write(f'# reads_total={tot_all}  reads_assigned={tot_assigned}  '
                 f'assignable_pct={100.0*tot_assigned/tot_all:.2f}\n' if tot_all else '#\n')
        fh.write('# unassigned reads still count in full toward the library-level report\n')
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r[c]) for c in cols) + '\n')

    # Within-condition spread: the point of a replicated run.
    by_cond = collections.defaultdict(list)
    for r in rows:
        by_cond[r['condition']].append(r['reads_assigned'])
    sys.stderr.write(f'per_sample_table: {len(rows)} samples, {len(by_cond)} conditions\n')
    for cond, v in sorted(by_cond.items()):
        mean = sum(v) / len(v) if v else 0
        cv = (sum((x - mean) ** 2 for x in v) / len(v)) ** 0.5 / mean * 100 if mean else 0
        sys.stderr.write(f'  {cond:<22} n={len(v)}  reads={v}  CV={cv:.1f}%\n')


if __name__ == '__main__':
    main()

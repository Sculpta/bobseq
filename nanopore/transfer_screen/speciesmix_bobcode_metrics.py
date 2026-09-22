#!/usr/bin/env python3
"""Species-mixing bobcode barcode accuracy (splint, polydT, 5'TSO, 3'TSO): a few metrics per ONT barcode.

Implements bobcode_general_rules.md. Everything that differs between datasets
(chemistry, codes, species, construct anchor, input FASTQ, barcodes to score) is
read from bobcode_run_rules.tsv; nothing run-specific is hardcoded here.

Per ONT barcode:  raw FASTQ -> trim ONT adapter/native barcode -> STARlong
(dense combined human+mouse index) -> classify every read -> score barcodes.

Output (one row per barcode, <out-dir>/bobcode_metrics.tsv):
  read categories (all mapped reads): mRNA, rRNA, MT, other
  human mRNA reads, mouse mRNA reads (protein-coding, species-unique)
  after filtering: barcode accuracy on human mRNA, mouse mRNA, all mRNA (+ n, 95% CI)
Per-barcode JSON (<out-dir>/<run>/<barcode>/metrics.json) holds the same values
plus every filter-step count, the code x species table and provenance.

Requirements: Python >= 3.8 (standard library only), samtools, STARlong 2.7.11b.

Example:
  python3 speciesmix_bobcode_metrics.py \\
      --fastq-dir /path/to/repository_fastqs --out-dir results \\
      --starlong /path/to/STARlong --star-index /path/to/STAR_index \\
      --gtf /path/to/combined_genome.gtf --rdna-bed /path/to/rdna_loci.bed
  --fastq-dir reads the repository FASTQ folder (one file per dataset); --data-dir
  instead reads the original run folders the run rules point to (relative paths).
  Add --run <run_folder> and/or --barcode BC06 to process a subset.
  --bam <file> scores an existing BAM instead of aligning (one run + one barcode).
"""
import argparse
import bisect
import csv
import datetime
import glob
import gzip
import hashlib
import json
import math
import os
import pickle
import platform
import re
import shutil
import subprocess
import sys
from collections import defaultdict

# ---------------------------------------------------------------------------
# Fixed by bobcode_general_rules.md (same for every dataset)
# ---------------------------------------------------------------------------
POST_BC_FLANK = 'CAGCACCT'        # §1 trim: 5' end, remove through this
POST_BC_FLANK_RC = 'AGGTGCTG'     # §1 trim: 3' end, remove from this
TRIM_SEARCH_BP = 80               # §1 window searched at each end

STAR_PARAMS = [                   # §2
    '--outSAMtype', 'BAM', 'Unsorted',
    '--outSAMattributes', 'NH', 'HI', 'AS', 'nM',
    '--outSAMunmapped', 'Within',
    '--outSAMmultNmax', '18446744073709551615',
    '--outFilterMultimapNmax', '10000',
    '--outFilterScoreMinOverLread', '0.3',
    '--outFilterMatchNminOverLread', '0.3',
    '--outFilterMismatchNmax', '1000',
    '--outFilterMismatchNoverLmax', '0.3',
    '--alignIntronMax', '1000000',
    '--alignSJDBoverhangMin', '1',
    '--alignEndsType', 'Local',
    '--seedPerReadNmax', '100000',
]

RRNA_BIOTYPES = {'rRNA', 'Mt_rRNA', 'rRNA_pseudogene'}         # §4
RP_GENE = re.compile(r'^M?RP[LS](?!6K)', re.I)                  # §4 ribosomal-protein genes (not RPS6K kinases)
LOWCX_ENTROPY_MIN = 0.65                                        # §5
UNIT_MAX_MM = 1                                                 # §6 tolerance used only to find units
UNIT_MIN_GAP = 27                                               # §6 hits closer than this = one unit
CODE_RE = re.compile(r'[ACGT]{4,12}')                            # §6 code lengths in use: 4, 6, 7, 8, 10, 12
DIRECTION_MIN_READS = 20                                        # §8 codes with fewer scored reads get no note
MIN_SCORED_READS = 100                                          # §7 datasets with fewer scored mRNA reads are excluded

CIG = re.compile(r'(\d+)([MIDNSHP=X])')
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
CATS = ('mRNA', 'rRNA', 'MT', 'other')

TSV_COLS = [
    'chemistry', 'run', 'barcode', 'plot_no', 'sample_note', 'construct', 'human_code', 'mouse_code',
    'reads_total', 'reads_mapped',
    'mRNA_reads', 'rRNA_reads', 'MT_reads', 'other_reads',
    'human_mRNA_reads', 'mouse_mRNA_reads',
    'human_mRNA_accuracy_pct', 'human_mRNA_scored',
    'mouse_mRNA_accuracy_pct', 'mouse_mRNA_scored',
    'total_mRNA_accuracy_pct', 'total_mRNA_scored',
    'total_ci95_lo_pct', 'total_ci95_hi_pct',
    'status',
]


def rc(s):
    return s.translate(COMP)[::-1]


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def wilson(k, n, z=1.96):
    if not n:
        return None, None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return round(100 * (c - h), 3), round(100 * (c + h), 3)


def pct(k, n):
    return round(100 * k / n, 3) if n else None


# ---------------------------------------------------------------------------
# Rules files
# ---------------------------------------------------------------------------
def load_run_rules(path):
    runs = []
    with open(path, newline='') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if not r.get('run_folder'):
                continue
            bcs = re.findall(r'(BC\d{2})(?:\(#(\d+)\))?', r['barcodes_scored (plot #)'])   # plot # optional (splint: not plotted)
            h, m, a = r['human_code'].strip(), r['mouse_code'].strip(), r['anchor_12nt'].strip()
            if not bcs:
                continue                 # a run listed only to record why it is not scored
            for s in (h, m):
                if not CODE_RE.fullmatch(s):
                    sys.exit(f"run rules: {r['run_folder']}: code '{s}' is not 4-12 nt of ACGT")
            if not re.fullmatch(r'[ACGT]{12}', a):
                sys.exit(f"run rules: {r['run_folder']}: anchor '{a}' is not 12 nt")
            if h == m or h.startswith(m) or m.startswith(h):
                sys.exit(f"run rules: {r['run_folder']}: codes {h}/{m} cannot be told apart")
            runs.append({
                'chemistry': r['chemistry'].strip(), 'date': r['seq_date'].strip().replace('-', '')[2:],
                'run': r['run_folder'].strip(), 'construct': r['construct'].strip(),
                'anchor': a, 'human_code': h, 'mouse_code': m,
                'fastq_pattern': r['raw_fastq_input'].strip(),
                'max_reads': int(r['max_reads']) if (r.get('max_reads') or '').strip().isdigit() else None,   # 'all' = no cap
                'notes': dict(tuple(x.split(':', 1)) for x in (y.strip() for y in (r.get('sample_notes') or '').split('|')) if ':' in x),
                'barcodes': [b for b, _ in bcs], 'plot_no': dict(bcs),
            })
    return runs


# ---------------------------------------------------------------------------
# References: exon index (from GTF) and rDNA loci
# ---------------------------------------------------------------------------
def _merge(ivs):
    ivs.sort()
    st, en = [], []
    for s, e in ivs:
        if st and s <= en[-1]:
            en[-1] = max(en[-1], e)
        else:
            st.append(s); en.append(e)
    return st, en


def load_exon_index(gtf, cache_dir):
    """{'rrna'|'rp'|'pc': {chrom: (starts, ends)}} of merged exon intervals
    (0-based half-open). pc = protein_coding genes that are NOT ribosomal-protein."""
    s = os.stat(gtf)
    key = f'{os.path.abspath(gtf)}|{s.st_size}|{int(s.st_mtime)}|{RP_GENE.pattern}|{sorted(RRNA_BIOTYPES)}|v1'
    os.makedirs(cache_dir, exist_ok=True)
    cache = os.path.join(cache_dir, f'exon_index_{hashlib.sha1(key.encode()).hexdigest()[:12]}.pkl')
    if os.path.exists(cache):
        with open(cache, 'rb') as f:
            return pickle.load(f)
    print(f'  parsing GTF (once; cached to {cache}) ...', file=sys.stderr, flush=True)
    raw = {'rrna': defaultdict(list), 'rp': defaultdict(list), 'pc': defaultdict(list)}
    rp_genes = set()
    bt_re = re.compile(r'gene_biotype "([^"]+)"')
    nm_re = re.compile(r'gene_name "([^"]+)"')
    op = gzip.open if gtf.endswith('.gz') else open
    with op(gtf, 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            p = line.split('\t', 8)
            if len(p) < 9 or p[2] != 'exon':
                continue
            b = bt_re.search(p[8])
            bt = b.group(1) if b else ''
            iv = (int(p[3]) - 1, int(p[4]))
            if bt in RRNA_BIOTYPES:
                raw['rrna'][p[0]].append(iv)
            elif bt == 'protein_coding':
                n = nm_re.search(p[8])
                name = n.group(1) if n else ''
                if RP_GENE.match(name):
                    raw['rp'][p[0]].append(iv); rp_genes.add(name)
                else:
                    raw['pc'][p[0]].append(iv)
    idx = {k: {c: _merge(v) for c, v in d.items()} for k, d in raw.items()}
    idx['_rp_genes'] = sorted(rp_genes)
    with open(cache, 'wb') as f:
        pickle.dump(idx, f)
    return idx


def load_rdna(path):
    d = defaultdict(list)
    with open(path) as f:
        for line in f:
            if not line.strip() or line[0] == '#':
                continue
            p = line.split('\t')
            d[p[0]].append((int(p[1]), int(p[2])))
    return {c: _merge(v) for c, v in d.items()}


def overlaps(index, chrom, s, e):
    """Any merged interval on chrom overlapping [s, e)."""
    t = index.get(chrom)
    if not t:
        return False
    st, en = t
    i = bisect.bisect_right(en, s)
    return i < len(st) and st[i] < e


# ---------------------------------------------------------------------------
# Per-read helpers
# ---------------------------------------------------------------------------
def ref_blocks(pos1, cigar):
    """Aligned reference blocks (0-based half-open), split at N (introns)."""
    cur = bs = pos1 - 1
    out = []
    for n, op in CIG.findall(cigar):
        n = int(n)
        if op in 'M=XD':
            cur += n
        elif op == 'N':
            if cur > bs:
                out.append((bs, cur))
            cur += n; bs = cur
    if cur > bs:
        out.append((bs, cur))
    return out


def aligned_query(seq, cigar):
    """Read bases under M/=/X (soft clips and insertions skipped)."""
    i = 0; out = []
    for n, op in CIG.findall(cigar):
        n = int(n)
        if op in 'M=X':
            out.append(seq[i:i + n]); i += n
        elif op in 'IS':
            i += n
    return ''.join(out)


def dinuc_entropy(s):
    n = len(s)
    if n < 3:
        return 1.0
    di = defaultdict(int)
    for i in range(n - 1):
        di[s[i:i + 2]] += 1
    t = n - 1
    return -sum((c / t) * math.log2(c / t) for c in di.values()) / 4.0


def ham(a, b):
    return sum(x != y for x, y in zip(a, b))


def find_units(seq, anchor, h, m):
    """Barcode units: exact anchor + a candidate within UNIT_MAX_MM of either code
    (each code compared over its own length), both strands; hits < UNIT_MIN_GAP
    apart on a strand collapse. -> [(position, dh, dm)]; a distance is 99 when the
    read ends before that code's length."""
    units = []
    for s in (seq, rc(seq)):
        last = -10 ** 9
        j = s.find(anchor)
        while j >= 0:
            p = j + len(anchor)
            ch, cm = s[p:p + len(h)], s[p:p + len(m)]
            dh = ham(ch, h) if len(ch) == len(h) else 99
            dm = ham(cm, m) if len(cm) == len(m) else 99
            if min(dh, dm) <= UNIT_MAX_MM and p - last >= UNIT_MIN_GAP:
                units.append((p, dh, dm)); last = p
            j = s.find(anchor, j + 1)
    return units


def categorize(prim_chrom, prim_blocks, any_rdna, exidx):
    """§4, first match wins."""
    if any_rdna:
        return 'rRNA'
    if prim_chrom.endswith('_MT'):
        return 'MT'
    if any(overlaps(exidx['rrna'], prim_chrom, s, e) for s, e in prim_blocks):
        return 'rRNA'
    if any(overlaps(exidx['rp'], prim_chrom, s, e) for s, e in prim_blocks):
        return 'other'                       # ribosomal-protein mRNA
    if any(overlaps(exidx['pc'], prim_chrom, s, e) for s, e in prim_blocks):
        return 'mRNA'
    return 'other'                           # lncRNA / other ncRNA / intronic / intergenic


# ---------------------------------------------------------------------------
# Steps
# ---------------------------------------------------------------------------
def trim_read(seq, qual):
    head = seq[:TRIM_SEARCH_BP]
    i5 = head.find(POST_BC_FLANK)
    start = i5 + len(POST_BC_FLANK) if i5 >= 0 else 0
    ts = max(0, len(seq) - TRIM_SEARCH_BP)
    i3 = seq[ts:].rfind(POST_BC_FLANK_RC)
    end = ts + i3 if i3 >= 0 else len(seq)
    if end <= start:
        return '', ''
    return seq[start:end], qual[start:end]


def repository_fastq(fastq_dir, run, bc):
    """§1: the dataset's single file in the repository FASTQ folder,
    <fastq_dir>/<chemistry>_<construct>/<YYMMDD>_<BC>.fastq.gz."""
    return os.path.join(fastq_dir, f"{run['chemistry'].replace(chr(39), '')}_{run['construct']}",
                        f"{run['date']}_{bc}.fastq.gz")


def input_fastqs(data_dir, run, bc, fastq_dir=None):
    """§1: every file matched by the run's pattern(s); ' ; ' separates patterns, each
    relative to the run folder unless absolute. 'barcodeNN' -> 'barcode06' etc.
    With fastq_dir (the repository FASTQ folder) the dataset's one file is used instead."""
    if fastq_dir:
        f = repository_fastq(fastq_dir, run, bc)
        if not os.path.exists(f):
            raise FileNotFoundError(f'no input FASTQ: {f}')
        return [f]
    files = []
    for pat in (x.strip() for x in run['fastq_pattern'].split(';')):
        if not pat:
            continue
        pat = os.path.expandvars(pat.replace('barcodeNN', 'barcode' + bc[2:]))   # e.g. $MINKNOW_DATA/...
        if '$' in pat:
            raise FileNotFoundError(f'unset environment variable in input pattern: {pat} '
                                    '(set MINKNOW_DATA to the MinKNOW output folder, or use --fastq-dir)')
        full = pat if os.path.isabs(pat) else os.path.join(data_dir, run['run'], pat)
        hits = sorted(glob.glob(full))
        if not hits:
            raise FileNotFoundError(f'no input FASTQ: {full}')
        files += hits
    return files


def trim_inputs(files, out_fq, max_reads=None):
    """Concatenate + trim; a read ID seen twice (pooled sources) is kept once.
    max_reads (run rules): stop after that many reads, in file order (§1)."""
    n_in = n_out = 0
    seen = set()
    with open(out_fq + '.part', 'w') as fo:
        for path in files:
            op = gzip.open if path.endswith('.gz') else open
            with op(path, 'rt') as f:
                while True:
                    h = f.readline()
                    if not h:
                        break
                    s = f.readline().rstrip('\n'); f.readline(); q = f.readline().rstrip('\n')
                    rid = h[1:].split(None, 1)[0]
                    if rid in seen:
                        continue
                    if max_reads is not None and n_in >= max_reads:
                        break
                    seen.add(rid)
                    n_in += 1
                    ts, tq = trim_read(s, q)
                    if ts:
                        fo.write(f'{h}{ts}\n+\n{tq}\n'); n_out += 1
    os.replace(out_fq + '.part', out_fq)
    return n_in, n_out


def count_reads(files, max_reads=None):
    """Distinct read IDs across the input files (§1), capped at max_reads."""
    seen = set()
    for path in files:
        op = gzip.open if path.endswith('.gz') else open
        with op(path, 'rt') as f:
            for i, line in enumerate(f):
                if i % 4 == 0:
                    seen.add(line[1:].split(None, 1)[0])
                    if max_reads is not None and len(seen) >= max_reads:
                        return len(seen)
    return len(seen)


def star_version(starlong):
    try:
        return subprocess.run([starlong, '--version'], capture_output=True, text=True).stdout.strip()
    except OSError:
        return None


def run_star(a, fq, prefix):
    busy = subprocess.run(['pgrep', '-x', 'STARlong'], capture_output=True, text=True).stdout.split()
    if busy:
        sys.exit(f'ABORT: a STARlong is already running (pid {" ".join(busy)}). '
                 'Only one may run at a time (each loads the ~53 GB index).')
    cmd = [a.starlong, '--runThreadN', str(a.threads), '--genomeDir', a.star_index,
           '--readFilesIn', fq, '--outFileNamePrefix', prefix] + STAR_PARAMS
    with open(prefix + 'console.txt', 'w') as log:
        rcode = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode
    bam = prefix + 'Aligned.out.bam'
    if rcode != 0 or not os.path.exists(prefix + 'Log.final.out') or not os.path.exists(bam):
        with open(prefix + 'console.txt') as f:
            tail = ''.join(f.readlines()[-15:])
        raise RuntimeError(f'STARlong failed (exit {rcode}):\n{tail}')
    return bam


def score_bam(bam, run, exidx, rdna):
    h, m = run['human_code'], run['mouse_code']
    code_species = {h: 'human', m: 'mouse'}
    c = defaultdict(int)                                   # counts
    acc = {'human': [0, 0], 'mouse': [0, 0]}               # aligned species -> [scored, correct]
    cross = {h: {'human': 0, 'mouse': 0}, m: {'human': 0, 'mouse': 0}}
    seen = set()

    def handle(recs):
        c['reads_in_bam'] += 1
        mapped = [r for r in recs if not int(r[1]) & 4]
        if not mapped:
            c['unmapped'] += 1
            return
        sps = set(); any_rdna = False; prim = None
        for r in mapped:
            ch = r[2]
            sps.add('human' if ch.startswith('HUMAN_') else 'mouse' if ch.startswith('MOUSE_') else ch)
            bl = ref_blocks(int(r[3]), r[5])
            if not any_rdna and bl and overlaps(rdna, ch, bl[0][0], bl[-1][1]):
                any_rdna = True
            if prim is None and not int(r[1]) & 0x900:
                prim = (r, bl)
        if prim is None:                                   # no primary record: use the first
            c['no_primary_record'] += 1
            prim = (mapped[0], ref_blocks(int(mapped[0][3]), mapped[0][5]))
        species = 'human' if sps == {'human'} else 'mouse' if sps == {'mouse'} else 'ambiguous'
        pr, bl = prim
        cat = categorize(pr[2], bl, any_rdna, exidx)
        c['mapped'] += 1; c[cat] += 1; c[f'species_{species}'] += 1
        # ---- filtering steps (§3-§6); only species-unique mRNA reads go on ----
        if cat != 'mRNA':
            return
        if species == 'ambiguous':
            c['mRNA_ambiguous_species'] += 1
            return
        c[f'{species}_mRNA'] += 1
        seq = pr[9]
        if seq == '*':
            seq = next((r[9] for r in recs if r[9] != '*'), '')
        if dinuc_entropy(aligned_query(seq, pr[5])) < LOWCX_ENTROPY_MIN:
            c['f_low_complexity'] += 1
            return
        units = find_units(seq, run['anchor'], h, m)
        if not units:
            c['f_no_barcode_unit'] += 1
            return
        if len(units) > 1:
            c['f_multiple_units'] += 1
            return
        _pos, dh, dm = units[0]
        if dh and dm:
            c['f_code_not_exact'] += 1
            return
        code = h if dh == 0 else m
        cross[code][species] += 1
        acc[species][0] += 1
        acc[species][1] += int(code_species[code] == species)

    p = subprocess.Popen(['samtools', 'view', bam], stdout=subprocess.PIPE, text=True)
    cur, recs = None, []
    for line in p.stdout:
        f = line.split('\t', 11)
        if f[0] != cur:
            if cur is not None:
                handle(recs); seen.add(cur)
            if f[0] in seen:
                raise RuntimeError(f'BAM records of read {f[0]} are not consecutive')
            cur, recs = f[0], []
        recs.append(f)
    if cur is not None:
        handle(recs)
    if p.wait() != 0:
        raise RuntimeError(f'samtools view failed on {bam}')
    return c, acc, cross


def direction_check(run, cross):
    """§8: note (never withhold) when a code's declared species is not the majority
    of its scored reads. Accuracy is always reported."""
    problems = []
    for code, sp in ((run['human_code'], 'human'), (run['mouse_code'], 'mouse')):
        n = sum(cross[code].values())
        if n >= DIRECTION_MIN_READS and cross[code][sp] * 2 <= n:
            problems.append(f'{code} declared {sp} but {cross[code][sp]}/{n} scored reads align {sp}')
    return problems


def process_barcode(a, run, bc, exidx, rdna, prov):
    wd = os.path.join(a.out_dir, os.path.basename(run['run']), bc)
    os.makedirs(wd, exist_ok=True)
    mfile = os.path.join(wd, 'metrics.json')
    if os.path.exists(mfile) and not a.force:
        print(f'  {run["run"]} {bc}: metrics.json present, skipping (--force to redo)', file=sys.stderr)
        with open(mfile) as f:
            return json.load(f)
    files = input_fastqs(a.data_dir, run, bc, a.fastq_dir)
    info = {'input_fastqs': [os.path.relpath(x, a.fastq_dir or a.data_dir) for x in files], 'run_folder': run['run']}
    if a.bam:
        bam = a.bam
        info['mode'] = f'scored an existing BAM: {os.path.abspath(bam)}'
        info['reads_total'] = count_reads(files, run['max_reads'])
        info['reads_after_trim'] = None
    else:
        fq = os.path.join(wd, 'trimmed.fastq')
        tfile = os.path.join(wd, 'trim_counts.json')
        if os.path.exists(fq) and os.path.exists(tfile) and not a.force:
            with open(tfile) as f:
                n_in, n_out = json.load(f)
        else:
            print(f'  {run["run"]} {bc}: trimming {len(files)} FASTQ file(s)', file=sys.stderr, flush=True)
            n_in, n_out = trim_inputs(files, fq, run['max_reads'])
            with open(tfile, 'w') as f:
                json.dump([n_in, n_out], f)
        info['reads_total'], info['reads_after_trim'] = n_in, n_out
        prefix = os.path.join(wd, 'star_')
        bam = prefix + 'Aligned.out.bam'
        if not (os.path.exists(bam) and os.path.exists(prefix + 'Log.final.out')) or a.force:
            print(f'  {run["run"]} {bc}: STARlong on {n_out:,} reads', file=sys.stderr, flush=True)
            bam = run_star(a, fq, prefix)
        info['mode'] = 'trimmed + aligned by this script'
    print(f'  {run["run"]} {bc}: scoring', file=sys.stderr, flush=True)
    c, acc, cross = score_bam(bam, run, exidx, rdna)
    problems = direction_check(run, cross)
    n_tot = acc['human'][0] + acc['mouse'][0]
    k_tot = acc['human'][1] + acc['mouse'][1]
    lo, hi = wilson(k_tot, n_tot)
    row = {
        'chemistry': run['chemistry'], 'run': os.path.basename(run['run']), 'barcode': bc,
        'plot_no': run['plot_no'].get(bc, ''),
        'construct': run['construct'], 'human_code': run['human_code'], 'mouse_code': run['mouse_code'],
        'reads_total': info['reads_total'], 'reads_mapped': c['mapped'],
        'mRNA_reads': c['mRNA'], 'rRNA_reads': c['rRNA'], 'MT_reads': c['MT'], 'other_reads': c['other'],
        'human_mRNA_reads': c['human_mRNA'], 'mouse_mRNA_reads': c['mouse_mRNA'],
        'human_mRNA_accuracy_pct': pct(acc['human'][1], acc['human'][0]),
        'human_mRNA_scored': acc['human'][0],
        'mouse_mRNA_accuracy_pct': pct(acc['mouse'][1], acc['mouse'][0]),
        'mouse_mRNA_scored': acc['mouse'][0],
        'total_mRNA_accuracy_pct': pct(k_tot, n_tot),
        'total_mRNA_scored': n_tot,
        'total_ci95_lo_pct': lo, 'total_ci95_hi_pct': hi,
        'status': ('EXCLUDE: fewer than %d scored mRNA reads; ' % MIN_SCORED_READS if n_tot < MIN_SCORED_READS else '')
                  + ('ok' if not problems else 'note: ' + '; '.join(problems)),
    }
    out = {
        'metrics': row,
        'correct': {'human_mRNA': acc['human'][1], 'mouse_mRNA': acc['mouse'][1], 'total_mRNA': k_tot},
        'counts': dict(c), 'code_by_species': cross, **info,
        'provenance': {**prov, 'bam': os.path.abspath(bam),
                       'processed': datetime.datetime.now().isoformat(timespec='seconds')},
    }
    with open(mfile, 'w') as f:
        json.dump(out, f, indent=1)
    msg = f"total {row['total_mRNA_accuracy_pct']}% (n={n_tot:,})" + ('' if not problems else f"  [{row['status']}]")
    print(f'  {run["run"]} {bc}: {msg}', file=sys.stderr, flush=True)
    return out


def write_tsv(a, runs):
    rows, done = [], set()        # a folder can hold two rows' barcodes (e.g. 260519): read each file once
    owner = {(os.path.basename(r['run']), b): r for r in runs for b in r['barcodes']}
    for run in runs:
        for bc in sorted(set(run['barcodes']) | {os.path.basename(p) for p in
                         glob.glob(os.path.join(a.out_dir, os.path.basename(run['run']), 'BC??'))}):
            f = os.path.join(a.out_dir, os.path.basename(run['run']), bc, 'metrics.json')
            if os.path.exists(f) and f not in done:
                done.add(f)
                with open(f) as fh:
                    row = json.load(fh)['metrics']
                # from the run rules, current version (plot numbers can change after scoring)
                own = owner.get((os.path.basename(run['run']), bc), run)   # the row that lists this barcode
                row['plot_no'] = own['plot_no'].get(bc) or 'not plotted'
                row['sample_note'] = own['notes'].get(bc, '').strip() or 'condition variant'
                rows.append(row)
    path = os.path.join(a.out_dir, 'bobcode_metrics.tsv')
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, TSV_COLS, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow({k: ('' if v is None else v) for k, v in r.items()})
    return path, len(rows)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    env = os.environ.get
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run-rules', default=os.path.join(here, 'bobcode_run_rules.tsv'))
    ap.add_argument('--general-rules', default=os.path.join(here, 'bobcode_general_rules.md'),
                    help='only checksummed into the provenance record')
    ap.add_argument('--data-dir', help='folder the run rules\' run_folder paths are relative to (original MinKNOW layout)')
    ap.add_argument('--fastq-dir', help='the repository FASTQ folder (<chemistry>_<construct>/<YYMMDD>_<BC>.fastq.gz); '
                    'use instead of --data-dir')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--starlong', default=env('SPECIESMIX_STARLONG', 'STARlong'))
    ap.add_argument('--star-index', default=env('SPECIESMIX_STAR_INDEX'))
    ap.add_argument('--gtf', default=env('SPECIESMIX_GTF'))
    ap.add_argument('--rdna-bed', default=env('SPECIESMIX_RDNA_BED', os.path.join(here, 'rdna_loci.bed')))
    ap.add_argument('--threads', type=int, default=1,
                    help='STARlong threads (default 1: multi-threaded STARlong can hang under Rosetta)')
    ap.add_argument('--run', action='append', help='run_folder (as in the run rules, or its last path part) to process; repeatable; default all')
    ap.add_argument('--chemistry', action='append', help="only this chemistry (splint, polydT, 5'TSO, 3'TSO); repeatable")
    ap.add_argument('--exclude-run', action='append', help='run_folder (or its last path part) to leave out; repeatable')
    ap.add_argument('--barcode', action='append', help='e.g. BC06 (repeatable); default the scored barcodes')
    ap.add_argument('--bam', help='score this existing BAM instead of trimming + aligning')
    ap.add_argument('--force', action='store_true', help='redo barcodes that already have results')
    a = ap.parse_args()

    if not (a.data_dir or a.fastq_dir):
        ap.error('give --fastq-dir (repository FASTQs) or --data-dir (original run folders)')
    for need in ('gtf', 'rdna_bed') + (() if a.bam else ('star_index',)):
        if not getattr(a, need):
            ap.error(f'--{need.replace("_", "-")} is required')
    if not shutil.which('samtools'):
        sys.exit('samtools not found on PATH')

    runs = load_run_rules(a.run_rules)
    if a.run:
        want = set(a.run)
        pick = [r for r in runs if r['run'] in want or os.path.basename(r['run']) in want]
        found = {r['run'] for r in pick} | {os.path.basename(r['run']) for r in pick}
        if want - found:
            sys.exit(f'not in run rules: {sorted(want - found)}')
        runs = pick
    if a.chemistry:
        runs = [r for r in runs if r['chemistry'] in set(a.chemistry)]
    if a.exclude_run:
        drop = set(a.exclude_run)
        runs = [r for r in runs if r['run'] not in drop and os.path.basename(r['run']) not in drop]
    if a.barcode:
        # A requested barcode goes to the run row(s) that score it; only a barcode no
        # selected row scores (e.g. an excluded negative) is run under every selected row.
        jobs = []
        for bc in a.barcode:
            owners = [r for r in runs if bc in r['barcodes']] or runs
            jobs += [(r, bc) for r in owners]
    else:
        jobs = [(r, bc) for r in runs for bc in r['barcodes']]
    if a.bam and len(jobs) != 1:
        sys.exit('--bam needs exactly one run and one barcode (--run ... --barcode ...)')

    prov = {
        'script': os.path.basename(__file__), 'script_sha256': sha256(__file__),
        'run_rules_sha256': sha256(a.run_rules),
        'general_rules_sha256': sha256(a.general_rules) if os.path.exists(a.general_rules) else None,
        'gtf': os.path.abspath(a.gtf), 'gtf_bytes': os.path.getsize(a.gtf),
        'rdna_bed_sha256': sha256(a.rdna_bed),
        'star_index': os.path.abspath(a.star_index) if a.star_index else None,
        'starlong_version': None if a.bam else star_version(a.starlong),
        'star_threads': a.threads,
        'samtools_version': subprocess.run(['samtools', '--version'], capture_output=True,
                                           text=True).stdout.split('\n')[0],
        'python': platform.python_version(),
    }
    exidx = load_exon_index(a.gtf, os.path.join(a.out_dir, '_cache'))
    prov['ribosomal_protein_genes_excluded'] = len(exidx['_rp_genes'])
    rdna = load_rdna(a.rdna_bed)

    failed = []
    for run, bc in jobs:
        try:
            out = process_barcode(a, run, bc, exidx, rdna, prov)
        except (FileNotFoundError, RuntimeError) as e:
            failed.append(f"{run['run']} {bc}: {e}")
            print(f"  {run['run']} {bc}: ERROR {e}", file=sys.stderr)
    path, n = write_tsv(a, load_run_rules(a.run_rules))
    print(f'wrote {path} ({n} barcodes)', file=sys.stderr)
    if failed:
        print('PROBLEMS:\n  ' + '\n  '.join(failed), file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()

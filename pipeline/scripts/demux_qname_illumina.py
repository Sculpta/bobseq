#!/usr/bin/env python3
"""Split one Illumina library into per-sample lanes by the bobcode prep wrote into the QNAME.

The Illumina counterpart of the ONT demultiplexer. ONT has to demultiplex BEFORE alignment
because the bobcode is found by an anchor search on the read; on Illumina the bobcode is
positional and prep_reads_illumina.py has already written it into the read name, so the whole
library can be aligned ONCE and split afterwards. For a 24-plex that is one STAR pass instead of 24.

    demux_qname_illumina.py <prepped.fastq[.gz]> <star.bam> <outdir> --run-json <run.json> --barcode <bc>
    demux_qname_illumina.py <prepped.fastq[.gz]> - <outdir> --run-json <run.json> [--barcode <bc>] [--r1 <mate_R1.fq[.gz]>]
        (BAM '-' = prepped-only: split the reads BEFORE alignment, for the per-sample mode that
         aligns every sample on its own; the stats then carry no bam_records fields.
         --r1: the R1 mate FASTQ prep wrote in the same order; split in lockstep into
         mateR1_<sample>.fastq, names checked record by record. An empty mate file (poly-dT
         chemistry) gives empty mateR1 files, which the alignment step treats as single-end.)
    demux_qname_illumina.py merge-stats --out merged.json chunk1.json chunk2.json ...
        (sum the stats of chunked demuxes; counts add, percentages are recomputed)

Writes, for every bobcode in the run's tso_species_map (empty files for absent ones):
    <outdir>/samples/<sample>/prepped_<sample>.fastq
    <outdir>/samples/<sample>/<sample>_star_combined.bam
and <outdir>/<bc>_demux_stats.json in the schema per_sample_table.py already reads.
Sample names are sanitised exactly as per_sample_table.fname does, so the BAM stems match the
table's lookup by construction.
"""
import argparse, gzip, json, os, subprocess, sys, collections

def safe(name):
    keep = ''.join(ch if (ch.isalnum() or ch in '-_.') else '_' for ch in name).strip('_')
    while '__' in keep:
        keep = keep.replace('__', '_')
    return keep

def opener(p): return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)

def code_of(qname):
    """<origid>_<umi1>_<umi2>_<bobcode>_<grun>_<T|F>  -> bobcode, or None."""
    p = qname.split()[0].split('_')
    return p[-3] if len(p) >= 6 else None

def merge_stats(parts):
    """One demux_stats dict from per-chunk dicts: every count adds, percentages are
    recomputed with the formulas below, the label/species fields must agree."""
    base = parts[0]
    out = {k: base[k] for k in ('ont_barcode', 'architecture', 'platform', 'split_on')}
    for k in ('reads_total', 'reads_assigned', 'reads_unclassified', 'reads_undeclared',
              'bam_records_total', 'bam_records_unassigned', 'mates_assigned'):
        if any(k in q for q in parts):
            out[k] = sum(int(q.get(k, 0)) for q in parts)
    nt, na = out['reads_total'], out['reads_assigned']
    out['assigned_pct'] = round(100.0 * na / nt, 3) if nt else 0.0
    per = {}
    for q in parts:
        for c, r in q['per_bobcode'].items():
            e = per.setdefault(c, {'sample': r['sample'], 'species': r['species'], 'reads': 0})
            if (e['sample'], e['species']) != (r['sample'], r['species']):
                sys.exit(f'demux merge-stats: chunks disagree on {c}: {e} vs {r}')
            e['reads'] += int(r.get('reads', 0))
            if 'bam_records' in r:
                e['bam_records'] = e.get('bam_records', 0) + int(r['bam_records'])
    for c, e in per.items():
        e['pct_of_total'] = round(100.0 * e['reads'] / nt, 4) if nt else 0.0
        e['pct_of_assigned'] = round(100.0 * e['reads'] / na, 4) if na else 0.0
    out['per_bobcode'] = {c: per[c] for c in sorted(per, key=lambda k: -per[k]['reads'])}
    und = collections.Counter()
    for q in parts:
        und.update({k: int(v) for k, v in (q.get('undeclared') or {}).items()})
    out['undeclared'] = dict(und.most_common(20))
    return out


if len(sys.argv) > 1 and sys.argv[1] == 'merge-stats':
    mp = argparse.ArgumentParser(prog='demux_qname_illumina.py merge-stats')
    mp.add_argument('parts', nargs='+'); mp.add_argument('--out', required=True)
    ma = mp.parse_args(sys.argv[2:])
    merged = merge_stats([json.load(open(x)) for x in ma.parts])
    json.dump(merged, open(ma.out, 'w'), indent=1)
    sys.stderr.write(f'demux merge-stats: {len(ma.parts)} chunk(s) -> {ma.out}: '
                     f'{merged["reads_total"]:,} reads, {merged["reads_assigned"]:,} assigned\n')
    sys.exit(0)

ap = argparse.ArgumentParser()
ap.add_argument('fastq'); ap.add_argument('bam'); ap.add_argument('outdir')
ap.add_argument('--run-json', required=True); ap.add_argument('--barcode', default='')
ap.add_argument('--r1', default='', help='R1 mate FASTQ from prep --r1-out (same order as the prepped reads)')
a = ap.parse_args()
cfg = json.load(open(a.run_json))
species = {k.upper(): v for k, v in (cfg.get('tso_species_map') or {}).items()}
labels = {k.upper(): v for k, v in (cfg.get('bobcode_labels') or {}).items()}
if not species:
    sys.exit('demux_qname: run json has no tso_species_map')
sample_of = {c: safe(labels.get(c, c)) for c in species}
if len(set(sample_of.values())) != len(sample_of):
    sys.exit(f'demux_qname: two bobcodes sanitise to the same sample name: {sample_of}')

sdir = os.path.join(a.outdir, 'samples')
for s in sample_of.values():
    os.makedirs(os.path.join(sdir, s), exist_ok=True)
fq = {c: open(os.path.join(sdir, s, f'prepped_{s}.fastq'), 'w') for c, s in sample_of.items()}
fq1 = {c: open(os.path.join(sdir, s, f'mateR1_{s}.fastq'), 'w') for c, s in sample_of.items()} if a.r1 else {}
r1h = opener(a.r1) if a.r1 else None
counts = collections.Counter(); undeclared = collections.Counter()
n_total = n_assigned = n_mates = 0
with opener(a.fastq) as fh:
    while True:
        h = fh.readline()
        if not h: break
        s = fh.readline(); p = fh.readline(); q = fh.readline()
        n_total += 1
        m = None
        if r1h is not None:
            h1 = r1h.readline()
            if h1:
                m = (h1, r1h.readline(), r1h.readline(), r1h.readline())
                if h1.split()[0] != h.split()[0]:
                    sys.exit(f'demux_qname: R1 mate out of step at read {n_total}: {h.split()[0]} vs {h1.split()[0]}')
            elif n_mates or n_total == 1:
                # an empty mate file is the single-end chemistry; a mate file that ENDS early is a bug
                if n_mates:
                    sys.exit(f'demux_qname: R1 mate FASTQ ended after {n_mates} reads, prepped has more')
                r1h.close(); r1h = None
        c = code_of(h[1:])
        if c in fq:
            fq[c].write(h); fq[c].write(s); fq[c].write(p); fq[c].write(q)
            counts[c] += 1; n_assigned += 1
            if m:
                fq1[c].write(''.join(m)); n_mates += 1
        elif c:
            undeclared[c] += 1
        else:
            counts['unclassified'] += 1
if r1h is not None:
    if r1h.readline():
        sys.exit('demux_qname: R1 mate FASTQ has more reads than the prepped FASTQ')
    r1h.close()
for f in list(fq.values()) + list(fq1.values()): f.close()

# BAM: one pass, one SAM writer per sample, then samtools converts each. Skipped in the
# prepped-only mode ('-'), where every sample is aligned on its own afterwards.
PREPPED_ONLY = (a.bam == '-')
bam_recs = collections.Counter(); bam_total = bam_unassigned = 0
if not PREPPED_ONLY:
  hdr = subprocess.run(['samtools', 'view', '-H', a.bam], capture_output=True, text=True, check=True).stdout
  sam = {c: open(os.path.join(sdir, s, f'{s}.sam'), 'w') for c, s in sample_of.items()}
  for f in sam.values(): f.write(hdr)
  rd = subprocess.Popen(['samtools', 'view', a.bam], stdout=subprocess.PIPE, text=True)
  for line in rd.stdout:
      bam_total += 1
      c = code_of(line[:line.index('\t')])
      if c in sam:
          sam[c].write(line); bam_recs[c] += 1
      else:
          bam_unassigned += 1
  rd.wait()
  for c, f in sam.items():
      f.close()
      s = sample_of[c]; sp = os.path.join(sdir, s, f'{s}.sam'); bp = os.path.join(sdir, s, f'{s}_star_combined.bam')
      subprocess.run(['samtools', 'view', '-b', '-o', bp, sp], check=True); os.remove(sp)

stats = {
    'ont_barcode': a.barcode, 'architecture': cfg.get('tso_arch', ''),
    'platform': cfg.get('platform', 'illumina'),
    'split_on': ('qname bobcode before alignment (prepped reads; every sample aligned on its own)'
                 if PREPPED_ONLY else 'qname bobcode after one pooled alignment'),
    'reads_total': n_total, 'reads_assigned': n_assigned,
    **({'mates_assigned': n_mates} if a.r1 else {}),
    'reads_unclassified': counts.get('unclassified', 0),
    'reads_undeclared': sum(undeclared.values()),
    **({} if PREPPED_ONLY else {'bam_records_total': bam_total, 'bam_records_unassigned': bam_unassigned}),
    'assigned_pct': round(100.0 * n_assigned / n_total, 3) if n_total else 0.0,
    'per_bobcode': {
        c: {'sample': labels.get(c, c), 'species': species.get(c, ''), 'reads': counts.get(c, 0),
            **({} if PREPPED_ONLY else {'bam_records': bam_recs.get(c, 0)}),
            'pct_of_total': round(100.0 * counts.get(c, 0) / n_total, 4) if n_total else 0.0,
            'pct_of_assigned': round(100.0 * counts.get(c, 0) / n_assigned, 4) if n_assigned else 0.0}
        for c in sorted(species, key=lambda k: -counts.get(k, 0))},
    'undeclared': dict(undeclared.most_common(20)),
}
json.dump(stats, open(os.path.join(a.outdir, f'{a.barcode or "library"}_demux_stats.json'), 'w'), indent=1)
sys.stderr.write(f'demux_qname: {n_total:,} reads, {n_assigned:,} assigned ({stats["assigned_pct"]}%) '
                 f'across {len(species)} samples; undeclared {sum(undeclared.values()):,}\n')

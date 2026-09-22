#!/usr/bin/env python3
"""Library metrics for the preprint: read length and insert length.
insert: library_metrics.py insert <sample> <pairs.bam> <out.json>
   Deduplicated BOBseq pair BAM (both mates): for every pair whose two mates are unique (MAPQ 255) primary alignments on the same
   canonical protein-coding transcript, the insert length on the mature transcript = |t(mate 1 5' end) - t(mate 2 5' end)| + 1
   (each mate's 5' end is one end of the cDNA fragment; exon-aware, so introns do not inflate it); the genomic template length (STAR
   TLEN) is kept beside it. Output: histogram (1-nt bins to 3,000) and quantiles per well.
readlen: library_metrics.py readlen <set>  -> aligned read length per method (and per mate for BOBseq) from positions_canonical npz.
"""

from settings import COVERAGE
import sys, os, json, numpy as np

A = COVERAGE
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
mode = sys.argv[1]
if mode == 'insert':
    import pysam, coverage_lib as cl

    sample, bam, out = sys.argv[2:5]
    tx, idx = cl.build()
    H = {}
    ins = []
    tl = []
    n_pairs = n_used = 0
    for a in pysam.AlignmentFile(bam, 'rb').fetch(until_eof=True):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality != 255:
            continue
        if not a.reference_name.startswith(cl.PREFIX):
            continue  # combined-genome BAM: the transcript index uses bare human chromosome names
        q = a.query_name
        p5 = a.reference_end - 1 if a.is_reverse else a.reference_start
        r = cl.position(idx, a.reference_name[6:], p5)
        if r is None:
            continue
        g, t, L, gs = r
        m = H.pop(q, None)
        if m is None:
            H[q] = (g, t, abs(a.template_length))
            continue
        n_pairs += 1
        if m[0] == g:
            ins.append(abs(m[1] - t) + 1)
            tl.append(m[2])
            n_used += 1
    ins = np.array(ins, int)
    tl = np.array(tl, int)
    h = np.bincount(np.minimum(ins, 3000), minlength=3001) if len(ins) else np.zeros(3001, int)
    if not len(ins):
        print(sample, 'no pairs positioned')
        sys.exit(0)
    json.dump(
        {
            'sample': sample,
            'pairs_both_mates_unique_on_canonical': n_pairs,
            'pairs_same_transcript': n_used,
            'insert_median': float(np.median(ins)),
            'insert_mean': float(ins.mean()),
            'insert_quantiles': {str(p): float(np.percentile(ins, p)) for p in (5, 10, 25, 50, 75, 90, 95)},
            'tlen_genomic_median': float(np.median(tl)),
            'hist_1nt_to_3000': h.tolist(),
        },
        open(out, 'w'),
    )
    print(
        sample,
        f'pairs {n_pairs:,}, same transcript {n_used:,}, insert median {np.median(ins):.0f} nt (p10 {np.percentile(ins, 10):.0f}, p90 {np.percentile(ins, 90):.0f}); genomic TLEN median {np.median(tl):.0f}',
    )
elif mode == 'readlen':
    from palette import METH

    SET = sys.argv[2]
    sfx = '' if SET == 'native' else '_50nt'
    out = {}
    S = [l.rstrip('\n').split('\t') for l in open(f'{A}/samples{sfx}.tsv') if l.strip()]
    for lab, bam, m in S:
        z = np.load(f'{A}/positions_canonical{sfx}/{lab}.npz')
        al = z['alen'].astype(int)
        mate = z['mate']
        for key, sel in ((m, np.ones(len(al), bool)),) + (
            ((f'{m} mate 1', mate == 1), (f'{m} mate 2', mate == 2)) if (mate > 0).any() else ()
        ):
            v = al[sel]
            if not len(v):
                continue
            d = out.setdefault(key, {'hist': np.zeros(302, int), 'n': 0})
            d['hist'] += np.bincount(np.minimum(v, 301), minlength=302)
            d['n'] += len(v)
    with open(
        f'{A}/aligned_bases_{SET}.tsv', 'w'
    ) as f:  # per sample: reads, mean aligned length, total aligned bases, fragments (BOBseq: mate-1 records = one per pair), bases per fragment
        f.write('method\tsample\treads\tmean_aligned_length\ttotal_aligned_bases\tfragments\tbases_per_fragment\n')
        for lab, bam, m in S:
            z = np.load(f'{A}/positions_canonical{sfx}/{lab}.npz')
            al = z['alen'].astype(np.int64)
            mate = z['mate']
            nf = int((mate == 1).sum()) if (mate > 0).any() else len(al)
            f.write(f'{m}\t{lab}\t{len(al)}\t{al.mean():.3f}\t{al.sum()}\t{nf}\t{al.sum() / max(nf, 1):.3f}\n')
    res = {}
    for k, d in out.items():
        h = d['hist']
        cum = np.cumsum(h) / h.sum()
        q = lambda p: int(np.searchsorted(cum, p))
        res[k] = {
            'n_reads': int(d['n']),
            'median': q(0.5),
            'p10': q(0.1),
            'p90': q(0.9),
            'mean': float((np.arange(302) * h).sum() / h.sum()),
            'hist_0_to_301': h.tolist(),
        }
        print(
            f"{k:22s} n {d['n']:>11,}  aligned length median {res[k]['median']} (p10 {res[k]['p10']}, p90 {res[k]['p90']}, mean {res[k]['mean']:.1f})"
        )
    json.dump(res, open(f'{A}/read_length_{SET}.json', 'w'))

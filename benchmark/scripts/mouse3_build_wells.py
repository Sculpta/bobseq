#!/usr/bin/env python3
"""Mouse wells of the 24-plex (RAW 264.7 controls 1-3) for the per-sample supplement. From each raw pair BAM (pre-dedup):
<arm>/mouse3/<well>/Aligned.out.bam   mate-1 records, read names '<id>_<bobcode>_<umi>' (what the native arm BAM holds; read_table.py input), as supp21_build_ris500.py
per_sample_bams/bobseq_pe_native_mouse/<well>.bam   primary mapped records (-F 0x904) on MOUSE contigs, bare contig names, sorted + indexed (the mouse
                                      counterpart of bobseq_pe_native; Picard, positions and the composition rules read it)
<arm>/mouse3/human_contig_records.json  mate-1 primary mapped records on HUMAN contigs, split by rDNA locus overlap: GRCm39 carries almost no rDNA, so STAR
                                      places mouse rRNA on the human rDNA copies (rdna_loci.bed header); those records are rRNA of the mouse well.
"""

from settings import WORK, PROJECT, RDNA_BED
import os, sys, json, collections, pysam

P = WORK
PR = PROJECT
A = f'{P}/benchmark_uniform/bob57_24plex_pe_native/mouse3'
PB = f'{P}/benchmark_uniform/per_sample_bams/bobseq_pe_native_mouse'
WELLS = ['RAW_control_1', 'RAW_control_2', 'RAW_control_3']
os.makedirs(PB, exist_ok=True)
rd = collections.defaultdict(list)
for l in open(RDNA_BED):
    if l[0] == '#' or not l.strip() or l.startswith(('track', 'browser')):
        continue
    c, s, e = l.split()[:3]
    rd[c].append((int(s), int(e)))


def in_rdna(c, s, e):
    return any(a < e and b > s for a, b in rd.get(c, ()))


hum = {}
for w in WELLS:
    src = f'{PR}/data/bams/pairs_raw/{w}_pairs_raw.bam'
    os.makedirs(f'{A}/{w}', exist_ok=True)
    b = pysam.AlignmentFile(src, 'rb')
    out = pysam.AlignmentFile(f'{A}/{w}/Aligned.out.bam', 'wb', template=b)
    refs = [(i, n) for i, n in enumerate(b.references) if n.startswith('MOUSE_')]
    remap = {i: k for k, (i, n) in enumerate(refs)}
    hdr = {'HD': {'VN': '1.4', 'SO': 'unsorted'}, 'SQ': [{'SN': n[6:], 'LN': b.lengths[i]} for i, n in refs]}
    tmp = f'{PB}/{w}.unsorted.bam'
    mb = pysam.AlignmentFile(tmp, 'wb', header=hdr)
    n1 = n2 = nm = 0
    h = collections.Counter()
    for a in b.fetch(until_eof=True):
        if not (a.flag & 0x904) and a.reference_id in remap:
            d = a.to_dict()
            d['ref_name'] = a.reference_name[6:]
            nr = a.next_reference_name
            d['next_ref_name'] = (
                '=' if a.next_reference_id == a.reference_id else (nr[6:] if nr and nr.startswith('MOUSE_') else '*')
            )
            if d['next_ref_name'] == '*':
                d['next_ref_pos'] = '0'
                d['length'] = '0'
            mb.write(pysam.AlignedSegment.from_dict(d, mb.header))
            nm += 1
        if a.is_read2:
            n2 += 1
            continue
        if not (a.flag & 0x904) and a.reference_name.startswith('HUMAN_'):
            h['rdna' if in_rdna(a.reference_name, a.reference_start, a.reference_end) else 'other'] += 1
        rid, rest = a.query_name.split('__', 1)
        umi, bobcode = rest.split('_')[:2]
        a.query_name = f'{rid}_{bobcode}_{umi}'
        n1 += 1
        out.write(a)
    out.close()
    mb.close()
    pysam.sort('-@', '4', '-o', f'{PB}/{w}.bam', tmp)
    pysam.index(f'{PB}/{w}.bam')
    os.remove(tmp)
    hum[w] = dict(h)
    print(
        w,
        f'mate-1 records {n1:,} (mate-2 skipped {n2:,}); mouse-contig primary records {nm:,}; mate-1 primary on human contigs: rDNA loci {h["rdna"]:,}, other {h["other"]:,}',
        flush=True,
    )
json.dump(hum, open(f'{A}/human_contig_records.json', 'w'), indent=1)
print('MOUSE3 BUILD DONE')

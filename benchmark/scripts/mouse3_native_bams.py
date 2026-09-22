#!/usr/bin/env python3
"""Mouse wells (RAW 264.7 controls 1-3): Picard CollectRnaSeqMetrics, rule-based read composition, canonical-transcript positions and insert length,
with the same code as the human wells and BM_SPECIES_PREFIX=MOUSE_ (set here before the imports). Mouse refFlat = protein-coding transcripts, no
ribosomal-protein genes, no MT, bare contig names (the rule of picard_profile.py). Composition: the rules on the mouse-contig BAM, plus the well's
mate-1 primary records on HUMAN rDNA loci counted as rRNA (GRCm39 has almost no rDNA; STAR places mouse rRNA on the human copies); the few other
human-contig records are left out, as mouse-contig records are for the human wells."""

from settings import CODE, WORK, PROJECT, COVERAGE, GTF, REFFLAT
import os, re, sys, json, subprocess, collections

os.environ['BM_SPECIES_PREFIX'] = 'MOUSE_'
os.environ['BM_COVSET'] = 'native'
H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, H)
P = WORK
PR = PROJECT
PB = f'{P}/benchmark_uniform/per_sample_bams'
A = COVERAGE
MET = f'{P}/preprint_figures_native/coverage_architecture/picard_metrics/native_mouse'
os.makedirs(MET, exist_ok=True)
REFFLAT = f'{A}/picard_refflat_pc_noRP_noMT_bare_mouse.txt'
PICARD = 'picard'
RIBO = re.compile(r'^(RP[LS]|RPLP|MRP[LS])[^K]*$', re.I)
WELLS = ['RAW_control_1', 'RAW_control_2', 'RAW_control_3']
step = sys.argv[1]
if step == 'picard':
    if not os.path.exists(REFFLAT):
        attr = re.compile(r'(\S+) "([^"]*)"')
        tx = {}
        ex = collections.defaultdict(list)
        cds = collections.defaultdict(list)
        for l in open(GTF):
            if l[0] == '#' or not l.startswith('MOUSE_'):
                continue
            f = l.rstrip('\n').split('\t')
            if f[2] not in ('exon', 'CDS'):
                continue
            a = dict(attr.findall(f[8]))
            if a.get('transcript_biotype') != 'protein_coding':
                continue
            chrom = f[0][6:]
            if chrom == 'MT':
                continue
            g = a.get('gene_name', a.get('gene_id', ''))
            t = a['transcript_id']
            if RIBO.match(g):
                continue
            tx[t] = (g, chrom, f[6])
            (ex if f[2] == 'exon' else cds)[t].append((int(f[3]) - 1, int(f[4])))
        with open(REFFLAT, 'w') as o:
            for t, (g, chrom, strand) in tx.items():
                e = sorted(ex[t])
                c = sorted(cds[t])
                cs, ce = (c[0][0], c[-1][1]) if c else (e[-1][1], e[-1][1])
                o.write(
                    f'{g}\t{t}\t{chrom}\t{strand}\t{e[0][0]}\t{e[-1][1]}\t{cs}\t{ce}\t{len(e)}\t{",".join(str(s) for s, _ in e)},\t{",".join(str(x) for _, x in e)},\n'
                )
        print('mouse refFlat written:', sum(1 for _ in open(REFFLAT)), 'transcripts', flush=True)
    procs = []
    for n in WELLS:
        cmd = [
            PICARD,
            'CollectRnaSeqMetrics',
            '-I',
            f'{PB}/bobseq_pe_native_mouse/{n}.bam',
            '-O',
            f'{MET}/{n}.RNA_Metrics.txt',
            '--REF_FLAT',
            REFFLAT,
            '--STRAND_SPECIFICITY',
            'NONE',
            '--MINIMUM_LENGTH',
            '1000',
            '--VALIDATION_STRINGENCY',
            'LENIENT',
        ]
        procs.append((n, subprocess.Popen(cmd, stdout=open(f'{MET}/{n}.log', 'w'), stderr=subprocess.STDOUT)))
    for n, p in procs:
        print(n, 'picard', 'done' if p.wait() == 0 else 'FAILED', flush=True)
elif step == 'composition':
    import composition_rules as cr

    assert cr.PREFIX == 'MOUSE_' and cr.IDXF.endswith('_mouse.npz')
    hum = json.load(open(f'{P}/benchmark_uniform/bob57_24plex_pe_native/mouse3/human_contig_records.json'))
    res = {}
    for n in WELLS:
        _, name, cnt = cr.classify_bam(('BOBseq', n, f'bobseq_pe_native_mouse/{n}.bam', True))
        cnt = dict(cnt)
        cnt['rRNA_on_mouse_contigs'] = cnt.get('rRNA', 0)
        cnt['rRNA'] = cnt.get('rRNA', 0) + hum[n].get('rdna', 0)
        cnt['human_contig_non_rdna_excluded'] = hum[n].get('other', 0)
        res[name] = cnt
        print(name, cnt, flush=True)
    os.makedirs(f'{PR}/results/supplement_21', exist_ok=True)
    json.dump(res, open(f'{PR}/results/supplement_21/composition_rules_mouse3_native.json', 'w'), indent=1)
elif step == 'positions':
    procs = [
        (
            k,
            subprocess.Popen(
                [
                    'python3',
                    f'{CODE}/build_positions.py',
                    f'bobseq-raw-control-{k}',
                    f'{PB}/bobseq_pe_native_mouse/RAW_control_{k}.bam',
                    f'{A}/positions_canonical/bobseq-raw-control-{k}.npz',
                ],
                stdout=open(f'{A}/positions_canonical/bobseq-raw-control-{k}.log', 'w'),
                stderr=subprocess.STDOUT,
            ),
        )
        for k in (1, 2, 3)
    ]
    for k, p in procs:
        print(
            'positions',
            k,
            'exit',
            p.wait(),
            open(f'{A}/positions_canonical/bobseq-raw-control-{k}.log').read().strip()[-160:],
            flush=True,
        )
elif step == 'insert':
    O = f'{PR}/results/library_metrics/insert_per_well'
    procs = [
        (
            n,
            subprocess.Popen(
                [
                    'python3',
                    f'{CODE}/library_metrics.py',
                    'insert',
                    n,
                    f'{PR}/data/bams/pipeline_dedup_pe/dedup_pairs_6NHH/{n}_star_pairs.bam',
                    f'{O}/{n}.json',
                ],
                stdout=open(f'{O}/{n}.log', 'w'),
                stderr=subprocess.STDOUT,
            ),
        )
        for n in WELLS
    ]
    for n, p in procs:
        print('insert', n, 'exit', p.wait(), open(f'{O}/{n}.log').read().strip()[-200:], flush=True)
print(f'MOUSE3 {step} DONE')

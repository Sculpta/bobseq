#!/usr/bin/env python3
"""Per-SAMPLE PDF report for a pooled multiplex run.

The library report (BarcodeSum) has one column per ONT barcode, which for a pooled run is
one column per PCR primer set -- not per sample. This renders the same run with the SAMPLE
as the unit, so the experiment is readable: conditions, replicates, and the two accuracy
views side by side. Platform-aware: on Illumina (run json platform illumina_*) the unit is
the bobcode read from R2 in one library, the funnel starts from the raw pairs (prep stats),
the duplicate rate is definition D from the dedup sidecars, and the ONT-only structural
pages are not drawn.

Inputs are the two published TSVs, so this cannot invent numbers: everything shown is
computed upstream by per_sample_table.py and per_sample_composition.py.

CLI: per_sample_report.py --per-sample <tsv> --composition <tsv> --sheet <tsv> --out <pdf>
"""
import argparse, collections, json, os, statistics, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 9, 'axes.linewidth': 0.8})
HUM, MOU, MUT, INK = '#C2268F', '#3D9B35', '#5E6B68', '#101615'
CATS = [('mRNA (protein-coding)', '#C2268F'), ('Ribosomal-protein mRNA', '#F0A8D8'),
        ('Other ncRNA (lncRNA/NMD/ret-intron/pseudo)', '#2E7FC2'), ('rRNA', '#3D9B35'),
        ('Mitochondrial', '#E8820C'), ('pre-mRNA / intronic', '#6A5ACD'),
        ('Intergenic / genomic', '#A08A00'), ('Unmapped', '#B8BEBC')]


def fname(s):
    """Same sanitiser the demux step uses for filenames and column names."""
    k = ''.join(c if (c.isalnum() or c in '-_.') else '_' for c in s)
    while '__' in k:
        k = k.replace('__', '_')
    return k.strip('_')


def read_tsv(p):
    rows = [l.rstrip('\n').split('\t') for l in open(p) if not l.startswith('#')]
    return [dict(zip(rows[0], r)) for r in rows[1:]]


def _sankey(ax, columns, flows, INK, MUT, gap_frac=0.035, node_w=0.012, label_fs=7.6):
    """Minimal Sankey: columns = list of lists of (node_id, label, value, colour); flows = list of
    (src_id, dst_id, value, colour). Node heights are proportional to value within a shared scale;
    each flow occupies consecutive segments of its source and target node (in list order)."""
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch, Rectangle
    ncol = len(columns)
    tot = max(sum(v for _, _, v, _ in col) for col in columns)
    scale = 1.0 / (tot * (1 + gap_frac * max(len(col) for col in columns)))
    pos = {}                                  # node_id -> (x, y0, y1, colour)
    for ci, col in enumerate(columns):
        x = ci / (ncol - 1) * (1 - node_w) if ncol > 1 else 0
        used = sum(v for _, _, v, _ in col) * scale + gap_frac * (len(col) - 1)
        y = 1 - (1 - used) / 2                # centre the column vertically
        for nid, lab, v, c in col:
            h = v * scale
            pos[nid] = (x, y - h, y, c)
            ax.add_patch(Rectangle((x, y - h), node_w, h, color=c, lw=0, zorder=4))
            y -= h + gap_frac
    src_off = {nid: 0.0 for nid in pos}; dst_off = {nid: 0.0 for nid in pos}
    for src, dst, v, c in flows:
        if v <= 0 or src not in pos or dst not in pos:
            continue
        h = v * scale
        xs, s0, s1, _ = pos[src]; xd, d0, d1, _ = pos[dst]
        ys_top = s1 - src_off[src]; yd_top = d1 - dst_off[dst]
        src_off[src] += h; dst_off[dst] += h
        x0, x1 = xs + node_w, xd
        cx0, cx1 = x0 + (x1 - x0) * 0.45, x0 + (x1 - x0) * 0.55
        verts = [(x0, ys_top), (cx0, ys_top), (cx1, yd_top), (x1, yd_top), (x1, yd_top - h), (cx1, yd_top - h), (cx0, ys_top - h), (x0, ys_top - h), (x0, ys_top)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
        ax.add_patch(PathPatch(Path(verts, codes), facecolor=c, edgecolor='none', alpha=0.45, zorder=2))
    # labels to the right of each node; within a column, labels are pushed apart so none overlap
    # (a thin leader joins a moved label to its node)
    min_gap = label_fs * 0.0026
    def _gap(lab): return min_gap * (1 + lab.count('\n'))
    for ci, col in enumerate(columns):
        items = [(nid, lab, (pos[nid][1] + pos[nid][2]) / 2) for nid, lab, v, c in col if lab]
        items.sort(key=lambda t: -t[2])            # top to bottom
        ys = [t[2] for t in items]
        for i in range(1, len(ys)):                 # push down anything too close to the one above
            g = (_gap(items[i - 1][1]) + _gap(items[i][1])) / 2
            if ys[i - 1] - ys[i] < g:
                ys[i] = ys[i - 1] - g
        shift = max(0.0, -0.01 - ys[-1]) if ys else 0.0   # keep the stack inside the axes
        ys = [y + shift for y in ys]
        for (nid, lab, y_node), y in zip(items, ys):
            x = pos[nid][0]
            if abs(y - y_node) > 1e-6:
                ax.plot([x + node_w, x + node_w + 0.006], [y_node, y], color=MUT, lw=0.5, zorder=5)
            ax.text(x + node_w + 0.007, y, lab, va='center', ha='left', fontsize=label_fs, color=INK, zorder=6,
                    bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.8))
    ax.set_xlim(-0.01, 1.3); ax.set_ylim(-0.03, 1.02); ax.axis('off')
    return pos


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--per-sample', required=True)
    ap.add_argument('--composition', required=True)
    ap.add_argument('--demux-stats', nargs='*', default=None,
                    help='demux_stats.json per ONT barcode. Supplies the run-level '
                         'funnel: what fraction of raw reads reached a sample at all, '
                         'and what was dropped on the way.')
    ap.add_argument('--master-summary', default=None,
                    help='sample-scope master_summary.tsv from analyze_samples. Supplies '
                         'the structural/chemistry metrics the analysis computes '
                         'and this script does not: TSO backbone, chimeras, G-run, polyT, '
                         'empty products, the duplicate funnel, insert sizes.')
    ap.add_argument('--run-json', required=True,
                    help='plate order comes from bobcode_labels, which preserves the '
                         'order the samples were declared in; no external sheet needed')
    ap.add_argument('--assets', default=None,
                    help='dir holding 24plex_bobcode_set_v1.tsv, for the BCnn display id')
    ap.add_argument('--run-name', default='')
    ap.add_argument('--prep-stats', nargs='*', default=None,
                    help='Illumina: library prep_stats_*.json (raw pairs, bobcode called, usable '
                         'after the poly(A) trim). Without it the funnel starts at the demux '
                         'total, which on Illumina is already the prepped reads.')
    ap.add_argument('--sidecars', nargs='*', default=None,
                    help='Illumina: per-sample *_dup_prefilter.json (mapped before dedup, kept '
                         'after dedup), so duplicates removed by dedup are not shown as unmapped')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    # Illumina prep and dedup accounting (both optional; the ONT module passes neither).
    # Prep drops reads without a bobcode BEFORE demux, so the demux total is not the raw
    # count and its unclassified count is zero by construction; the prep stats hold the
    # real losses. The per-sample BAMs are post-dedup, so the composition's 'Unmapped'
    # is everything not kept (duplicates + unmapped); the sidecars separate the two.
    pr = {'n': 0, 'bobcode_called': 0, 'no_bobcode': 0, 'too_short_after_trim': 0, 'written': 0}
    for f in (args.prep_stats or []):
        if os.path.exists(f):
            d = json.load(open(f))
            for k in pr:
                pr[k] += int(d.get(k, 0) or 0)
    sc = {'mapped_with_umi': 0, 'unmapped_with_umi': 0, 'no_umi_call': 0, 'dedup_kept_aligned': 0}
    scs = {}                                   # per sample, keyed by the file's <sample>_ prefix
    for f in (args.sidecars or []):
        if os.path.exists(f):
            d = json.load(open(f))
            if d.get('empty_sample'):
                continue
            scs[os.path.basename(f).replace('_dup_prefilter.json', '')] = d
            for k in sc:
                sc[k] += int(d.get(k, 0) or 0)

    ps = {r['sample']: r for r in read_tsv(args.per_sample)}
    comp = {r['sample']: r for r in read_tsv(args.composition)}

    # master_summary is metric-per-ROW, sample-per-COLUMN; invert to {sample: {metric: v}}
    # run-level funnel, so the report accounts for every raw read rather than starting
    # at 'assigned' and leaving the reader to wonder what happened to the rest.
    dx = {'total': 0, 'assigned': 0, 'undeclared': 0, 'unclassified': 0, 'codes': {}}
    for f in (args.demux_stats or []):
        if not os.path.exists(f):
            continue
        d = json.load(open(f))
        dx['total'] += d.get('reads_total', 0)
        dx['assigned'] += d.get('reads_assigned', 0)
        dx['undeclared'] += d.get('reads_undeclared', 0)
        dx['unclassified'] += d.get('reads_unclassified', 0)
        for k, v in (d.get('undeclared_bobcodes') or {}).items():
            dx['codes'][k] = dx['codes'].get(k, 0) + v

    ms = {}
    if args.master_summary and os.path.exists(args.master_summary):
        rows = [l.rstrip('\n').split('\t') for l in open(args.master_summary)]
        hdr = rows[0]
        for j, col in enumerate(hdr[1:], 1):
            ms[col] = {r[0]: (r[j] if j < len(r) else '') for r in rows[1:] if r and r[0]}

    def msnum(sample, metric, default=None):
        """First number in a master_summary cell. Cells mix formats -- '3.67%',
        '247 / 219', '1,557 (35.5%)' -- so pull the leading numeric token and let the
        caller decide what it means."""
        import re as _re
        # master_summary columns use the SANITISED sample name (the demux safe():
        # every non [A-Za-z0-9-_.] becomes '_'), so 'CHX dose 1ug/mL-1' is stored as
        # 'CHX_dose_1ug_mL-1'. Substituting spaces alone would silently blank every
        # sample whose name also contains a slash.
        row = ms.get(sample) or ms.get(fname(sample)) or {}
        v = row.get(metric)
        if not v:
            return default
        m = _re.search(r'-?[\d,]+\.?\d*', v.replace(',', ''))
        return float(m.group()) if m else default

    # Plate order comes from the run config: bobcode_labels is written in declaration
    # order, which is the order the samples were laid out. Keeping the report inside the
    # pipeline means it cannot drift from the config the run actually used.
    import json as _json
    cfg = _json.load(open(args.run_json))
    illumina = str(cfg.get('platform', '')).lower().startswith('illumina')
    # Samples are named by their label, or by the code itself when the run declares no
    # labels (demux_qname / per_sample_table do the same); without this fallback an
    # unlabeled run would have an empty sample order and divide by zero below.
    lab_by_seq = cfg.get('bobcode_labels') or {seq: seq for seq in (cfg.get('tso_species_map') or {})}
    name_by_seq = {}
    if args.assets:
        reg = os.path.join(args.assets, '24plex_bobcode_set_v1.tsv')
        if os.path.exists(reg):
            for ln in open(reg):
                f = ln.rstrip('\n').split('\t')
                if f[0] != 'ID' and len(f) > 2:
                    name_by_seq[f[2]] = f[1].replace('BOB', 'BC')
    order = [nm for seq, nm in lab_by_seq.items() if nm in ps]
    code = {nm: name_by_seq.get(seq, '') for seq, nm in lab_by_seq.items() if nm in ps}
    labels = [f"{code[s]:<6} {s}" if code[s] else s for s in order]
    y = list(range(len(order)))[::-1]          # top-to-bottom in plate order
    sp = [ps[s]['species'] for s in order]
    col = [HUM if x == 'human' else MOU for x in sp]
    run = args.run_name or 'run'
    n_tot = sum(int(ps[s]['reads_assigned']) for s in order)

    # Per-stage run totals, so the funnel page and the page subtitles all quote the
    # same arithmetic. Every stage is a share of the RAW reads that entered the run.
    def _pct(s_, cat):
        return float(comp[s_][cat]) if s_ in comp and cat in comp[s_] else 0.0
    n_unmap = sum(int(ps[s]['reads_assigned']) * _pct(s, 'Unmapped') / 100.0 for s in order)
    n_map = n_tot - n_unmap
    n_mrna = sum(int(ps[s]['reads_assigned'])
                 * (_pct(s, 'mRNA (protein-coding)') + _pct(s, 'Ribosomal-protein mRNA'))
                 / 100.0 for s in order)
    n_rrna = sum(int(ps[s]['reads_assigned']) * _pct(s, 'rRNA') / 100.0 for s in order)
    raw = pr['n'] or dx['total'] or n_tot
    pc = lambda v: 100.0 * v / raw if raw else 0.0
    # what the composition's 'Unmapped' means depends on what the BAM holds
    unmapped_label = 'removed by dedup or unmapped' if sc['dedup_kept_aligned'] else 'unmapped'

    with PdfPages(args.out) as pdf:
        # ---- page 1: where the reads went ----------------------------------
        # The funnel exists so no stage is silently skipped: a reader can see that a
        # read either reached a sample and got classified, or was dropped, and where.
        # The first three rows are sequential stages; the last two are subsets of the
        # mapped reads, drawn lighter so nobody reads them as further attrition.
        undecl = (f" + {dx['undeclared']:,} bobcode not in this pool" if dx['undeclared'] else '')
        if pr['n']:
            # Illumina: prep accounting first, then (if the sidecars came along) mapping
            # before dedup and what dedup kept; rRNA/mRNA are subsets of the kept reads.
            stages = [('raw read pairs sequenced', raw, INK, 1.0),
                      ('carried a readable bobcode', pr['bobcode_called'], HUM, 1.0),
                      ('usable insert (>= 30 nt after poly(A) trim)', pr['written'], HUM, 0.8)]
            drops = ['', f"lost: {pr['no_bobcode']:,} no readable bobcode" + undecl,
                     f"lost: {pr['too_short_after_trim']:,} bobcode ran into poly(A)"]
            if sc['dedup_kept_aligned']:
                n_map = sc['dedup_kept_aligned']
                stages += [('mapped to the genome (before dedup)', sc['mapped_with_umi'], '#6E8FA6', 1.0),
                           ('kept after deduplication (definition D)', n_map, '#6E8FA6', 0.7)]
                drops += [f"lost: {sc['unmapped_with_umi'] + sc['no_umi_call']:,} unmapped",
                          f"lost: {sc['mapped_with_umi'] - n_map:,} duplicates removed"]
            else:
                stages += [('kept after deduplication and mapped', n_map, '#6E8FA6', 1.0)]
                drops += [f'lost: {n_unmap:,.0f} {unmapped_label}']
            subset_of = 'reads kept after deduplication'
        else:
            stages = [('raw reads sequenced', raw, INK, 1.0),
                      ('carried a readable bobcode', dx['assigned'] or n_tot, HUM, 1.0),
                      ('mapped to the genome', n_map, '#6E8FA6', 1.0)]
            drops = ['', f"lost: {dx['unclassified']:,} no readable bobcode" + undecl,
                     f'lost: {n_unmap:,.0f} {unmapped_label}']
            subset_of = 'mapped reads'
        stages += [('\u2514 of which rRNA', n_rrna, '#7FA98C', 0.45),
                   ('\u2514 of which mRNA (incl. RP)', n_mrna, '#B4508E', 0.45)]
        drops += ['', '']
        fig, ax = plt.subplots(figsize=(8.5, 3.4 + 0.4 * len(stages)))   # 5 stages = the ONT page
        yy = list(range(len(stages)))[::-1]
        def _wash(hexc, keep):
            r, g, b = (int(hexc[i:i + 2], 16) for i in (1, 3, 5))
            return '#%02X%02X%02X' % tuple(int(v * keep + 255 * (1 - keep))
                                           for v in (r, g, b))
        ax.barh(yy, [v for _, v, _, _ in stages],
                color=[c if al >= 1.0 else _wash(c, al) for _, _, c, al in stages],
                height=0.62, zorder=3)
        for j, (lab, v, _, _) in zip(yy, stages):
            ax.text(v + raw * 0.012, j, f'{v:,.0f}   {pc(v):.1f}%', va='center',
                    fontsize=9, color=INK, fontweight='bold')
        ax.set_yticks(yy)
        ax.set_yticklabels([l for l, _, _, _ in stages], fontsize=9.5)
        ax.set_xlim(0, raw * 1.32)
        ax.set_xlabel('reads')
        ax.set_title(f'Where the reads went\n{run}', loc='left', fontsize=12,
                     fontweight='bold', color=INK)
        for j, d in zip(yy, drops):
            if d:
                ax.text(raw * 0.012, j - 0.40, d, ha='left', va='center',
                        fontsize=7.5, color='#B04A3A', style='italic')
        ax.grid(True, axis='x', ls=':', lw=0.5, color='#C8CFCD', zorder=0)
        for s_ in ('top', 'right', 'left'):
            ax.spines[s_].set_visible(False)
        note = ('Percentages are of raw reads. The two lighter bars are subsets of the '
                f'{subset_of}, not further stages: '
                + (f'rRNA is {100.0*n_rrna/n_map:.0f}% and mRNA {100.0*n_mrna/n_map:.0f}% '
                   + ('of what mapped.' if not pr['n'] else f'of the {subset_of}.')
                   if n_map else 'nothing mapped.'))
        if dx['codes']:
            top = sorted(dx['codes'].items(), key=lambda kv: -kv[1])[:3]
            note += '\nBobcodes seen but not declared in this pool: ' + ', '.join(
                f'{k} ({v:,})' for k, v in top) + '  — contamination or misread.'
        fig.text(0.01, 0.005, note, fontsize=7.5, color=MUT)
        pdf.savefig(fig, bbox_inches='tight'); plt.close(fig)


        # ---- page 1b: read attrition (Sankey), Illumina only ------------------------
        # Every band is a count from the same files as the funnel: prep stats (raw, bobcode called
        # per code, no bobcode), demux stats (reads with an insert per sample), sidecars (mapped,
        # unmapped, kept after dedup) and the composition table (what the kept reads are).
        by_code = {}
        for f in (args.prep_stats or []):
            if os.path.exists(f):
                for k, v in (json.load(open(f)).get('by_bobcode') or {}).items():
                    by_code[k.upper()] = by_code.get(k.upper(), 0) + int(v)
        code_of = {}
        for f in (args.demux_stats or []):
            if os.path.exists(f):
                for k, v in (json.load(open(f)).get('per_bobcode') or {}).items():
                    code_of[v.get('sample', k)] = k.upper()
        if pr['n'] and sc['dedup_kept_aligned'] and by_code and code_of:
            raw = pr['n']
            def _k(s_): return s_.replace(' ', '_').replace('/', '_')
            sample_nodes, flows = [], []
            cols_sample = []
            for s_ in order:
                code = code_of.get(s_, ''); called = by_code.get(code, 0); prepped = int(ps[s_]['reads_assigned'])
                c = HUM if ps[s_]['species'] == 'human' else MOU
                cols_sample.append((f's:{s_}', f'{s_}  {called/1e6:.2f} M', called, c))
                flows.append(('raw', f's:{s_}', called, c))
                flows.append((f's:{s_}', 'prepped', prepped, c))
                flows.append((f's:{s_}', 'noins', max(called - prepped, 0), c))
            cols_sample.append(('nobc', f"no bobcode\n{pr['no_bobcode']/1e6:.2f} M ({100*pr['no_bobcode']/raw:.1f}%)", pr['no_bobcode'], '#C8CFCD'))
            flows.append(('raw', 'nobc', pr['no_bobcode'], '#C8CFCD'))
            prepped_t = sum(int(ps[s_]['reads_assigned']) for s_ in order); noins_t = sum(max(by_code.get(code_of.get(s_, ''), 0) - int(ps[s_]['reads_assigned']), 0) for s_ in order)
            mapped_t, unm_t, kept_t = sc['mapped_with_umi'], sc['unmapped_with_umi'] + sc['no_umi_call'], sc['dedup_kept_aligned']
            dup_t = max(mapped_t - kept_t, 0)
            flows += [('prepped', 'mapped', mapped_t, '#6E8FA6'), ('prepped', 'unmapped', unm_t, '#B0BCC0'),
                      ('mapped', 'kept', kept_t, '#6E8FA6'), ('mapped', 'dups', dup_t, '#B0BCC0')]
            # composition of the kept reads, summed over samples, scaled onto the kept node
            comp_counts = {}
            for s_ in order:
                n_ = int(ps[s_]['reads_assigned'])
                for cat, c in CATS:
                    if cat != 'Unmapped':
                        comp_counts[cat] = comp_counts.get(cat, 0.0) + n_ * _pct(s_, cat) / 100.0
            tot_c = sum(comp_counts.values()) or 1.0
            comp_cols = []
            for cat, c in CATS:
                if cat == 'Unmapped':
                    continue
                v = comp_counts[cat] / tot_c * kept_t
                short = cat.split(' (')[0]
                comp_cols.append((f'c:{cat}', f'{short}  {v/1e6:.2f} M ({100*v/raw:.1f}%)', v, c)); flows.append(('kept', f'c:{cat}', v, c))   # one line: these are the right-most labels
            def _lab(name, v): return f'{name}\n{v/1e6:.2f} M ({100*v/raw:.1f}%)'
            columns = [[('raw', _lab('raw pairs', raw), raw, INK)],
                       cols_sample,
                       [('prepped', _lab('prepped (insert >= 30 nt)', prepped_t), prepped_t, '#6E8FA6'), ('noins', _lab('too short after poly(A) trim', noins_t), noins_t, '#D9C8A6')],
                       [('mapped', _lab('mapped', mapped_t), mapped_t, '#6E8FA6'), ('unmapped', _lab('unmapped', unm_t), unm_t, '#B0BCC0')],
                       [('kept', _lab('kept after dedup', kept_t), kept_t, '#6E8FA6'), ('dups', _lab('duplicates removed', dup_t), dup_t, '#B0BCC0')],
                       comp_cols]
            fig, ax = plt.subplots(figsize=(13.5, 9.5))
            _sankey(ax, columns, flows, INK, MUT, gap_frac=0.006 if len(order) > 8 else 0.03, label_fs=6.6 if len(order) > 8 else 8)
            ax.set_title(f'Read attrition, raw pairs to classified kept reads\n{run}', loc='left', fontsize=12, fontweight='bold', color=INK)
            fig.text(0.01, 0.005, 'Band height = reads; percentages are of raw pairs. Samples are coloured by declared species (magenta human, green mouse). '
                     'Stages: bobcode called (prep) -> insert kept after the poly(A) trim (prepped, what demux sees) -> mapped by STAR -> kept after\n'
                     'definition-D deduplication -> composition of the kept reads (per-sample composition table, scaled onto the kept total). '
                     'Unmapped and duplicate reads are pooled over samples; per-sample values are in the tables of this report.',
                     fontsize=7.5, color=MUT)
            pdf.savefig(fig, bbox_inches='tight'); plt.close(fig)

        # ---- page 1c: the same picture in MOLECULES (duplicates already removed) ----------
        # Every band is a
        # definition-D molecule from the sidecars: per sample kept molecules (aligned + unmapped
        # kept) -> aligned / unmapped -> composition of the aligned molecules. No read-level
        # stage appears, so the picture is what the library delivers after deduplication.
        if illumina and scs and sc['dedup_kept_aligned']:
            mol_cols, flows2 = [], []
            tot_al = tot_un = 0
            for s_ in order:
                d = scs.get(s_.replace(' ', '_').replace('/', '_'), {})
                al = int(d.get('dedup_kept_aligned', 0) or 0); un = max(int(d.get('dedup_kept', 0) or 0) - al, 0)
                c = HUM if ps[s_]['species'] == 'human' else MOU
                mol_cols.append((f'm:{s_}', f'{s_}  {(al + un)/1e6:.2f} M', al + un, c))
                flows2.append((f'm:{s_}', 'al', al, c)); flows2.append((f'm:{s_}', 'un', un, c))
                tot_al += al; tot_un += un
            tot_m = tot_al + tot_un or 1
            comp_counts = {}
            for s_ in order:
                d = scs.get(s_.replace(' ', '_').replace('/', '_'), {}); al = int(d.get('dedup_kept_aligned', 0) or 0)
                for cat, c in CATS:
                    if cat != 'Unmapped':
                        comp_counts[cat] = comp_counts.get(cat, 0.0) + al * _pct(s_, cat) / 100.0
            tot_c = sum(comp_counts.values()) or 1.0
            comp_cols2 = []
            for cat, c in CATS:
                if cat == 'Unmapped':
                    continue
                v = comp_counts[cat] / tot_c * tot_al; short = cat.split(' (')[0]
                comp_cols2.append((f'cm:{cat}', f'{short}  {v/1e6:.2f} M ({100*v/tot_m:.1f}%)', v, c)); flows2.append(('al', f'cm:{cat}', v, c))
            columns2 = [mol_cols,
                        [('al', f'aligned molecules\n{tot_al/1e6:.2f} M ({100*tot_al/tot_m:.1f}%)', tot_al, '#6E8FA6'),
                         ('un', f'unmapped molecules\n{tot_un/1e6:.2f} M ({100*tot_un/tot_m:.1f}%)', tot_un, '#B0BCC0')],
                        comp_cols2]
            fig, ax = plt.subplots(figsize=(11, 9.5))
            _sankey(ax, columns2, flows2, INK, MUT, gap_frac=0.006 if len(order) > 8 else 0.03, label_fs=6.6 if len(order) > 8 else 8)
            ax.set_title(f'Molecules after deduplication, per sample to composition\n{run}', loc='left', fontsize=12, fontweight='bold', color=INK)
            fig.text(0.01, 0.005, f'Band height = definition-D molecules (one per bobcode + position + fragment length + UMI), {tot_m/1e6:.2f} M in total; percentages are of that total. '
                     'Samples coloured by declared species (magenta human, green mouse).\nUnmapped molecules are unmapped reads collapsed by bobcode + UMI. '
                     'Composition of the aligned molecules from the per-sample composition table (computed on the deduplicated BAMs), scaled onto the aligned total.',
                     fontsize=7.5, color=MUT)
            pdf.savefig(fig, bbox_inches='tight'); plt.close(fig)

        # ---- page 2: reads per sample --------------------------------------
        fig, ax = plt.subplots(figsize=(8.5, 9))
        v = [int(ps[s]['reads_assigned']) for s in order]
        ax.barh(y, v, color=col, height=0.72, zorder=3)
        for yy, vv in zip(y, v):
            ax.text(vv + max(v) * 0.01, yy, f'{vv:,}', va='center', fontsize=7.5, color=MUT)
        ax.axvline(n_tot / len(order), ls='--', lw=1.1, color=MUT, zorder=4)
        ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8, family='monospace')
        ax.set_xlabel('reads assigned to this sample')
        ax.set_title(f'Reads per sample\n{run}', loc='left', fontsize=12, fontweight='bold', color=INK)
        ax.grid(True, axis='x', ls=':', lw=0.5, color='#C8CFCD', zorder=0)
        for s_ in ('top', 'right'):
            ax.spines[s_].set_visible(False)
        h = [plt.Rectangle((0, 0), 1, 1, color=HUM), plt.Rectangle((0, 0), 1, 1, color=MOU)]
        ax.legend(h, ['human', 'mouse'], frameon=False, ncol=2, loc='lower right', fontsize=8.5)
        fig.text(0.01, 0.005,
                 f'{len(order)} samples, {n_tot:,} of {raw:,} raw reads assigned to a sample '
                 f'({pc(n_tot):.1f}%). Dashed = even-split mean. '
                 + ('One Illumina library: the sample is the bobcode read at the start of R2.'
                    if illumina else
                    'Both ONT primer sets pooled — the ONT barcode is the PCR primer set, not the sample.'),
                 fontsize=7.5, color=MUT)
        pdf.savefig(fig, bbox_inches='tight'); plt.close(fig)

        # ---- page 3: composition -------------------------------------------
        # Illumina: the BAM is post-dedup, so 'Unmapped' in the table is everything dedup
        # removed plus the true unmapped (~75% of prepped reads on a deep run). Plotting
        # that squashes the biology into a quarter of the axis, so the page shows the
        # composition of the KEPT reads and states the removed share in the subtitle;
        # the table keeps the all-reads denominator the validator requires.
        kept_view = bool(sc['dedup_kept_aligned'])
        def _share(s_, cat):
            v = float(comp[s_][cat])
            if not kept_view:
                return v
            kept = 100.0 - float(comp[s_]['Unmapped'])
            return 0.0 if cat == 'Unmapped' or kept <= 0 else 100.0 * v / kept
        fig, ax = plt.subplots(figsize=(8.5, 9))
        left = [0.0] * len(order)
        for cat, c in CATS:
            vals = [_share(s, cat) for s in order]
            ax.barh(y, vals, left=left, color=c, height=0.72, zorder=3)
            left = [a + b for a, b in zip(left, vals)]
        resid = [max(0.0, 100 - l) for l in left]
        if max(resid) > 0.5:                   # only real leftovers get a bar
            ax.barh(y, resid, left=left, color='#D8DCDA', height=0.72, zorder=3)
        ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8, family='monospace')
        ax.set_xlabel("% of the sample's reads kept after deduplication" if kept_view
                      else "% of the sample's assigned reads"); ax.set_xlim(0, 100)
        if kept_view:
            sub = (f'of the kept reads {100.0*n_rrna/n_map:.0f}% rRNA, {100.0*n_mrna/n_map:.0f}% mRNA; '
                   f'{100.0*n_unmap/n_tot:.0f}% of prepped reads were {unmapped_label}')
        else:
            sub = (f'overall {100.0*n_rrna/n_tot:.0f}% rRNA, {100.0*n_mrna/n_tot:.0f}% mRNA, '
                   f'{100.0*n_unmap/n_tot:.0f}% {unmapped_label}')
        ax.set_title(f'Read composition per sample\n{run}   ·   {sub}',
                     loc='left', fontsize=12, fontweight='bold', color=INK)
        shown = [(c, col_) for c, col_ in CATS if not (kept_view and c == 'Unmapped')]
        hh = [plt.Rectangle((0, 0), 1, 1, color=c) for _, c in shown]
        if max(resid) > 0.5:
            hh.append(plt.Rectangle((0, 0), 1, 1, color='#D8DCDA'))
        ax.legend(hh, [c.split(' (')[0] for c, _ in shown]
                  + (['other'] if max(resid) > 0.5 else []),
                  frameon=False, ncol=4, loc='upper center',
                  bbox_to_anchor=(0.5, -0.045), fontsize=7.5)
        for s_ in ('top', 'right'):
            ax.spines[s_].set_visible(False)
        pdf.savefig(fig, bbox_inches='tight'); plt.close(fig)

        # ---- page 3: the two accuracy views --------------------------------
        fig, axes = plt.subplots(1, 2, figsize=(10, 9), sharey=True)
        for ax, key, ttl in ((axes[0], 'pct_correct_all_reads', 'ALL reads (rRNA in)'),
                             (axes[1], 'pct_correct_mrna_only', 'mRNA only')):
            vals = [float(comp[s][key]) for s in order]
            # Stacked human/mouse rather than a single bar: the wrong-species fraction is
            # the point, and a lone "% correct" makes the reader infer it.
            for yy, s_, v in zip(y, order, vals):
                exp = ps[s_]['species']
                right, wrong = v, 100.0 - v
                hu = right if exp == 'human' else wrong
                ax.barh(yy, hu, color=HUM, height=0.72, zorder=3)
                ax.barh(yy, 100 - hu, left=hu, color=MOU, height=0.72, zorder=3)
            for yy, vv in zip(y, vals):
                ax.text(101.5, yy, f'{vv:.1f}', va='center', fontsize=7.5, color=MUT)
            ax.set_xlim(0, 112); ax.set_xlabel('% of species-called reads')
            ax.set_title(ttl, fontsize=11, fontweight='bold', color=INK)
            ax.grid(True, axis='x', ls=':', lw=0.5, color='#C8CFCD', zorder=0)
            for s_ in ('top', 'right'):
                ax.spines[s_].set_visible(False)
        axes[0].set_yticks(y); axes[0].set_yticklabels(labels, fontsize=8, family='monospace')
        fig.suptitle(f'Barcode accuracy per sample — {run}', x=0.01, ha='left',
                     fontsize=12, fontweight='bold', color=INK)
        fig.text(0.01, 0.005,
                 'Ground truth is the experiment: HEK is human, RAW is mouse. '
                 'rDNA is in GRCh38 but largely absent from GRCm39, so conserved rRNA '
                 'force-maps human and depresses the\nmouse samples on the LEFT panel. '
                 'The right panel excludes rRNA, MT and intergenic and is the view that '
                 'reflects barcode swapping.', fontsize=7.5, color=MUT)
        pdf.savefig(fig, bbox_inches='tight'); plt.close(fig)

        # ---- page 4: library-quality metrics that ARE per-sample properties ----
        # Everything here is a property of the sample's own molecules, so splitting it
        # per sample is meaningful. Deliberately ABSENT: flowcell/optical duplicates (a
        # property of the flowcell) and swap direction (cross-sample by definition).
        def _num(s_, key):
            v = comp.get(s_, {}).get(key, comp.get(s_, {}).get('ont_error_pct' if key == 'aln_error_pct' else key))
            try:
                return float(v)
            except (TypeError, ValueError):
                return 0.0
        if illumina:
            # the BAMs are post-dedup, so the duplicate rate comes from the dedup sidecar
            # (definition D, the library report's definition); 'insert' is per-read aligned length
            panels = [('insert_median', 'aligned R2 length (nt, median)', _num),
                      ('aln_error_pct', 'alignment error rate (%)', _num),
                      ('dup_pct_defD_aligned', 'duplicate rate, definition D (%)',
                       lambda s_, k: float(scs.get(s_.replace(' ', '_').replace('/', '_'), {}).get(k, 0) or 0))]
        else:
            panels = [('insert_median', "insert length (bp, median)", _num),
                      ('aln_error_pct', 'ONT error rate (%)', _num),
                      ('dup_pct_mrna', 'duplicate rate, mRNA (%)', _num)]
        fig, axes = plt.subplots(1, 3, figsize=(12, 9), sharey=True)
        for ax, (key, ttl, get) in zip(axes, panels):
            vals = [get(s, key) for s in order]
            ax.barh(y, vals, color=col, height=0.72, zorder=3)
            hi = max(vals) if vals else 1
            for yy, vv in zip(y, vals):
                ax.text(vv + hi * 0.02, yy, f'{vv:g}', va='center', fontsize=7.5, color=MUT)
            ax.set_xlim(0, hi * 1.22)
            ax.set_title(ttl, fontsize=11, fontweight='bold', color=INK)
            ax.grid(True, axis='x', ls=':', lw=0.5, color='#C8CFCD', zorder=0)
            for s_ in ('top', 'right'):
                ax.spines[s_].set_visible(False)
        axes[0].set_yticks(y); axes[0].set_yticklabels(labels, fontsize=8, family='monospace')
        fig.suptitle(f'Library quality per sample — {run}', x=0.01, ha='left',
                     fontsize=12, fontweight='bold', color=INK)
        fig.text(0.01, 0.005,
                 ('Duplicate rate is definition D from each sample\'s dedup sidecar (same 5\' '
                  'position + UMI, aligned reads), the library report\'s definition. Aligned '
                  'length is per read\nand capped by the read length; error rate = substitutions '
                  '+ indels over aligned bases (STAR nM). Flowcell duplicates and swap direction '
                  'are library-level: see the library report.')
                 if illumina else
                 ('Duplicate rate is over protein-coding mRNA. The all-reads rate is much '
                  'higher but inflated: GC-rich rDNA basecalls near-identically across\n'
                  'INDEPENDENT molecules. Flowcell duplicates and swap direction are NOT '
                  'shown here -- they are library-level and cross-sample properties.'),
                 fontsize=7.5, color=MUT)
        pdf.savefig(fig, bbox_inches='tight'); plt.close(fig)

        # ---- page 4b (paired-end Illumina only): insert length per sample -------------
        # From the dedup sidecars' `insert_size` block, measured on the FULL pair BAM before
        # dedup (one value per sequenced fragment whose two mates aligned to one contig):
        # |TLEN| minus the introns either mate spans, plus the outer soft clips (R2 beyond
        # the G-run, R1's hexamer). Absent on the poly-dT chemistry, where R1 is not aligned.
        def _side(s_):
            return scs.get(s_.replace(' ', '_').replace('/', '_'), {})
        ins_blocks = {s_: _side(s_).get('insert_size') for s_ in order}
        ins_blocks = {k: v for k, v in ins_blocks.items() if v and v.get('n')}
        if illumina and ins_blocks:
            fig, axes = plt.subplots(1, 2, figsize=(12, 9), gridspec_kw={'width_ratios': [1, 1.3]})
            ax = axes[0]
            meds = [ins_blocks.get(s_, {}).get('median') or 0 for s_ in order]
            q1s = [ins_blocks.get(s_, {}).get('q1') or 0 for s_ in order]
            q3s = [ins_blocks.get(s_, {}).get('q3') or 0 for s_ in order]
            ax.barh(y, meds, color=col, height=0.72, zorder=3)
            ax.errorbar(meds, y, xerr=[[m - a_ for m, a_ in zip(meds, q1s)], [b_ - m for m, b_ in zip(meds, q3s)]],
                        fmt='none', ecolor=INK, elinewidth=0.8, capsize=2, zorder=4)
            hi = max(q3s) if q3s else 1
            for yy_, s_ in zip(y, order):
                bl = ins_blocks.get(s_)
                if bl:
                    ax.text(bl['q3'] + hi * 0.02, yy_,
                            f"{bl['median']:g}  ({bl['pct_shorter_than_r1']:.0f}% < R1)",
                            va='center', fontsize=7.2, color=MUT)
            ax.set_xlim(0, hi * 1.45)
            ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8, family='monospace')
            ax.set_title('insert length (nt): median, IQR', fontsize=11, fontweight='bold', color=INK)
            ax.set_xlabel('cDNA insert, nt')
            ax.grid(True, axis='x', ls=':', lw=0.5, color='#C8CFCD', zorder=0)
            for s_ in ('top', 'right'):
                ax.spines[s_].set_visible(False)
            ax = axes[1]
            lib = {}
            for s_, bl in ins_blocks.items():
                h = {int(k): int(v) for k, v in (bl.get('hist') or {}).items()}
                n_ = sum(h.values()) or 1
                xs = sorted(h)
                ax.plot(xs, [h[x_] / n_ for x_ in xs], lw=0.7, alpha=0.45, color=col[order.index(s_)] if isinstance(col, list) else MUT, zorder=2)
                for k_, v_ in h.items():
                    lib[k_] = lib.get(k_, 0) + v_
            n_ = sum(lib.values()) or 1
            xs = sorted(lib)
            ax.plot(xs, [lib[x_] / n_ for x_ in xs], lw=1.8, color=INK, zorder=3, label='library (all samples)')
            r1len = next((bl.get('r1_read_len_mode') for bl in ins_blocks.values() if bl.get('r1_read_len_mode')), None)
            if r1len:
                ax.axvline(r1len, ls=':', lw=1, color=MUT)
                ax.text(r1len, ax.get_ylim()[1] * 0.98, f' R1 read length {r1len}: shorter inserts read into the adapter',
                        fontsize=7.5, color=MUT, va='top')
            cap = next((bl.get('cap') for bl in ins_blocks.values() if bl.get('cap')), None)
            # x-range: the 99th percentile of the library histogram (the cap bin holds the
            # non-overlapping, intron-spanning fragments and would stretch the axis to the cap)
            body = [x_ for x_ in xs if not cap or x_ < cap]           # the cap bin is not a length
            n_body = sum(lib[x_] for x_ in body) or 1
            acc = 0; x995 = max(body) if body else 1500
            for x_ in body:
                acc += lib[x_]
                if acc >= 0.99 * n_body:
                    x995 = x_; break
            ax.set_xlim(0, max((x995 // 100) * 100, 300, (r1len or 0) * 2))   # 99th percentile, rounded down to 100 nt (600 on the 24-plex)
            if cap and cap in lib:
                ax.text(0.98, 0.55, f'{100 * lib[cap] / n_:.1f}% of fragments at the {cap} nt cap\n(mates not overlapping; off scale)',
                        transform=ax.transAxes, ha='right', va='center', fontsize=7.5, color=MUT)
            ax.set_xlabel('cDNA insert, nt'); ax.set_ylabel('fraction of fragments')
            ax.set_title('insert length distribution (thin: per sample; bold: library)', fontsize=11, fontweight='bold', color=INK)
            ax.legend(frameon=False, fontsize=8, loc='center right')
            ax.grid(True, ls=':', lw=0.5, color='#C8CFCD', zorder=0)
            for s_ in ('top', 'right'):
                ax.spines[s_].set_visible(False)
            n_tot = sum(bl.get('n', 0) for bl in ins_blocks.values())
            ov = sum(bl.get('n_mates_overlap', 0) for bl in ins_blocks.values())
            fig.suptitle(f'Insert length per sample (paired-end) \u2014 {run}', x=0.01, ha='left',
                         fontsize=12, fontweight='bold', color=INK)
            fig.text(0.01, 0.005,
                     f'Per sequenced fragment, before dedup, both mates aligned to one contig: {n_tot:,} fragments; '
                     f'{100 * ov / max(n_tot, 1):.1f}% with overlapping mates. Insert = |TLEN| minus the introns either mate spans, '
                     'plus the outer soft clips\n(R2 beyond the non-templated G-run, R1 including the priming hexamer). Fragments '
                     f'whose mates do not overlap can hide an intron between them and are counted as read. Bins capped at {cap} nt.',
                     fontsize=7.5, color=MUT)
            pdf.savefig(fig, bbox_inches='tight'); plt.close(fig)

        # ---- pages 5-7: the structural metrics, from the sample-scope analysis ----
        # Row-oriented on purpose. The library report puts each unit in a COLUMN, which
        # is right for 2-4 ONT barcodes and unreadable at 24 samples; the numbers are the
        # same, only the axis changes.
        PANEL_SETS = [
            ('Read structure', [
                ('% perfect TSO backbone (full 18bp)', 'perfect TSO backbone (%)'),
                ('% reads with polyT', 'reads with polyT (%)'),
                ('% reads with \u22652 TSO (chimera)', 'TSO-TSO chimeras (%)')]),
            ('Template switching and artifacts', [
                ('G-run mean / median (non-templated G after TSO)', 'G-run mean (bp)'),
                ('No template switch (G-run = 0)', 'no template switch (%)'),
                ('Empty products % (of all reads)', 'empty products (%)')]),
            ('Duplicate funnel', [
                ('Dup funnel \u2014 same-gene dup %', 'same-gene dup (%)'),
                ('Dup funnel \u2014 + UMI dup %', '+ UMI dup (%)'),
                ('Dup funnel \u2014 + same strand = PCR dup %', 'PCR dup (%)')]),
        ]
        # A panel set is suppressed when its own denominator is too small to support a
        # percentage. On a shallow run only 2-5 reads per sample may carry a callable UMI,
        # so the duplicate funnel would plot percentages computed on n=3. Showing that is worse
        # than showing nothing, so it is skipped and named in the caption instead.
        GUARDS = {'Duplicate funnel':
                  ('Dup funnel \u2014 mRNA reads (non-ribo) with callable UMI', 50)}
        suppressed = []
        # Illumina: these are ONT read-structure metrics (TSO backbone, polyT, the ONT
        # duplicate funnel) and the Illumina module passes no per-sample master summary,
        # so the pages are not drawn and nothing is reported as suppressed.
        for title, panels in ([] if illumina else PANEL_SETS):
            guard = GUARDS.get(title)
            if guard:
                metric, floor = guard
                ns = [msnum(s, metric) for s in order]
                ns = [n for n in ns if n is not None]
                med = statistics.median(ns) if ns else 0
                if med < floor:
                    suppressed.append(f'{title} (median n={med:g} per sample, needs \u2265{floor})')
                    continue
            vals_by_panel = [[msnum(s, m) for s in order] for m, _ in panels]
            if all(all(v is None for v in vv) for vv in vals_by_panel):
                continue                                   # metric absent for this chemistry
            fig, axes = plt.subplots(1, len(panels), figsize=(4 * len(panels), 9), sharey=True)
            axes = axes if len(panels) > 1 else [axes]
            for ax, (metric, lab), vals in zip(axes, panels, vals_by_panel):
                v = [0.0 if x is None else x for x in vals]
                ax.barh(y, v, color=col, height=0.72, zorder=3)
                hi = max(v) if any(v) else 1.0
                for yy, vv, raw in zip(y, v, vals):
                    ax.text(vv + hi * 0.02, yy, '-' if raw is None else f'{vv:g}',
                            va='center', fontsize=7.5, color=MUT)
                ax.set_xlim(0, hi * 1.25)
                ax.set_title(lab, fontsize=10.5, fontweight='bold', color=INK)
                ax.grid(True, axis='x', ls=':', lw=0.5, color='#C8CFCD', zorder=0)
                for s_ in ('top', 'right'):
                    ax.spines[s_].set_visible(False)
            axes[0].set_yticks(y)
            axes[0].set_yticklabels(labels, fontsize=8, family='monospace')
            fig.suptitle(f'{title} \u2014 {run}', x=0.01, ha='left', fontsize=12,
                         fontweight='bold', color=INK)
            fig.text(0.01, 0.005, 'From the sample-scope analysis (the same code that '
                                  'produces the library report), one row per sample.',
                     fontsize=7.5, color=MUT)
            pdf.savefig(fig, bbox_inches='tight'); plt.close(fig)

        # ---- last page: replicate consistency, the point of the design ----------
        conds = collections.OrderedDict()
        for s in order:
            conds.setdefault(ps[s]['condition'], []).append(s)
        fig, axes = plt.subplots(1, 2, figsize=(10, 5.2))
        names = list(conds)
        means = [statistics.mean(int(ps[x]['reads_assigned']) for x in conds[c]) for c in names]
        cvs = [(statistics.pstdev([int(ps[x]['reads_assigned']) for x in conds[c]]) /
                statistics.mean([int(ps[x]['reads_assigned']) for x in conds[c]]) * 100)
               if len(conds[c]) > 1 else 0.0 for c in names]
        yy = list(range(len(names)))[::-1]
        axes[0].barh(yy, means, color=MUT, alpha=0.75, height=0.66, zorder=3)
        axes[0].set_title('mean reads per condition', fontsize=10.5, fontweight='bold', color=INK)
        axes[1].barh(yy, cvs, color=[('#A5342B' if c > 20 else '#0B6E63') for c in cvs],
                     height=0.66, zorder=3)
        axes[1].set_title('replicate CV (%)  — red above 20%', fontsize=10.5,
                          fontweight='bold', color=INK)
        for ax in axes:
            ax.set_yticks(yy); ax.grid(True, axis='x', ls=':', lw=0.5, color='#C8CFCD', zorder=0)
            for s_ in ('top', 'right'):
                ax.spines[s_].set_visible(False)
        axes[0].set_yticklabels(names, fontsize=9)
        axes[1].set_yticklabels([])
        for a, v in zip(yy, cvs):
            axes[1].text(v + 0.5, a, f'{v:.0f}%', va='center', fontsize=8, color=MUT)
        fig.suptitle(f'Replicate consistency \u2014 {run}', x=0.01, ha='left',
                     fontsize=12, fontweight='bold', color=INK)
        _note = ('CV between the replicates of each condition. This is what a replicated '
                 'design is for, and nothing else in the run reports it.')
        if suppressed:
            _note += ('\nSuppressed for insufficient depth, not omitted: '
                      + '; '.join(suppressed) + '.')
        fig.text(0.01, 0.005, _note, fontsize=7.5, color=MUT)
        pdf.savefig(fig, bbox_inches='tight'); plt.close(fig)

    sys.stderr.write(f'per_sample_report: {len(order)} samples -> {args.out}\n')


if __name__ == '__main__':
    main()

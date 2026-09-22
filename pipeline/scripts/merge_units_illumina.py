#!/usr/bin/env python3
"""Merge per-sample analysis results into ONE library-level unit (per-sample mode).

The Illumina pipeline processes every sample of a pooled library on its own (prep chunks ->
demultiplex by the bobcode prep wrote into the read name -> STAR, dedup and the analyzer per
sample). The library-level report is then rendered from a unit assembled here, by the same
report code and in the same single-column format as a pooled run, so the numbers are the
pooled numbers: every quantity the report shows is a count, a histogram, or a ratio of counts,
and counts add. Where a per-unit field is a derived value (a percentage, a median, a top-N
list) it is recomputed here from the merged counts with the same formula the per-unit code
uses; the per-unit modules store the counts those formulas need
(fullmol_funnel *_correct_n, star_gene_reads, rt_extent dist_counts).

Fields that are NOT exactly mergeable are the legacy ones the Illumina report does not print
(the position-only duplicate groups, the UMI-only dup rate, the insert-length medians of the
ONT structure detectors); they are summed or averaged and listed under `merged_from.approx`.

    merge_units_illumina.py --out <library>_speciesmix_results.json --unit <library name>
        [--prep-stats merged_prep_stats.json] <sample1_speciesmix_results.json> <sample2 ...>
"""
import argparse, json, os, sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def _is_int(v):
    return isinstance(v, int) and not isinstance(v, bool)


def _sum_counts(vals):
    """Recursive sum of count structures: ints add, dicts merge by key, equal-length lists of
    ints add elementwise. Non-numeric leaves must agree (first wins, disagreement noted)."""
    vals = [v for v in vals if v is not None]
    if not vals:
        return None
    v0 = vals[0]
    if all(_is_int(v) for v in vals):
        return sum(vals)
    if all(isinstance(v, float) for v in vals):
        return sum(vals)
    if all(isinstance(v, dict) for v in vals):
        keys = []
        for v in vals:
            for k in v:
                if k not in keys:
                    keys.append(k)
        return {k: _sum_counts([v[k] for v in vals if k in v]) for k in keys}
    if all(isinstance(v, list) for v in vals):
        if all(len(v) == len(v0) for v in vals) and all(all(_is_int(x) or isinstance(x, float) for x in v) for v in vals):
            return [sum(v[i] for v in vals) for i in range(len(v0))]
        out = []
        for v in vals:
            out.extend(x for x in v if x not in out)
        return out
    return v0


def _pct(a, b, nd=2):
    return round(100.0 * a / b, nd) if b else None


def _hist_stats(counter):
    """(q25, median, q75, pct<150, pct<500, pct>2000) from {distance: count}, with the
    per-unit rule: quantile p = the value at index floor(p*n) of the sorted distances."""
    items = sorted((int(k), int(v)) for k, v in counter.items())
    n = sum(v for _, v in items)
    if not n:
        return {'q25': None, 'median': None, 'q75': None, 'pct_lt150': None, 'pct_lt500': None, 'pct_gt2000': None}
    def q(p):
        idx = min(n - 1, int(p * n)); c = 0
        for k, v in items:
            c += v
            if c > idx:
                return k
        return items[-1][0]
    return {'q25': q(0.25), 'median': q(0.5), 'q75': q(0.75),
            'pct_lt150': round(100 * sum(v for k, v in items if k < 150) / n, 1),
            'pct_lt500': round(100 * sum(v for k, v in items if k < 500) / n, 1),
            'pct_gt2000': round(100 * sum(v for k, v in items if k > 2000) / n, 1)}


def merge_fullmol(parts):
    parts = [p for p in parts if p]
    if not parts:
        return None
    n = sum(p['n_total'] for p in parts)
    steps = []
    for i, st0 in enumerate(parts[0]['steps']):
        ss = [p['steps'][i] for p in parts]
        if any(s['key'] != st0['key'] for s in ss):
            sys.exit('merge_units: funnel steps differ between samples')
        cum = sum(s['n'] for s in ss); sc = sum(s['acc_n'] for s in ss); co = sum(s['acc_correct'] for s in ss)
        steps.append({'key': st0['key'], 'label': st0['label'], 'n': cum, 'pct': round(100 * cum / n, 1),
                      'acc': round(100 * co / sc, 2) if sc else None, 'acc_n': sc, 'acc_correct': co})
    out = {'n_total': n, 'steps': steps,
           'correct_vs_len': _sum_counts([p.get('correct_vs_len') for p in parts]),
           'n_ambiguous_excluded': sum(p.get('n_ambiguous_excluded', 0) for p in parts)}
    for pfx in ('ms', 'hu'):
        sn = sum(p.get(f'{pfx}_n', 0) for p in parts)
        if all(f'{pfx}_correct_n' in p for p in parts):
            sc = sum(p[f'{pfx}_correct_n'] for p in parts)
        else:                                   # older per-unit JSON: back out the count from the %
            sc = sum(round((p.get(f'{pfx}_correct') or 0) * p.get(f'{pfx}_n', 0) / 100) for p in parts)
        out[f'{pfx}_correct'] = round(100 * sc / sn, 2) if sn else None
        out[f'{pfx}_n'] = sn
        out[f'{pfx}_correct_n'] = sc
    return out


def merge_funnel_spec(parts):
    parts = [p for p in parts if p]
    if not parts:
        return None
    out = {}
    for k in parts[0]:
        if k.startswith('_'):
            continue
        c = sum(p[k]['correct'] for p in parts if k in p); n = sum(p[k]['n'] for p in parts if k in p)
        nr = sum(p[k].get('n_reach', 0) for p in parts if k in p)
        out[k] = {'correct': c, 'n': n, 'n_reach': nr, 'label': parts[0][k]['label'], 'specificity': _pct(c, n)}
    for pfx in ('ms', 'hu'):
        n = sum(p.get(f'_{pfx}_n', 0) for p in parts)
        if all(f'_{pfx}_correct_n' in p for p in parts):
            c = sum(p[f'_{pfx}_correct_n'] for p in parts)
        else:
            c = sum(round((p.get(f'_{pfx}_correct') or 0) * p.get(f'_{pfx}_n', 0) / 100) for p in parts)
        out[f'_{pfx}_correct'] = _pct(c, n); out[f'_{pfx}_n'] = n; out[f'_{pfx}_correct_n'] = c
    return out


def merge_filter_funnel(parts):
    parts = [p for p in parts if p]
    if not parts:
        return None
    out = {}
    steps = []
    for i, s0 in enumerate(parts[0].get('steps', [])):
        ss = [p['steps'][i] for p in parts]
        n = sum(s['n'] for s in ss); nb = sum(s['n_barcoded'] for s in ss); sw = sum(s['n_swap'] for s in ss)
        steps.append({'n': n, 'n_barcoded': nb, 'n_swap': sw,
                      'acc': (100.0 * (nb - sw) / nb) if nb else None, 'ci_lo': None, 'ci_hi': None,
                      'label': s0['label']})
    out['steps'] = steps
    gh = _sum_counts([p.get('grun_hist') for p in parts]) or {}
    tot = sum(gh.values())
    out['grun_hist'] = {str(k): v for k, v in sorted(((int(k), v) for k, v in gh.items()))}
    out['grun_n'] = tot
    out['grun_mean'] = (sum(int(k) * v for k, v in gh.items()) / tot) if tot else None
    med = None
    if tot:
        c = 0
        for k, v in sorted((int(k), v) for k, v in gh.items()):
            c += v
            if c >= (tot + 1) // 2:
                med = k; break
    out['grun_median'] = med
    for k in ('ms_correct', 'hu_correct', 'terminal_correct'):
        vals = [(p.get(k), p.get('n_final', 0)) for p in parts if p.get(k) is not None]
        w = sum(n for _, n in vals)
        out[k] = (sum(v * n for v, n in vals) / w) if w else parts[0].get(k)
    out['n_final'] = sum(p.get('n_final', 0) for p in parts)
    return out


def merge_gene_stats(units):
    """(star_transcript_diversity, star_top_genes, star_gene_reads) from the per-unit
    star_gene_reads tables, with _star_gene_stats' own formulas."""
    tables = [u.get('star_gene_reads') for u in units]
    if not all(tables):
        return None, None, None
    by_sp = {}
    for t in tables:
        for sp, genes in t.items():
            d = by_sp.setdefault(sp, {})
            for g, (n, m) in genes.items():
                e = d.setdefault(g, [0, 0]); e[0] += n; e[1] += m
    mrna_total = Counter()
    for sp, genes in by_sp.items():
        for g, (n, _m) in genes.items():
            mrna_total[g] += n
    total = sum(mrna_total.values())
    top = {sp: {'total': sum(n for n, _ in genes.values()),
                'top': [[g, n, round(m / n, 2)] for g, (n, m) in sorted(genes.items(), key=lambda kv: (-kv[1][0], kv[0]))[:7]]}
           for sp, genes in by_sp.items() if sp in ('human', 'mouse')}
    div = None
    if total:
        counts = sorted(mrna_total.values(), reverse=True)
        simpson = 1.0 - sum((c / total) ** 2 for c in counts)
        topk = lambda k: round(sum(counts[:k]) / total * 100, 1)
        div = {'total_mapped': total, 'unique_genes': len(mrna_total),
               'unique_genes_human': len(by_sp.get('human', {})), 'unique_genes_mouse': len(by_sp.get('mouse', {})),
               'rrna_masked': sum((u.get('star_transcript_diversity') or {}).get('rrna_masked', 0) for u in units),
               'simpson': round(simpson, 4), 'top1_pct': topk(1), 'top5_pct': topk(5), 'top10_pct': topk(10)}
    return div, top, by_sp


def merge_rt_extent(parts):
    parts = [p for p in parts if p and p.get('n')]
    if not parts:
        return None
    if not all('dist_counts' in p for p in parts):
        sys.exit('merge_units: rt_extent without dist_counts (per-unit code predates the merge)')
    dc = Counter()
    for p in parts:
        dc.update({int(k): int(v) for k, v in p['dist_counts'].items()})
    out = {k: sum(p[k] for p in parts) for k in ('n', 'n_primary', 'n_unique', 'n_in_gene')}
    for k in ('binw', 'maxd', 'fbinw', 'fmax', 'anchor', 'scope'):
        out[k] = parts[0][k]
    out['hist'] = [sum(p['hist'][i] for p in parts) for i in range(len(parts[0]['hist']))]
    out['fhist'] = [sum(p['fhist'][i] for p in parts) for i in range(len(parts[0]['fhist']))]
    out['dist_counts'] = {str(k): v for k, v in sorted(dc.items())}
    out.update(_hist_stats(dc))
    return out


def merge_sidecar(parts):
    import dedup_reads_illumina as dr
    return dr.merge_sidecars(parts)


def merge_swap(parts):
    parts = [p for p in parts if p]
    if not parts:
        return None
    out = {'per_category': {}}
    for p in parts:
        for cat, v in p['per_category'].items():
            e = out['per_category'].setdefault(cat, {'correct': 0, 'n': 0})
            e['correct'] += v['correct']; e['n'] += v['n']
    for e in out['per_category'].values():
        e['specificity'] = _pct(e['correct'], e['n'])
    c = sum(p['overall_rrna_excluded']['correct'] for p in parts); n = sum(p['overall_rrna_excluded']['n'] for p in parts)
    out['overall_rrna_excluded'] = {'correct': c, 'n': n, 'specificity': _pct(c, n)}
    for k in parts[0]:
        if k not in out:
            out[k] = _sum_counts([p.get(k) for p in parts])
    return out


def merge_detected(parts):
    parts = [p for p in parts if p]
    if not parts:
        return None
    out = {'architecture': parts[0].get('architecture'),
           'n_scanned': sum(p.get('n_scanned', 0) for p in parts),
           'arch_counts': _sum_counts([p.get('arch_counts') for p in parts]),
           'matched_by_arch': _sum_counts([p.get('matched_by_arch') for p in parts])}
    by_seq = {}
    for p in parts:
        for b in p.get('barcodes') or []:
            e = by_seq.setdefault(b['seq'], dict(b)); 
            if e is not b:
                e['count'] = e.get('count', 0) + b.get('count', 0)
    for e in by_seq.values():
        e['frac'] = (e.get('count', 0) / out['n_scanned']) if out['n_scanned'] else 0.0
    out['barcodes'] = sorted(by_seq.values(), key=lambda b: -b.get('count', 0))
    out['n_matched'] = sum(p.get('n_matched', 0) for p in parts)
    out['frac_matched'] = (out['n_matched'] / out['n_scanned']) if out['n_scanned'] else 0.0
    out['anchor_signature'] = parts[0].get('anchor_signature')
    return out


# per-unit keys whose merged value is simply the recursive count-sum (all report-facing)
_SUM_KEYS = ('star_taxonomy', 'star_crosstab', 'star_crosstab_mtrep', 'star_crosstab_clean',
             'star_grun_accuracy', 'tso7_retention', 'tso7_retention_no_rrna', 'structure_counts',
             'bio_groups', 'seq_ref_counts', 'tso_barcode_counts', 'bob_barcode_counts',
             'read_tso_barcodes', 'read_bob_barcodes', 'star_len_hist', 'flank_freq',
             'star_pca_gene_counts', 'tso_multiplicity', 'illumina_qc', 'read_composition',
             'total_reads', 'insert_count')
# scalar counts at the top level (reads_with_* etc.): summed when every unit has an int
_CONST_KEYS = ('construct_type', 'star_taxonomy_order', 'star_taxonomy_scope', 'star_len_hist_labels',
               'rt_inferred', 'fastq_file')
# legacy derived values without exact merge: summed/averaged and flagged
_APPROX_KEYS = ('duplicates', 'star_duplicates', 'len_mean', 'len_median', 'len_min', 'len_max',
                'bio_group_insert_stats', 'structure_insert_stats', 'empty_products', 'polydt_less')


def merge_units(units, name, prep_stats=None):
    out = {'barcode_name': name}
    approx = []
    for k in _SUM_KEYS:
        vals = [u.get(k) for u in units if k in u]
        if vals:
            out[k] = _sum_counts(vals)
    for k in _CONST_KEYS:
        vals = [u.get(k) for u in units if k in u]
        if vals:
            out[k] = vals[0]
    # every other top-level int (reads_with_*, reads_multi_*, reads_complete, ...) adds
    for k in units[0]:
        if k in out or k in _APPROX_KEYS:
            continue
        vals = [u.get(k) for u in units]
        if all(_is_int(v) for v in vals):
            out[k] = sum(vals)
    for k in _APPROX_KEYS:
        vals = [u.get(k) for u in units if u.get(k) is not None]
        if not vals:
            continue
        if k in ('len_min',):
            out[k] = min(vals)
        elif k in ('len_max',):
            out[k] = max(vals)
        elif k in ('len_mean', 'len_median'):
            w = [u.get('insert_count', 0) or 0 for u in units if u.get(k) is not None]
            out[k] = (sum(v * n for v, n in zip(vals, w)) / sum(w)) if sum(w) else vals[0]
        else:
            out[k] = _sum_counts(vals)
        approx.append(k)
    out['fullmol_funnel'] = merge_fullmol([u.get('fullmol_funnel') for u in units])
    out['star_funnel_spec'] = merge_funnel_spec([u.get('star_funnel_spec') for u in units])
    out['filter_funnel'] = merge_filter_funnel([u.get('filter_funnel') for u in units])
    div, top, gene_reads = merge_gene_stats(units)
    if div is not None:
        out['star_transcript_diversity'] = div; out['star_top_genes'] = top; out['star_gene_reads'] = gene_reads
    out['rt_extent'] = merge_rt_extent([u.get('rt_extent') for u in units])
    out['dup_funnel'] = merge_sidecar([u.get('dup_funnel') for u in units])
    out['star_swap'] = merge_swap([u.get('star_swap') for u in units])
    out['detected_bobcodes'] = merge_detected([u.get('detected_bobcodes') for u in units])
    out['star_summary'] = merge_star_summary(units, out)
    if prep_stats is not None:
        out['prep_stats'] = prep_stats
    out['merged_from'] = {'units': [u.get('barcode_name') for u in units], 'approx': approx}
    return out


def merge_star_summary(units, merged):
    """star_summary is a bag of derived values; recompute the ones with a count formula from
    the merged structures, weight the rest by mapped reads and flag them."""
    ss = [u.get('star_summary') or {} for u in units]
    if not any(ss):
        return None
    tax = merged.get('star_taxonomy') or {}
    mapped = sum(sum(v.values()) for v in tax.values() if isinstance(v, dict))
    out = {}
    def cat_n(cat):
        return sum((tax.get(cat) or {}).values())
    for k in ss[0]:
        vals = [s.get(k) for s in ss if k in s]
        if all(_is_int(v) for v in vals):
            out[k] = sum(vals)
        elif all(isinstance(v, dict) for v in vals):
            out[k] = _sum_counts(vals)
        elif all(isinstance(v, str) for v in vals):
            out[k] = vals[0]
        else:
            w = [sum(sum(v.values()) for v in (u.get('star_taxonomy') or {}).values() if isinstance(v, dict)) for u in units]
            pairs = [(v, n) for v, n in zip(vals, w) if v is not None]
            out[k] = (sum(v * n for v, n in pairs) / sum(n for _, n in pairs)) if pairs and sum(n for _, n in pairs) else (vals[0] if vals else None)
    # exact recomputations
    if mapped:
        out['mrna_pct'] = round(100 * cat_n('mRNA (protein-coding)') / mapped, 2)
        out['rrna_pct'] = round(100 * cat_n('rRNA') / mapped, 2)
        out['mt_pct'] = round(100 * cat_n('Mitochondrial') / mapped, 2)
    mr = tax.get('mRNA (protein-coding)') or {}
    if mr:
        n_m = sum(mr.values())
        out['mrna_hu_n'] = mr.get('human', 0); out['mrna_ms_n'] = mr.get('mouse', 0); out['mrna_amb_n'] = mr.get('ambiguous', 0)
        out['hu_mrna_pct'] = 100.0 * mr.get('human', 0) / n_m if n_m else None
        out['ms_mrna_pct'] = 100.0 * mr.get('mouse', 0) / n_m if n_m else None
    sw = merged.get('star_swap') or {}
    if sw.get('overall_rrna_excluded'):
        out['pct_correct'] = sw['overall_rrna_excluded']['specificity']
    if all('terminal_n' in x and 'aln_n' in x for x in ss):
        tn = sum(x['terminal_n'] for x in ss); an = sum(x['aln_n'] for x in ss)
        out['terminal_n'] = tn; out['aln_n'] = an
        out['terminal_pct'] = round(tn / an * 100, 1) if an else None
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('units', nargs='+', help='per-sample *_speciesmix_results.json')
    ap.add_argument('--out', required=True); ap.add_argument('--unit', required=True, help='library unit name, e.g. barcode01')
    ap.add_argument('--prep-stats', default=None, help='merged prep_stats.json of the library')
    a = ap.parse_args()
    units = [json.load(open(p)) for p in a.units]
    prep = json.load(open(a.prep_stats)) if a.prep_stats else None
    out = merge_units(units, a.unit, prep)
    json.dump(out, open(a.out, 'w'), indent=1)
    print(f'merge_units: {len(units)} unit(s) -> {a.out}; approx fields: {out["merged_from"]["approx"]}')


if __name__ == '__main__':
    main()

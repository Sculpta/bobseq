#!/usr/bin/env python3
"""Log.final.selected.out for an arm whose STAR run covered more samples than the benchmark uses (BOBseq: 18 aligned, 17 kept).
Input reads = sum of per-sample input counts of the kept wells (matched arm: prep_stats.json 'well:<bobcode>' written counts;
native arm: the pipeline's per-sample STAR logs); mapped / uniquely mapped = primary mapped rows / NH=1 rows of the arm's
read table (kept wells only). '% of reads unmapped: too short' is copied from the full STAR log (a library-level rate).
    make_selected_log.py <arm_dir> <input_counts.tsv: well \t input_reads> <perwell_mapped.tsv: well \t mapped \t unique>
"""

import sys, os

D, inp, pw = sys.argv[1:4]
n_wells = int(next(l.split('\t')[1] for l in open(f'{D}/stats_E_U/basic/summary.tsv') if l.startswith('n_wells\t')))
rows = [l.rstrip('\n').split('\t') for l in open(pw) if '\t' in l and not l.startswith('DONE')]
keep = {r[0] for r in rows}  # the wells of the read table = the benchmark samples
assert len(keep) == n_wells, f'{len(keep)} wells in the read table vs n_wells {n_wells} in summary.tsv'
n_in = sum(int(l.split('\t')[1]) for l in open(inp) if l.split('\t')[0] in keep)
assert sum(1 for l in open(inp) if l.split('\t')[0] in keep) == n_wells, 'input counts missing for some wells'
mapped = sum(int(r[1]) for r in rows)
uniq = sum(int(r[2]) for r in rows)
full = {l.split('|')[0].strip(): l.split('|')[1].strip() for l in open(f'{D}/Log.final.out') if '|' in l}
ts = full.get('% of reads unmapped: too short', '0%')
assert n_in >= mapped >= uniq > 0, (n_in, mapped, uniq)
with open(f'{D}/Log.final.selected.out', 'w') as f:
    f.write(
        f"# synthesized for the {len(keep)} benchmark samples (see make_selected_log.py); the arm's Log.final.out is the full STAR run\n"
    )
    f.write(f"                          Number of input reads |\t{n_in}\n")
    f.write(f"                   Uniquely mapped reads number |\t{uniq}\n")
    f.write(f"                        Uniquely mapped reads % |\t{100*uniq/n_in:.2f}%\n")
    f.write(f"        Number of reads mapped to multiple loci |\t{mapped-uniq}\n")
    f.write(f"             % of reads mapped to multiple loci |\t{100*(mapped-uniq)/n_in:.2f}%\n")
    f.write(f"                 % of reads unmapped: too short |\t{ts}\n")
print(
    f'{os.path.basename(D)}: {len(keep)} wells, input {n_in:,}, mapped {mapped:,} ({100*mapped/n_in:.1f}%), unique {uniq:,} ({100*uniq/n_in:.1f}%), too short {ts} (library rate)'
)

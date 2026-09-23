"""Stage psi: exon-7 PSI variants + c.840 T-fraction per sample and pooled per condition, with
Clopper-Pearson 95 % CIs. Fragment-level: incl = fragments with an ex6>ex7/ex7>ex8 junction, skip = ex6>ex8.

A skip fragment (ex6>ex8) carries no exon-7 sequence and is not assigned to a paralog: the skip pool is
"SMN exon-7 skipping" as a whole, and each paralog's PSI is its own inclusion against that shared pool.

  psi_pooled     (SMN1+SMN2 together)  incl_all  / (incl_all  + skip_all)       must rise under risdiplam
  psi_SMN2       SMN2 inclusion        incl_SMN2 / (incl_SMN2 + skip_all)
  psi_SMN1       SMN1 inclusion        incl_SMN1 / (incl_SMN1 + skip_all)       the same skips in the denominator
  t_frac         c.840 T / (T + C) among exon-7-covering fragments               skip-independent readout
Inclusion fragments whose exon-7 portion reaches no diagnostic base (incl_ambiguous) enter psi_pooled only.
Writes results/SMN2_psi_samples.tsv and results/SMN2_psi_conditions.tsv."""
import csv, os, sys
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

KEYS = ["incl_SMN2", "incl_SMN1", "incl_ambiguous", "skip_SMN2", "skip_SMN1", "skip_ambiguous", "c840_T", "c840_C"]
JKEYS = [f"{j}_{p}" for j in ("ex6>7", "ex7>8", "ex6>8") for p in ("SMN2", "SMN1", "ambiguous")]   # fragments per junction per paralog

def metrics(c):
    incl_all = c["incl_SMN2"] + c["incl_SMN1"] + c["incl_ambiguous"]; skip_all = c["skip_SMN2"] + c["skip_SMN1"] + c["skip_ambiguous"]
    out = {}
    def add(name, k, n):
        lo, hi = cp_ci(k, n); out[name] = round(k / n, 4) if n else ""; out[name + "_lo"] = round(lo, 4) if n else ""; out[name + "_hi"] = round(hi, 4) if n else ""; out[name + "_n"] = n
    add("psi_pooled", incl_all, incl_all + skip_all)
    add("psi_SMN2", c["incl_SMN2"], c["incl_SMN2"] + skip_all)
    add("psi_SMN1", c["incl_SMN1"], c["incl_SMN1"] + skip_all)
    add("t_frac", c["c840_T"], c["c840_T"] + c["c840_C"])
    return out

def main():
    rows = list(csv.DictReader(open(os.path.join(ANALYSIS, "results", "SMN2_samples.tsv")), delimiter="\t"))
    per_sample = []; pooled = defaultdict(lambda: defaultdict(int)); meta = {}
    for r in rows:
        c = {k: int(r[k]) for k in KEYS + JKEYS}
        per_sample.append({**{k: r[k] for k in ("acc_number", "set", "method", "sample_name", "condition", "in_benchmark", "unit")}, **c, **metrics(c)})
        grp = (r["set"], r["condition"] if r["method"] == "bobseq" else "all", r["unit"])
        if True:                                     # ALL wells of a condition pooled, Hek_control_2 and the Ris 500 mM wells included
            for k in KEYS + JKEYS: pooled[grp][k] += c[k]
            m = meta.setdefault(grp, {"set": grp[0], "condition": grp[1], "unit": grp[2], "method": r["method"], "n_samples": 0})
            m["n_samples"] += 1
    per_cond = [{**meta[g], **dict(pooled[g]), **metrics(pooled[g])} for g in sorted(pooled)]
    write_tsv(os.path.join(ANALYSIS, "results", "SMN2_psi_samples.tsv"), per_sample, list(per_sample[0].keys()))
    write_tsv(os.path.join(ANALYSIS, "results", "SMN2_psi_conditions.tsv"), per_cond, list(per_cond[0].keys()))
    # compact junction-count table: fragments per junction per paralog (skips ex6>8 have no exon-7 base -> mostly unassigned)
    jc = []
    for r in per_cond + per_sample:
        if r.get("set") not in ("bobseq_pe_native", "bobseq_50nt"): continue
        jc.append({"level": "condition" if "n_samples" in r else "well", "set": r["set"], "condition": r["condition"],
                   "well": r.get("sample_name", f"pooled ({r.get('n_samples')} wells)"), "in_benchmark": r.get("in_benchmark", ""), "unit": r["unit"],
                   **{f"{p}_{j}": r[f"{j}_{p}"] for p in ("SMN2", "SMN1", "ambiguous") for j in ("ex6>7", "ex7>8", "ex6>8")},
                   "skip_all_ex6>8": r["ex6>8_SMN2"] + r["ex6>8_SMN1"] + r["ex6>8_ambiguous"],
                   "psi_SMN2": r["psi_SMN2"], "psi_SMN1": r["psi_SMN1"], "psi_pooled": r["psi_pooled"]})
    jc.sort(key=lambda r: (r["set"], r["unit"], r["level"] != "condition", r["condition"], r["well"]))
    write_tsv(os.path.join(ANALYSIS, "results", "SMN2_junction_counts.tsv"), jc, list(jc[0].keys()))
    # significance vs Hek control (BOBseq PE): Fisher's exact test on the pooled junction reads (inclusion vs shared skips),
    # plus a Welch t-test on the per-well PSIs (n = 3 vs 3) for the replicate-level view
    from scipy.stats import fisher_exact, ttest_ind
    tests = []
    for unit in ("raw", "dedup"):
        ctl = next(r for r in per_cond if r["set"] == "bobseq_pe_native" and r["condition"] == "Hek control" and r["unit"] == unit)
        for r in per_cond:
            if r["set"] != "bobseq_pe_native" or r["unit"] != unit or r["condition"] == "Hek control": continue
            for para in ("SMN2", "SMN1"):
                sk = lambda x: x["skip_SMN2"] + x["skip_SMN1"] + x["skip_ambiguous"]
                tab = [[r[f"incl_{para}"], sk(r)], [ctl[f"incl_{para}"], sk(ctl)]]
                odds, p_f = fisher_exact(tab)
                pw = lambda c: [float(x[f"psi_{para}"]) for x in per_sample if x["set"] == "bobseq_pe_native" and x["unit"] == unit and x["condition"] == c and x[f"psi_{para}"] != ""]
                a, b = pw(r["condition"]), pw("Hek control")
                p_t = ttest_ind(a, b, equal_var=False).pvalue if len(a) > 1 and len(b) > 1 else float("nan")
                tests.append({"unit": unit, "paralog": para, "contrast": f"{r['condition']} vs Hek control", "incl_treated": tab[0][0], "skip_treated": tab[0][1],
                              "incl_control": tab[1][0], "skip_control": tab[1][1], "psi_treated": r[f"psi_{para}"], "psi_control": ctl[f"psi_{para}"],
                              "dpsi": round(float(r[f"psi_{para}"]) - float(ctl[f"psi_{para}"]), 4) if r[f"psi_{para}"] != "" else "",
                              "odds_ratio": round(odds, 3), "p_fisher": p_f, "n_wells_treated": len(a), "n_wells_control": len(b), "p_welch_perwell": p_t})
    write_tsv(os.path.join(ANALYSIS, "results", "SMN2_tests.tsv"), tests, list(tests[0].keys()))
    print("\nvs Hek control (BOBseq PE): Fisher exact on pooled junction reads | Welch t on per-well PSI")
    for t in tests:
        if t["unit"] == "dedup": print(f"  {t['paralog']}  {t['contrast']:<32} {t['incl_treated']:>3}/{t['skip_treated']:<3} vs {t['incl_control']:>3}/{t['skip_control']:<3}  dPSI {t['dpsi']:+.2f}  Fisher p = {t['p_fisher']:.2g}  | Welch p = {t['p_welch_perwell']:.2g}")
    print(f"{'set / condition':<36}{'unit':<6}{'n':>3}{'incl':>6}{'skip':>6}{'PSI pooled [CI]':>24}{'SMN2 PSI [CI]':>24}{'SMN1 PSI [CI]':>24}{'T-frac [CI]':>22}{'T+C':>6}"
          f"   SMN2 6>7/7>8/6>8   SMN1 6>7/7>8/6>8   amb 6>7/7>8/6>8")
    for r in per_cond:
        if r["unit"] != "dedup" and r["method"] != "bobseq": continue
        if r["unit"] == "raw" and r["set"] != "bobseq_pe_native": continue
        f = lambda k: f"{r[k]:.2f} [{r[k+'_lo']:.2f}-{r[k+'_hi']:.2f}]" if r[k] != "" else "   -"
        print(f"{(r['set']+' / '+r['condition'])[:36]:<36}{r['unit']:<6}{r['n_samples']:>3}{r['psi_pooled_n']-(r['skip_SMN2']+r['skip_SMN1']+r['skip_ambiguous']):>6}{r['skip_SMN2']+r['skip_SMN1']+r['skip_ambiguous']:>6}{f('psi_pooled'):>24}{f('psi_SMN2'):>24}{f('psi_SMN1'):>24}{f('t_frac'):>22}{r['t_frac_n']:>6}"
              + "".join(f"   {r[f'ex6>7_{p}']:>4}/{r[f'ex7>8_{p}']:>4}/{r[f'ex6>8_{p}']:>4}" for p in ("SMN2", "SMN1", "ambiguous")))

if __name__ == "__main__":
    main()

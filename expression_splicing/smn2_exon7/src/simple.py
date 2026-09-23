"""Stage simple: the BLAST-native exon-7 inclusion index -- reads only, no junctions, no skip attribution.

Every read over exon 7 carries c.840 (T = SMN2, C = SMN1); every read over exon 8 c.*239 carries A = SMN2 /
G = SMN1. Exon 8 is the constitutive last exon, so per paralog  ex7 reads / ex8 reads  is its exon-7 inclusion
(times a coverage-shape factor that is the same for both paralogs in the same well). Normalising to SMN1
(constitutively included) cancels that factor:

    index(SMN2) = (SMN2 ex7 / SMN2 ex8) / (SMN1 ex7 / SMN1 ex8)         ~ SMN2 exon-7 inclusion, SMN1 = 1

95 % CI by parametric bootstrap of the four counts (Poisson). Writes results/SMN2_index_samples.tsv and
results/SMN2_index_conditions.tsv."""
import csv, os, sys
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
KEYS = ["c840_T", "c840_C", "c239_A", "c239_G"]

def index(c, rng=None, n=4000):
    T, C, A, G = (c[k] for k in KEYS)
    if min(T, C, A, G) == 0: return {"index": "", "index_lo": "", "index_hi": "", "ex7_per_ex8_SMN2": "", "ex7_per_ex8_SMN1": ""}
    idx = (T / A) / (C / G)
    b = (rng.poisson(T, n) / rng.poisson(A, n)) / (rng.poisson(C, n) / rng.poisson(G, n)) if rng is not None else np.array([idx])
    b = b[np.isfinite(b)]
    return {"index": round(idx, 3), "index_lo": round(float(np.percentile(b, 2.5)), 3), "index_hi": round(float(np.percentile(b, 97.5)), 3),
            "ex7_per_ex8_SMN2": round(T / A, 4), "ex7_per_ex8_SMN1": round(C / G, 4)}

def main():
    rng = np.random.default_rng(0)
    rows = list(csv.DictReader(open(os.path.join(ANALYSIS, "results", "SMN2_samples.tsv")), delimiter="\t"))
    per = []; pooled = defaultdict(lambda: defaultdict(int)); meta = {}
    for r in rows:
        c = {k: int(r[k]) for k in KEYS}
        per.append({**{k: r[k] for k in ("acc_number", "set", "method", "sample_name", "condition", "in_benchmark", "unit")}, **c, **index(c, rng)})
        grp = (r["set"], r["condition"] if r["method"] == "bobseq" else "all", r["unit"])
        if True:                                     # ALL wells pooled
            for k in KEYS: pooled[grp][k] += c[k]
            m = meta.setdefault(grp, {"set": grp[0], "condition": grp[1], "unit": grp[2], "method": r["method"], "n_samples": 0}); m["n_samples"] += 1
    cond = [{**meta[g], **dict(pooled[g]), **index(pooled[g], rng)} for g in sorted(pooled)]
    write_tsv(os.path.join(ANALYSIS, "results", "SMN2_index_samples.tsv"), per, list(per[0].keys()))
    write_tsv(os.path.join(ANALYSIS, "results", "SMN2_index_conditions.tsv"), cond, list(cond[0].keys()))
    print(f"{'set / condition':<36}{'unit':<6}{'n':>3}{'ex7 T/C':>10}{'ex8 A/G':>12}{'SMN2 ex7/ex8':>13}{'SMN1 ex7/ex8':>13}{'index [95% CI]':>22}")
    for r in cond:
        if r["method"] != "bobseq" and r["unit"] != "dedup": continue
        if r["set"] == "bobseq_50nt" and r["unit"] == "raw": continue
        f = f"{r['index']:.2f} [{r['index_lo']:.2f}-{r['index_hi']:.2f}]" if r["index"] != "" else "-"
        print(f"{(r['set']+' / '+r['condition'])[:36]:<36}{r['unit']:<6}{r['n_samples']:>3}{str(r['c840_T'])+'/'+str(r['c840_C']):>10}{str(r['c239_A'])+'/'+str(r['c239_G']):>12}{r['ex7_per_ex8_SMN2'] or '-':>13}{r['ex7_per_ex8_SMN1'] or '-':>13}{f:>22}")

if __name__ == "__main__":
    main()

"""Stage fragments: type every fragment in the two SMN loci.

Per fragment (read name; both mates together), from the region BAMs (ALL MAPQ, raw and dedup):
  splicing class w.r.t. exon 7 -- 'incl' (ex6>ex7 and/or ex7>ex8 junction), 'skip' (ex6>ex8),
                                  'body' (unspliced over exon 7), 'other'
  paralog call from the BLAST diagnostic sites the fragment covers (either mate, either locus):
                                  'SMN1' / 'SMN2' / 'ambiguous' (covers none) / 'conflict' (mixed)
  c840 base (T = SMN2, C = SMN1) when the fragment covers exon 7 +6.
Junctions and sites are recognised at whichever locus STAR placed the fragment, so a multimapper placed
at the wrong paralog is still classified, and its bases still say which paralog it came from.
Writes results/SMN2_fragments.tsv.gz and results/SMN2_samples.tsv (per sample x unit counts, incl. fragments per
junction ex6>7 / ex7>8 / ex6>8 per paralog)."""
import csv, gzip, os, sys
from collections import defaultdict
import pysam
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

def load_sites():
    sites = {"SMN2": {}, "SMN1": {}}
    with open(os.path.join(ANALYSIS, "results", "SMN_diagnostic_sites.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["type"] != "mismatch": continue
            sites["SMN2"][int(r["SMN2_pos"]) - 1] = (r["SMN2_allele"], r["SMN1_allele"])   # 0-based -> (SMN2 allele, SMN1 allele)
            sites["SMN1"][int(r["SMN1_pos"]) - 1] = (r["SMN2_allele"], r["SMN1_allele"])
    return sites

def locus_of(pos):
    for para, (c, s, e) in LOCI.items():
        if s - 1 <= pos < e: return para
    return None

def type_bam(path, sites):
    frags = defaultdict(lambda: {"incl": 0, "skip": 0, "body": 0, "other": 0, "v2": 0, "v1": 0, "c840": "", "c239": "", "locus": set(), "j67": 0, "j78": 0, "j68": 0})
    for r in pysam.AlignmentFile(path, "rb").fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        para = locus_of(r.reference_start)
        if para is None: continue
        f = frags[r.query_name]; f["locus"].add(para)
        intr = introns(para); e7 = EXONS[para]["ex7"]
        pos = r.reference_start; blocks = []; juncs = []
        for op, ln in r.cigartuples:
            if op in (0, 7, 8): blocks.append((pos, pos + ln)); pos += ln
            elif op == 2: pos += ln
            elif op == 3: juncs.append((pos, pos + ln)); pos += ln
        for j in juncs:
            if j == intr["ex6>ex8"]: f["skip"] += 1; f["j68"] += 1
            elif j == intr["ex6>ex7"]: f["incl"] += 1; f["j67"] += 1
            elif j == intr["ex7>ex8"]: f["incl"] += 1; f["j78"] += 1
            else: f["other"] += 1
        if any(b[0] < e7[1] and b[1] > e7[0] - 1 for b in blocks): f["body"] += 1
        seq = r.query_sequence; st = sites[para]; c840 = C840[para] - 1; c239 = C239[para] - 1
        for qp, rp in r.get_aligned_pairs(matches_only=True):
            if rp in st:
                b = seq[qp]
                if b == st[rp][0]: f["v2"] += 1
                elif b == st[rp][1]: f["v1"] += 1
            if rp == c840: f["c840"] = seq[qp]
            if rp == c239: f["c239"] = seq[qp]
    return frags

def classify(f):
    if f["skip"] and f["incl"]: cls = "conflict_splice"
    elif f["skip"]: cls = "skip"
    elif f["incl"]: cls = "incl"
    elif f["body"]: cls = "body"
    else: cls = "other"
    if f["v2"] and f["v1"]: para = "conflict"
    elif f["v2"]: para = "SMN2"
    elif f["v1"]: para = "SMN1"
    else: para = "ambiguous"
    return cls, para

def main():
    sites = load_sites(); d = os.path.join(ANALYSIS, ".cache", "region_bams")
    samples = load_samples(); out_rows = []; frag_rows = []
    for r in samples:
        acc = r["acc_number"]
        for unit in ("raw", "dedup"):
            frags = type_bam(f"{d}/{acc}.{unit}.smn.bam", sites)
            cnt = defaultdict(int); c840 = defaultdict(int); c239 = defaultdict(int); jn = defaultdict(int)
            for q, f in frags.items():
                cls, para = classify(f); cnt[(cls, para)] += 1
                for j in ("j67", "j78", "j68"):
                    if f[j]: jn[(j, para)] += 1                      # fragments carrying that junction (a fragment can carry j67 and j78)
                if f["c840"] in ("T", "C"): c840[f["c840"]] += 1
                if f["c239"] in ("A", "G"): c239[f["c239"]] += 1
                frag_rows.append((acc, unit, q, "+".join(sorted(f["locus"])), cls, para, f["v2"], f["v1"], f["c840"], f["c239"]))
            row = {"acc_number": acc, "set": r["set"], "method": r["method"], "sample_name": r["sample_name"],
                   "condition": r["condition"], "in_benchmark": r["in_benchmark"], "unit": unit, "fragments": len(frags),
                   "c840_T": c840["T"], "c840_C": c840["C"], "c239_A": c239["A"], "c239_G": c239["G"]}
            for cls in ("incl", "skip", "body", "other", "conflict_splice"):
                for para in ("SMN2", "SMN1", "ambiguous", "conflict"):
                    row[f"{cls}_{para}"] = cnt[(cls, para)]
            for j, lab in (("j67", "ex6>7"), ("j78", "ex7>8"), ("j68", "ex6>8")):
                for para in ("SMN2", "SMN1", "ambiguous", "conflict"):
                    row[f"{lab}_{para}"] = jn[(j, para)]
            out_rows.append(row)
    cols = list(out_rows[0].keys())
    write_tsv(os.path.join(ANALYSIS, "results", "SMN2_samples.tsv"), out_rows, cols)
    with gzip.open(os.path.join(ANALYSIS, "results", "SMN2_fragments.tsv.gz"), "wt") as fh:
        fh.write("acc_number\tunit\tqname\tlocus_placed\tclass\tparalog\tsites_SMN2\tsites_SMN1\tc840\tc239\n")
        for t in frag_rows: fh.write("\t".join(map(str, t)) + "\n")
    tot = defaultdict(int)
    for row in out_rows:
        if row["unit"] != "raw": continue
        for k, v in row.items():
            if isinstance(v, int): tot[k] += v
    print(f"fragments: {len(samples)} samples x 2 units; raw totals: incl SMN2/SMN1/amb/conf = "
          f"{tot['incl_SMN2']}/{tot['incl_SMN1']}/{tot['incl_ambiguous']}/{tot['incl_conflict']}; "
          f"skip = {tot['skip_SMN2']}/{tot['skip_SMN1']}/{tot['skip_ambiguous']}/{tot['skip_conflict']}; "
          f"body = {tot['body_SMN2']}/{tot['body_SMN1']}/{tot['body_ambiguous']}; c840 T/C = {tot['c840_T']}/{tot['c840_C']}")

if __name__ == "__main__":
    main()

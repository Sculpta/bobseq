"""Stage blast_loci: BLAST the SMN1 and SMN2 genomic loci against each other and list every difference.

Writes .cache/blast/{SMN1,SMN2}.fa, results/SMN_locus_blast.tsv (the HSPs) and
results/SMN_diagnostic_sites.tsv: one row per mismatch/indel with the coordinate in BOTH paralogs, the
allele in each, and the feature it sits in. The canonical five (exon 7 +6 c.840 C>T; intron 6 -45;
intron 7 +100, +214; exon 8 c.*239 G>A) must be among them.
"""
import os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

def fetch(para, env):
    ch, s, e = LOCI[para]
    fa = os.path.join(ANALYSIS, ".cache", "blast", f"{para}.fa")
    if not os.path.exists(fa):
        seq = subprocess.run(["samtools", "faidx", FASTA, f"{FASTA_PREFIX}{ch}:{s}-{e}"], capture_output=True, text=True, check=True, env=env).stdout
        lines = seq.split("\n"); open(fa, "w").write(f">{para} {ch}:{s}-{e}\n" + "".join(lines[1:]).upper() + "\n")
    return fa

def feature(para, pos):
    for name, (a, b) in EXONS[para].items():
        if a <= pos <= b: return name
    e = EXONS[para]
    if e["ex6"][1] < pos < e["ex7"][0]: return "intron6"
    if e["ex7"][1] < pos < e["ex8"][0]: return "intron7"
    return "other"

def main():
    env = gcs_env(); d = os.path.join(ANALYSIS, ".cache", "blast")
    q, s = fetch("SMN2", env), fetch("SMN1", env)
    out = subprocess.run(["blastn", "-query", q, "-subject", s, "-dust", "no", "-soft_masking", "false",
                          "-outfmt", "6 qstart qend sstart send length pident bitscore qseq sseq"],
                         capture_output=True, text=True, check=True, env=env).stdout.strip().split("\n")
    hsps = [l.split("\t") for l in out]
    with open(os.path.join(ANALYSIS, "results", "SMN_locus_blast.tsv"), "w") as fh:
        fh.write("qstart\tqend\tsstart\tsend\tlength\tpident\tbitscore\n")
        for h in hsps: fh.write("\t".join(h[:7]) + "\n")
    # walk ONLY the collinear paralog alignment (the top HSP, ~29 kb at 99.9 %); the other HSPs are
    # repeat elements matching each other at the wrong places (Alu etc.) and would fake thousands of sites
    rows = []; q0 = LOCI["SMN2"][1]; s0 = LOCI["SMN1"][1]
    top = max(hsps, key=lambda h: float(h[6]))
    for h in [top]:
        qs, qe, ss, se = map(int, h[:4]); qseq, sseq = h[7], h[8]
        qi, si = qs, ss
        for a, b in zip(qseq, sseq):
            if a != b:
                gq = q0 + qi - 1 if a != "-" else None; gs = s0 + si - 1 if b != "-" else None
                rows.append({"SMN2_pos": gq or "", "SMN1_pos": gs or "", "SMN2_allele": a, "SMN1_allele": b,
                             "type": "mismatch" if "-" not in (a, b) else "indel",
                             "feature_SMN2": feature("SMN2", gq) if gq else "", "feature_SMN1": feature("SMN1", gs) if gs else ""})
            if a != "-": qi += 1
            if b != "-": si += 1
    rows.sort(key=lambda r: (r["SMN2_pos"] or 0))
    write_tsv(os.path.join(ANALYSIS, "results", "SMN_diagnostic_sites.tsv"), rows,
              ["SMN2_pos", "SMN1_pos", "SMN2_allele", "SMN1_allele", "type", "feature_SMN2", "feature_SMN1"])
    c840 = [r for r in rows if r["SMN2_pos"] == C840["SMN2"]]
    assert c840 and c840[0]["SMN2_allele"] == "T" and c840[0]["SMN1_allele"] == "C", "c.840 C>T not found where expected"
    print(f"blast_loci: {len(hsps)} HSP(s), top identity {hsps[0][5]}% over {hsps[0][4]} nt; "
          f"{len(rows)} differences ({sum(r['type']=='mismatch' for r in rows)} mismatches); "
          f"c.840 OK (SMN2 T / SMN1 C); in exon7: {sum(r['feature_SMN2']=='ex7' for r in rows)}, "
          f"exon8: {sum(r['feature_SMN2']=='ex8' for r in rows)}, intron7: {sum(r['feature_SMN2']=='intron7' for r in rows)}, "
          f"intron6: {sum(r['feature_SMN2']=='intron6' for r in rows)}")

if __name__ == "__main__":
    main()

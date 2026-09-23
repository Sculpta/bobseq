"""Stage blast_reads: the literal read-level BLAST cross-check of the site-based paralog calls.

Every raw-fragment read in the two loci is blasted (blastn, megablast) against a 2-sequence subject
(SMN1 + SMN2 loci). Per fragment: best bitscore per paralog over both mates; call SMN1/SMN2 when one
paralog scores higher by >= 2 bits (one mismatch under megablast reward 1 / penalty -2 is ~3 bits), else
'tie'. Writes results/SMN2_blast_reads.tsv (per fragment) and prints the agreement with the
site-based calls of stage `fragments`."""
import csv, gzip, os, subprocess, sys
from collections import defaultdict, Counter
import pysam
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

def main():
    env = gcs_env(); d = os.path.join(ANALYSIS, ".cache"); rb = f"{d}/region_bams"
    subj = f"{d}/blast/SMN_both.fa"
    open(subj, "w").write(open(f"{d}/blast/SMN1.fa").read() + open(f"{d}/blast/SMN2.fa").read())
    q = f"{d}/blast/reads.fa"; n = 0
    with open(q, "w") as fh:
        for r in load_samples():
            for rec in pysam.AlignmentFile(f"{rb}/{r['acc_number']}.raw.smn.bam", "rb").fetch(until_eof=True):
                if rec.is_secondary or rec.is_supplementary or not rec.query_sequence: continue
                seq = rec.get_forward_sequence() or rec.query_sequence
                fh.write(f">{r['acc_number']}|{rec.query_name}|{1 if rec.is_read1 else 2}\n{seq}\n"); n += 1
    out = subprocess.run(["blastn", "-task", "megablast", "-query", q, "-subject", subj, "-dust", "no", "-soft_masking", "false",
                          "-max_hsps", "1", "-outfmt", "6 qseqid sseqid bitscore pident length"], capture_output=True, text=True, check=True, env=env).stdout
    best = defaultdict(lambda: {"SMN1": 0.0, "SMN2": 0.0})
    for line in out.strip().split("\n"):
        if not line: continue
        qid, sid, bits = line.split("\t")[:3]; acc, qn, _ = qid.split("|")
        b = best[(acc, qn)]; b[sid] = max(b[sid], float(bits))
    site = {}
    with gzip.open(os.path.join(ANALYSIS, "results", "SMN2_fragments.tsv.gz"), "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["unit"] == "raw": site[(r["acc_number"], r["qname"])] = r["paralog"]
    rows = []; agree = Counter()
    for (acc, qn), b in best.items():
        call = "SMN2" if b["SMN2"] - b["SMN1"] >= 2 else "SMN1" if b["SMN1"] - b["SMN2"] >= 2 else "tie"
        s = site.get((acc, qn), "?"); agree[(s, call)] += 1
        rows.append({"acc_number": acc, "qname": qn, "bits_SMN1": b["SMN1"], "bits_SMN2": b["SMN2"], "blast_call": call, "site_call": s})
    write_tsv(os.path.join(ANALYSIS, "results", "SMN2_blast_reads.tsv"), rows, ["acc_number", "qname", "bits_SMN1", "bits_SMN2", "blast_call", "site_call"])
    print(f"blast_reads: {n} reads / {len(best)} fragments blasted")
    print(f"{'site-based call':<16}{'BLAST SMN1':>11}{'BLAST SMN2':>11}{'BLAST tie':>10}")
    for s in ("SMN1", "SMN2", "ambiguous", "conflict"):
        print(f"{s:<16}{agree[(s,'SMN1')]:>11}{agree[(s,'SMN2')]:>11}{agree[(s,'tie')]:>10}")

if __name__ == "__main__":
    main()

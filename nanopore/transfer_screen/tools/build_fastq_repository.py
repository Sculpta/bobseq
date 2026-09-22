#!/usr/bin/env python3
"""Build the repository FASTQ folder: one file per dataset scored in
@RULES/bobcode_run_rules.tsv, copied from the original MinKNOW output.

  <this folder>/<chemistry>_<construct>/<YYMMDD>_BC<NN>.fastq.gz

A dataset's original chunk files (and, for 260507, both MinKNOW runs) are joined
byte-for-byte: gzip chunks are concatenated as-is (a valid multi-member gzip, reads
unchanged); plain .fastq sources are gzipped. Nothing is trimmed or filtered, except
that a run with max_reads in the run rules (260120) keeps only its first max_reads
reads (files in name order, reads in file order); the rest are not copied.
Every output is re-read and checked: read count equals the sources' total (or max_reads) and no
read ID occurs twice. manifest.tsv records sources, read count, size and MD5.

Usage (provenance: how the deposited FASTQ folder was built from the MinKNOW output):
  MINKNOW_DATA=/path/to/MinKNOW/data python3 tools/build_fastq_repository.py \
      --data-dir "/path/to/Species Mixing Projects" --out-dir fastq_repository [--force]
"""
import csv, gzip, hashlib, os, shutil, sys

import argparse

RULES_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repository root
sys.path.insert(0, RULES_DIR)
import speciesmix_bobcode_metrics as m               # same input rules as the scorer (§1)

_ap = argparse.ArgumentParser()
_ap.add_argument("--data-dir", required=True, help="folder the run rules' run_folder paths are relative to")
_ap.add_argument("--out-dir", required=True, help="the FASTQ repository folder to write")
_ap.add_argument("--force", action="store_true")
_A = _ap.parse_args()
PROJECT, HERE, FORCE = _A.data_dir, _A.out_dir, _A.force


def repo_name(run, bc):
    chem = run["chemistry"].replace("'", "")
    date = run["run_date"]
    return f"{chem}_{run['construct']}", f"{date}_{bc}.fastq.gz"


def read_ids(path):
    op = gzip.open if path.endswith(".gz") else open
    ids = []
    with op(path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                ids.append(line[1:].split(None, 1)[0])
    return ids


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    rules = os.path.join(RULES_DIR, "bobcode_run_rules.tsv")
    dates = {r["run_folder"]: r["seq_date"].replace("-", "")[2:]
             for r in csv.DictReader(open(rules), delimiter="\t")}
    runs = m.load_run_rules(rules)
    rows, problems = [], []
    for run in runs:
        run["run_date"] = dates[run["run"]]
        for bc in run["barcodes"]:
            sub, name = repo_name(run, bc)
            os.makedirs(os.path.join(HERE, sub), exist_ok=True)
            out = os.path.join(HERE, sub, name)
            srcs = m.input_fastqs(PROJECT, run, bc)
            src_reads = sum(len(read_ids(s)) for s in srcs)
            cap = run["max_reads"]
            expect = min(src_reads, cap) if cap else src_reads
            if FORCE or not os.path.exists(out):
                tmp = out + ".part"
                if cap:                                   # first `cap` reads only
                    seen, kept = set(), 0
                    with gzip.open(tmp, "wt") as fo:
                        for s in srcs:
                            op = gzip.open if s.endswith(".gz") else open
                            with op(s, "rt") as fi:
                                while kept < cap:
                                    rec = [fi.readline() for _ in range(4)]
                                    if not rec[0]:
                                        break
                                    rid = rec[0][1:].split(None, 1)[0]
                                    if rid in seen:
                                        continue
                                    seen.add(rid); kept += 1
                                    fo.write("".join(rec))
                            if kept >= cap:
                                break
                elif all(s.endswith(".gz") for s in srcs):
                    with open(tmp, "wb") as fo:
                        for s in srcs:
                            with open(s, "rb") as fi:
                                shutil.copyfileobj(fi, fo)
                elif not any(s.endswith(".gz") for s in srcs):
                    with gzip.open(tmp, "wb") as fo:
                        for s in srcs:
                            with open(s, "rb") as fi:
                                shutil.copyfileobj(fi, fo)
                else:
                    sys.exit(f"{run['run']} {bc}: mixed gz / plain sources")
                os.replace(tmp, out)
            ids = read_ids(out)
            if len(ids) != expect:
                problems.append(f"{sub}/{name}: {len(ids)} reads written, expected {expect}")
            if len(set(ids)) != len(ids):
                problems.append(f"{sub}/{name}: {len(ids) - len(set(ids))} duplicate read IDs")
            rows.append({
                "file": f"{sub}/{name}", "chemistry": run["chemistry"], "construct": run["construct"],
                "run_folder": run["run"], "ont_barcode": bc, "plot_no": run["plot_no"].get(bc, ""),
                "reads": len(ids), "bytes": os.path.getsize(out), "md5": md5(out),
                "n_source_files": len(srcs),
                "note": f"first {cap:,} of {src_reads:,} reads (rest not used, not uploaded)" if cap else "all reads",
                "source_files": " ; ".join(os.path.relpath(s, PROJECT) if not s.startswith(os.environ.get("MINKNOW_DATA", "\0"))
                                           else "$MINKNOW_DATA/" + os.path.relpath(s, os.environ["MINKNOW_DATA"]) for s in srcs),
            })
            print(f"  {sub}/{name}: {len(srcs)} source file(s), {len(ids):,} reads", flush=True)
    order = {"polydT": 0, "splint": 1, "5'TSO": 2, "3'TSO": 3}           # figure panel order
    rows.sort(key=lambda r: (order.get(r["chemistry"], 9), int(r["plot_no"]) if r["plot_no"] else 10 ** 6))
    with open(os.path.join(HERE, "manifest.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print(f"{len(rows)} files, {sum(r['reads'] for r in rows):,} reads, "
          f"{sum(r['bytes'] for r in rows) / 1e9:.2f} GB; manifest.tsv written")
    if problems:
        print("PROBLEMS:\n  " + "\n  ".join(problems)); sys.exit(1)


if __name__ == "__main__":
    main()

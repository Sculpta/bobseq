"""Stage region_bams: slice BOTH SMN loci out of every per-sample BAM (raw) and its UMI-deduplicated version, ALL MAPQ,
-> .cache/region_bams/<acc>.{raw,dedup}.smn.bam

Raw BAMs are read from the per-sample BAM store of the benchmark (locations.SAMPLE_BAMS, <set>/<sample>.bam). The
deduplicated BAM of a sample is looked up in locations.DEDUP_BAMS (<acc>.dedup.bam) and made there with
benchmark/scripts/dedup_umi_position.py when it does not exist yet."""
import os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

REG = [f"{c}:{s}-{e}" for c, s, e in LOCI.values()]

def dedup_bam(st, smp, acc, env):
    out = DEDUP_BAM.format(acc=acc)
    if not os.path.exists(out + ".bai"):
        os.makedirs(os.path.dirname(out), exist_ok=True)
        subprocess.run([sys.executable, DEDUP_SCRIPT, "--set", st, "--sample", smp, "--in", RAW_BAM.format(set=st, sample=smp), "--out", out],
                       check=True, env=env, capture_output=True)
    return out

def slice_one(job):
    acc, src, out, env = job
    if os.path.exists(out + ".bai"): return acc, "exists"
    tmp = out + ".tmp.bam"
    subprocess.run(["samtools", "view", "-b", "-o", tmp, src] + REG, check=True, env=env, capture_output=True)
    subprocess.run(["samtools", "sort", "-o", out, tmp], check=True, env=env, capture_output=True); os.remove(tmp)
    subprocess.run(["samtools", "index", out], check=True, env=env, capture_output=True)
    return acc, "ok"

def main():
    env = gcs_env(); d = os.path.join(ANALYSIS, ".cache", "region_bams"); os.makedirs(d, exist_ok=True)
    jobs = []
    for r in load_samples():
        acc, st, smp = r["acc_number"], r["set"], r["sample_name"]
        jobs.append((acc, RAW_BAM.format(set=st, sample=smp), f"{d}/{acc}.raw.smn.bam", env))
        jobs.append((acc, dedup_bam(st, smp, acc, env), f"{d}/{acc}.dedup.smn.bam", env))
    with ThreadPoolExecutor(4) as ex:
        res = list(ex.map(slice_one, jobs))
    print(f"region_bams: {sum(s=='ok' for _,s in res)} sliced, {sum(s=='exists' for _,s in res)} existed, of {len(jobs)}")

if __name__ == "__main__":
    main()

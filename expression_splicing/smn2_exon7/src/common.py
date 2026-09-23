"""Shared constants/helpers for the SMN2 exon-7 analysis.

Exon numbering follows the SMN literature (exon 7 = the 54-nt alternatively spliced exon), which is
Ensembl's exon 8 of the canonical transcripts (Ensembl splits the historical exon 2 into 2a/2b):
    lit ex6 = Ens ex7 (111 nt) | lit ex7 = Ens ex8 (54 nt) | lit ex8 = Ens ex9 (577 nt, last).
"""
import csv, os, subprocess, sys
from math import log
from scipy.stats import beta

HERE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.dirname(HERE)
sys.path.insert(0, os.path.dirname(ANALYSIS))
from locations import FASTA, FASTA_PREFIX, SAMPLE_BAMS, DEDUP_BAMS, DEDUP_SCRIPT, SAMPLES
RAW_BAM = os.path.join(SAMPLE_BAMS, "{set}", "{sample}.bam")          # the per-sample BAM store of the benchmark
DEDUP_BAM = os.path.join(DEDUP_BAMS, "{acc}.dedup.bam")               # UMI + position deduplicated (benchmark/scripts/dedup_umi_position.py)

# genomic loci (1-based inclusive, + strand, GRCh38), ~30 kb each, whole gene + margins
LOCI = {"SMN2": ("5", 70049000, 70078600), "SMN1": ("5", 70924500, 70953600)}
# literature exons 6/7/8 = Ensembl canonical exons 7/8/9 (1-based inclusive)
EXONS = {
    "SMN2": {"ex6": (70070641, 70070751), "ex7": (70076521, 70076574), "ex8": (70077019, 70077595)},
    "SMN1": {"ex6": (70946066, 70946176), "ex7": (70951941, 70951994), "ex8": (70952439, 70953015)},
}
# the splicing-relevant site: exon 7 +6 (c.840). SMN2 = T, SMN1 = C
C840 = {"SMN2": 70076526, "SMN1": 70951946}
# exon 8 (constitutive last exon) paralog site c.*239: SMN2 = A, SMN1 = G  (from the locus BLAST)
C239 = {"SMN2": 70077254, "SMN1": 70952674}


def introns(para):
    """0-based half-open intron intervals for the three junctions of interest, per paralog."""
    e = EXONS[para]
    return {"ex6>ex7": (e["ex6"][1], e["ex7"][0] - 1),
            "ex7>ex8": (e["ex7"][1], e["ex8"][0] - 1),
            "ex6>ex8": (e["ex6"][1], e["ex8"][0] - 1)}


def gcs_env():
    """The environment the external tools (samtools, blastn) run in."""
    return dict(os.environ)


def load_samples():
    with open(SAMPLES) as fh:
        return list(csv.DictReader(fh))


def cp_ci(k, n, alpha=0.05):
    """Clopper-Pearson interval for k successes in n."""
    if n == 0: return (float("nan"), float("nan"))
    lo = 0.0 if k == 0 else beta.ppf(alpha / 2, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - alpha / 2, k + 1, n - k)
    return (lo, hi)


def write_tsv(path, rows, cols):
    with open(path + ".tmp", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    os.replace(path + ".tmp", path)

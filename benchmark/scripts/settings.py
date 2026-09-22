"""Locations and frozen parameters shared by every benchmark script.

Every script in this directory imports its paths from here; nothing else is hard-coded. Override a location with the
environment variable named beside it.

    BOBSEQ_WORK      working directory of the benchmark (default <repo>/work); the arms, per-sample BAMs, coverage tables
                     and generated figures live under it, in the layout described in benchmark/README.md
    BOBSEQ_REF_DIR   the combined GRCh38.113 + GRCm39.113 reference built by reference/build_reference.sh
"""
import os

CODE = os.path.dirname(os.path.abspath(__file__))          # benchmark/scripts
ROOT = os.path.dirname(os.path.dirname(CODE))               # the repository
CONFIG = os.path.join(ROOT, 'config')
PIPELINE_DIR = os.path.join(ROOT, 'pipeline', 'scripts')    # the BOBseq Illumina processing modules

WORK = os.environ.get('BOBSEQ_WORK', os.path.join(ROOT, 'work'))
PROJECT = WORK                                              # data/ (staged inputs) and results/ sit beside the arms
UNIFORM = os.path.join(WORK, 'benchmark_uniform')           # one directory per arm
COVERAGE = os.path.join(WORK, 'coverage_architecture')      # per-sample position tables and the coverage sample sheet
FIGURES_NATIVE = os.path.join(WORK, 'preprint_figures_native')
FIGURES_50NT = os.path.join(WORK, 'preprint_figures_50nt')

REF_DIR = os.environ.get('BOBSEQ_REF_DIR', os.path.join(ROOT, 'reference', 'combined_GRCh38.113_GRCm39.113'))
GTF = os.path.join(REF_DIR, 'combined_genome.gtf')
REFFLAT = os.path.join(REF_DIR, 'combined_genome.refflat')
FASTA = os.path.join(REF_DIR, 'combined_genome.fa')
RDNA_BED = os.path.join(CONFIG, 'rdna_loci.bed')
RUN_JSON = os.path.join(CONFIG, 'run_24plex.json')          # the 24-plex run: bobcodes, sample labels, species

# Frozen benchmark parameters (the same values as settings.sh).
READLEN = 50            # nt kept from the cDNA read in the read-length matched arms (the minimum across methods)
UMI_MAXED = 1           # UMI edit distance tolerated within one molecule (definition D); 0 for the 6-nt BOBseq UMI
POS_TOL = 0             # bp window for definition D: 0 = reads must start at the same position
THREADS = int(os.environ.get('BOBSEQ_THREADS', 8))

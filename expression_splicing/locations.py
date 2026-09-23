"""Locations shared by the analyses in this directory. Everything that is not part of the repository is found through this
module and can be overridden with an environment variable.

  BOBSEQ_TABLES         the supplementary data tables of the preprint (S1_sample_metadata.csv to
                        S11_gene_counts_featurecounts.csv.gz), unpacked into one directory
                        default: expression_splicing/data/supplement_tables
  BOBSEQ_HUMAN_GTF      the Ensembl 113 human annotation, Homo_sapiens.GRCh38.113.gtf (plain or gzip)
                        default: the file reference/build_reference.sh downloads
  BOBSEQ_GENOME_FASTA   an indexed genome FASTA that holds GRCh38 (the SMN loci are cut from it)
                        default: the combined genome of reference/build_reference.sh, whose human contigs carry the prefix
                        BOBSEQ_GENOME_PREFIX (HUMAN_; set it to an empty string for a plain GRCh38 FASTA)
  BOBSEQ_SAMPLE_BAMS    the per-sample BAM store of the benchmark, <set>/<sample>.bam (benchmark/build_per_sample_bams.sh)
                        default: work/benchmark_uniform/per_sample_bams
  BOBSEQ_DEDUP_BAMS     where UMI-deduplicated per-sample BAMs are kept, <set>__<sample>.dedup.bam; made on demand with
                        benchmark/scripts/dedup_umi_position.py when missing
                        default: expression_splicing/data/dedup_bams
"""
import gzip
import os

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
REF = os.path.join(REPO, 'reference', 'combined_GRCh38.113_GRCm39.113')
TABLES = os.environ.get('BOBSEQ_TABLES', os.path.join(HERE, 'data', 'supplement_tables'))
GTF = os.environ.get('BOBSEQ_HUMAN_GTF', os.path.join(REF, 'dl', 'Homo_sapiens.GRCh38.113.gtf.gz'))
FASTA = os.environ.get('BOBSEQ_GENOME_FASTA', os.path.join(REF, 'combined_genome.fa'))
FASTA_PREFIX = os.environ.get('BOBSEQ_GENOME_PREFIX', '' if 'BOBSEQ_GENOME_FASTA' in os.environ else 'HUMAN_')
SAMPLE_BAMS = os.environ.get('BOBSEQ_SAMPLE_BAMS', os.path.join(REPO, 'work', 'benchmark_uniform', 'per_sample_bams'))
DEDUP_BAMS = os.environ.get('BOBSEQ_DEDUP_BAMS', os.path.join(HERE, 'data', 'dedup_bams'))
DEDUP_SCRIPT = os.path.join(REPO, 'benchmark', 'scripts', 'dedup_umi_position.py')
SAMPLES = os.path.join(HERE, 'sample_metadata.csv')      # the 119 sample-arms of the count matrices (= S1_sample_metadata.csv)


def table(prefix):
    """Path of the supplementary data table whose file name starts with <prefix>_, e.g. table('S6')."""
    hits = sorted(f for f in os.listdir(TABLES) if f.startswith(prefix + '_')) if os.path.isdir(TABLES) else []
    if not hits:
        raise SystemExit(f'supplementary data table {prefix}_* not found in {TABLES} (unpack the tables there or set BOBSEQ_TABLES)')
    return os.path.join(TABLES, hits[0])


def open_text(path):
    """Open a plain or gzip-compressed text file for reading."""
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)

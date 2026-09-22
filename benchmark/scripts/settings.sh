# Frozen benchmark settings. Every arm is processed through these and nothing else; changing a value here invalidates
# every arm, so re-run all of them together. The Python scripts read the same values from settings.py.
HERE_SETTINGS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BM_ROOT="$(dirname "$(dirname "$HERE_SETTINGS")")"
BM_CODE="$HERE_SETTINGS"
BM_CONFIG="$BM_ROOT/config"
BM_PIPELINE="$BM_ROOT/pipeline"
BM_WORK="${BOBSEQ_WORK:-$BM_ROOT/work}"
BM_UNIFORM="$BM_WORK/benchmark_uniform"
BM_COVERAGE="$BM_WORK/coverage_architecture"
BM_REF="${BOBSEQ_REF_DIR:-$BM_ROOT/reference/combined_GRCh38.113_GRCm39.113}"
BM_RUN_JSON="$BM_CONFIG/run_24plex.json"
BM_THREADS="${BOBSEQ_THREADS:-8}"
export BOBSEQ_WORK="$BM_WORK" BOBSEQ_REF_DIR="$BM_REF" PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg

BM_READLEN=50                     # nt kept from the cDNA read in the read-length matched arms: the minimum across methods
                                  #   DRUG-seq 52 -> 50, prime-seq 50 -> 50, BOBseq read 2 (150) -> 50

# STAR. Identical for every arm; do not add an arm-specific flag here.
BM_STAR_PARAMS="--outSAMtype BAM Unsorted \
--outSAMattributes NH HI AS nM \
--outSAMunmapped Within \
--outSAMmultNmax 18446744073709551615 \
--outFilterMultimapNmax 10000 \
--outFilterScoreMinOverLread 0.66 \
--outFilterMatchNminOverLread 0.66 \
--outFilterMismatchNmax 10 \
--outFilterMismatchNoverLmax 0.04 \
--alignIntronMax 1000000 \
--alignSJDBoverhangMin 3 \
--alignEndsType Local"

# Molecule definition (applied downstream, identical for every arm).
BM_UMI_MAXED=1                    # UMI edit distance tolerated within one molecule (definition D); an arm may override it in arm.env
BM_POS_TOL=0                      # bp window: 0 = reads must start at the SAME position
#
# Why 0 and not a window. These are 3'-end counting methods: reverse transcription primes at the poly(A) site, so every
# cDNA from one RNA starts at the same 3' end; two reads of one well with the same UMI at the same position are one
# molecule, at different positions two. The libraries are tagmented after cDNA synthesis, so a read start is where the
# transposase cut and sibling fragments of one cDNA start at different positions and are counted separately under this
# definition. Molecules under this definition are capture events, not RNA molecules. The +/-1 kb variant (definition E)
# is computed alongside so the size of the choice stays visible.

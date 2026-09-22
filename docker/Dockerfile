# Container with every tool the repository uses: the BOBseq Illumina pipeline, the benchmark and the figures.
#   docker build -t bobseq:0.1 -f docker/Dockerfile .
# Tools: STAR 2.7.11b (the official static build), samtools, Picard, deepTools (bamCoverage), pigz, mawk, gcc,
# Python 3.12 with pysam, numpy, scipy, matplotlib, reportlab, openpyxl, adjustText. Liberation Sans is registered
# under the name Arial for matplotlib (docker/sitecustomize.py), so the figures render with the metric-compatible face.
# The image carries tools only; the repository and the data are mounted at run time.
FROM debian:12-slim
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl unzip ca-certificates procps bzip2 gzip pigz mawk fonts-liberation \
    && rm -rf /var/lib/apt/lists/*
# STAR: the official precompiled static binary of the release
ENV STAR_VERSION=2.7.11b
RUN wget -q -O /tmp/STAR.zip https://github.com/alexdobin/STAR/archive/refs/tags/${STAR_VERSION}.zip && \
    unzip -q /tmp/STAR.zip -d /tmp && \
    cp /tmp/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin/STAR && chmod +x /usr/local/bin/STAR && \
    rm -rf /tmp/STAR.zip /tmp/STAR-${STAR_VERSION}
# conda-forge / bioconda tools through micromamba, one environment from environment.yml
ENV MAMBA_ROOT_PREFIX=/opt/conda
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /usr/local bin/micromamba
COPY environment.yml /tmp/environment.yml
RUN micromamba create -y -n bobseq -f /tmp/environment.yml && micromamba clean -a -y
ENV PATH="/opt/conda/envs/bobseq/bin:$PATH"
# Arial -> Liberation Sans alias for matplotlib (the figure scripts ask for Arial by name)
COPY docker/sitecustomize.py /opt/conda/envs/bobseq/lib/python3.12/site-packages/sitecustomize.py
# every tool must answer, or the build fails
RUN STAR --version && samtools --version | head -1 && picard CollectRnaSeqMetrics --version 2>&1 | tail -1 && \
    bamCoverage --version && pigz --version && mawk -W version 2>&1 | head -1 && gcc --version | head -1 && \
    python3 -c "import pysam, numpy, scipy, matplotlib, reportlab, openpyxl, adjustText; print('python deps ok')" && \
    python3 -c "import matplotlib; matplotlib.use('Agg'); from matplotlib import font_manager as fm, rcParams; rcParams['font.family']='Arial'; f=fm.findfont('Arial'); assert 'Liberation' in f, f; print('font alias ok:', f.split('/')[-1])"
WORKDIR /work

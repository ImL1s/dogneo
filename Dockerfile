# Dog mRNA SOS — Multi-stage Dockerfile
# FOR RESEARCH USE ONLY

# ==========================================================================
# Stage 1: Base bioinformatics tools
# ==========================================================================
FROM ubuntu:22.04 AS biotools

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl ca-certificates unzip \
    build-essential cmake zlib1g-dev libbz2-dev liblzma-dev \
    libncurses5-dev libcurl4-openssl-dev libssl-dev \
    openjdk-17-jre-headless \
    python3 python3-pip python3-venv \
    && rm -rf /var/lib/apt/lists/*

# BWA
RUN wget -q https://github.com/lh3/bwa/releases/download/v0.7.18/bwa-0.7.18.tar.bz2 \
    && tar xjf bwa-0.7.18.tar.bz2 && cd bwa-0.7.18 && make -j$(nproc) \
    && cp bwa /usr/local/bin/ && cd / && rm -rf bwa-*

# SAMtools
RUN wget -q https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 \
    && tar xjf samtools-1.20.tar.bz2 && cd samtools-1.20 \
    && ./configure && make -j$(nproc) && make install \
    && cd / && rm -rf samtools-*

# STAR
RUN wget -q https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz \
    && tar xzf 2.7.11b.tar.gz && cp STAR-2.7.11b/bin/Linux_x86_64/STAR /usr/local/bin/ \
    && rm -rf STAR-* 2.7.11b.tar.gz

# Salmon
RUN wget -q https://github.com/COMBINE-lab/salmon/releases/download/v1.10.3/salmon-1.10.3_linux_x86_64.tar.gz \
    && tar xzf salmon-*.tar.gz && cp salmon-*/bin/salmon /usr/local/bin/ \
    && cp -r salmon-*/lib/* /usr/local/lib/ && rm -rf salmon-*

# GATK
RUN wget -q https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip \
    && unzip -q gatk-*.zip && mv gatk-4.5.0.0 /opt/gatk \
    && ln -s /opt/gatk/gatk /usr/local/bin/gatk && rm -f gatk-*.zip

# SnpEff
RUN wget -q https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip \
    && unzip -q snpEff_latest_core.zip -d /opt/ && rm -f snpEff_latest_core.zip
ENV PATH="/opt/snpEff/exec:${PATH}"


# ==========================================================================
# Stage 2: Python application
# ==========================================================================
FROM biotools AS app

WORKDIR /app

# Copy project files
COPY pyproject.toml README.md LICENSE ./
COPY dogneo/ ./dogneo/

# Install Python dependencies
RUN python3 -m pip install --no-cache-dir --upgrade pip \
    && python3 -m pip install --no-cache-dir -e ".[dev]"

# Create non-root user
RUN groupadd -r dogneo && useradd -r -g dogneo -d /app dogneo \
    && chown -R dogneo:dogneo /app
USER dogneo

# Default results directory
RUN mkdir -p /app/results

ENTRYPOINT ["dogneo"]
CMD ["--help"]

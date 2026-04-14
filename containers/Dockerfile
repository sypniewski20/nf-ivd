FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# =========================
# system deps
# =========================
RUN apt-get update && apt-get install -y \
    wget curl git unzip build-essential \
    zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev \
    libcurl4-openssl-dev libssl-dev \
    python3 python3-pip python3-dev \
    openjdk-17-jre perl \
    && rm -rf /var/lib/apt/lists/*

# =========================
# copy manifest into image
# =========================
COPY build_manifest.tsv /opt/build-manifest.tsv

# =========================
# helper: extract values from TSV
# =========================
RUN apt-get update && apt-get install -y gawk && rm -rf /var/lib/apt/lists/*

# =========================
# SAMTOOLS
# =========================
RUN TOOL=samtools && \
    URL=$(awk -F'\t' -v t=$TOOL '$1==t {print $3}' /opt/build-manifest.tsv) && \
    SHA=$(awk -F'\t' -v t=$TOOL '$1==t {print $4}' /opt/build-manifest.tsv) && \
    FILE=$(basename $URL) && \
    cd /opt && \
    wget -O $FILE $URL && \
    echo "$SHA  $FILE" | sha256sum -c - && \
    tar -xjf $FILE && \
    cd samtools-* && \
    ./configure --prefix=/usr/local && make -j 4 && make install

# =========================
# HTSLIB
# =========================
RUN TOOL=htslib && \
    URL=$(awk -F'\t' -v t=$TOOL '$1==t {print $3}' /opt/build-manifest.tsv) && \
    SHA=$(awk -F'\t' -v t=$TOOL '$1==t {print $4}' /opt/build-manifest.tsv) && \
    FILE=$(basename $URL) && \
    cd /opt && \
    wget -O $FILE $URL && \
    echo "$SHA  $FILE" | sha256sum -c - && \
    tar -xjf $FILE && \
    cd htslib-* && \
    ./configure --prefix=/usr/local && make -j 4 && make install

# =========================
# BWA
# =========================
RUN TOOL=bwa && \
    URL=$(awk -F'\t' -v t=$TOOL '$1==t {print $3}' /opt/build-manifest.tsv) && \
    SHA=$(awk -F'\t' -v t=$TOOL '$1==t {print $4}' /opt/build-manifest.tsv) && \
    FILE=$(basename $URL) && \
    cd /opt && \
    wget -O $FILE $URL && \
    echo "$SHA  $FILE" | sha256sum -c - && \
    tar -xzf $FILE && \
    cd bwa-* && \
    make -j 4 && cp bwa /usr/local/bin/

# =========================
# DRAGMAP
# =========================
RUN TOOL=dragmap && \
    URL=$(awk -F'\t' -v t=$TOOL '$1==t {print $3}' /opt/build-manifest.tsv) && \
    SHA=$(awk -F'\t' -v t=$TOOL '$1==t {print $4}' /opt/build-manifest.tsv) && \
    FILE=$(basename $URL) && \
    cd /opt && \
    wget -O $FILE $URL && \
    echo "$SHA  $FILE" | sha256sum -c - && \
    tar -xzf $FILE && \
    cd DRAGMAP-* && \
    make -j 4 && cp build/dragmap /usr/local/bin/ || true

# =========================
# GATK
# =========================
RUN TOOL=gatk && \
    URL=$(awk -F'\t' -v t=$TOOL '$1==t {print $3}' /opt/build-manifest.tsv) && \
    SHA=$(awk -F'\t' -v t=$TOOL '$1==t {print $4}' /opt/build-manifest.tsv) && \
    FILE=$(basename $URL) && \
    cd /opt && \
    wget -O $FILE $URL && \
    echo "$SHA  $FILE" | sha256sum -c - && \
    unzip $FILE && \
    mv gatk-* /opt/gatk

ENV PATH="/opt/gatk:${PATH}"

# =========================
# python tools (stable pip layer)
# =========================
RUN pip3 install --upgrade pip && \
    pip3 install multiqc numpy pandas pysam cython matplotlib

CMD ["/bin/bash"]
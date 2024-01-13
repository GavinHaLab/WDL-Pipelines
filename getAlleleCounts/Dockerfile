# build as ../titancna:v1.17.1
FROM rocker/r-base:4.2.2

# Install dependencies for building from source
RUN apt-get update && \
    apt-get install -y \
    build-essential \
    libncurses5-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    wget \
    ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Download and install samtools v1.3.1
RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
    tar -xjf samtools-1.3.1.tar.bz2 && \
    cd samtools-1.3.1 && \
    ./configure --without-curses && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-1.3.1.tar.bz2 samtools-1.3.1

# Download and install bcftools v1.2
RUN wget https://github.com/samtools/bcftools/releases/download/1.2/bcftools-1.2.tar.bz2 && \
    tar -xjf bcftools-1.2.tar.bz2 && \
    cd bcftools-1.2 && \
    make && \
    make install && \
    cd .. && \
    rm -rf bcftools-1.2.tar.bz2 bcftools-1.2
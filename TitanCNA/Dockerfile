# build as ../titancna:v1.17.1
FROM rocker/r-base:4.2.2

# Install system dependencies
RUN apt-get update --fix-missing -y && \
    apt-get install -y cmake libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev && \
    rm -rf /var/lib/apt/lists/*

# Install remotes package
RUN R -e "install.packages('remotes', dependencies = TRUE, repos = 'http://cran.rstudio.com/')"

#Install with remotes to have Titan as an R package
RUN R -e "remotes::install_github('gavinha/TitanCNA', dependencies = TRUE)"

# Download Titan to have scripts
RUN wget https://github.com/gavinha/TitanCNA/archive/refs/tags/v1.17.1.zip && \
    unzip v1.17.1.zip && \
    mv TitanCNA-1.17.1/ TitanCNA/ && \
    rm v1.17.1.zip

# Override the R CMD for workflow managers
CMD ["/bin/bash"]


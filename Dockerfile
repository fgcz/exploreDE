FROM rocker/shiny-verse:4.4.3

# Install system dependencies (including python3-dev)
RUN apt-get update && apt-get install -y --no-install-recommends \
    locales git patch cmake ldap-utils \
    python3-pip python3-setuptools python3-wheel python3-venv python3-dev libmysqlclient-dev \
    libgit2-dev libicu-dev libpq-dev libglpk-dev libhdf5-dev libboost-dev libboost-system-dev libboost-filesystem-dev libomp-dev libgsl-dev \
  && rm -rf /var/lib/apt/lists/*

# Set the locale
RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8

# Create the app directory
WORKDIR /srv/shiny-server/exploreDE

# Copy the app files
COPY /inst/app/* /srv/shiny-server/

# Shiny Server runs as user 'shiny'
RUN chown -R shiny:shiny /srv/shiny-server

# Copy the custom shiny-server.conf
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

# Prepare shiny user that can access gstore and so:
RUN usermod -u 5008 shiny
RUN groupadd -g 9951 SG_DockerTRX
RUN usermod -a -G SG_DockerTRX shiny
RUN groupadd -g 10147 SG_Employees
RUN usermod -a -G SG_Employees shiny

# b-fabric certificate so the user within docker can access the ldap information
COPY chain.pem /usr/local/share/ca-certificates/ca.crt
ENV LDAPTLS_CACERT=/usr/local/share/ca-certificates/ca.crt
RUN update-ca-certificates

# Install R packages
RUN R -e 'install.packages(c("sortable", "waiter", "ggprism", "rstatix", "shinylogs", "fresh", "shinyBS", "scico", "chameleon", "pals", "qs", "qs2", "shinyalert", "openxlsx", "shinydashboard", "DT", "colourpicker", "writexl", "circlize", "remotes", "rsconnect", "jsonlite", "hdf5r", "ggpubr", "ggrepel", "ggbeeswarm", "GGally", "patchwork", "Matrix"), repos = "https://cloud.r-project.org/"); \
  BiocManager::install(c("DESeq2", "pathview", "reshape2", "vsn", "Rsubread", "preprocessCore", "wesanderson", "RCurl", "caTools", "matrixStats", "Repitools", "htmltools", "biomaRt", "grid", "gridExtra", "RColorBrewer", "WGCNA", "plyr", "pvclust", "parallel", "Biostrings", "Rsamtools", "Hmisc", "XML", "stringr", "GenomicAlignments", "GenomicFeatures", "GenomicRanges", "ShortRead", "Gviz", "gplots", "GO.db","GOstats", "annotate", "bitops", "edgeR", "limma", "S4Vectors", "VariantAnnotation", "rmarkdown", "plotly", "scran", "data.table", "kableExtra", "htmlwidgets" ,"webshot", "clusterProfiler", "dupRadar", "pheatmap", "taxize", "SingleCellExperiment", "SummarizedExperiment", "scater", "DropletUtils", "shiny", "heatmaply", "readxl","readr", "dplyr", "shinycssloaders", "shinyjs", "slingshot","Rmagic", "reticulate", "viridis", "Seurat", "tidyverse", "httr", "jsonlite", "xml2", "zip", "pcaMethods", "ComplexHeatmap"),update = FALSE, ask = FALSE)'
RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/Rmagic/Rmagic_2.0.3.tar.gz"); remotes::install_github("uzh/ezRun", dependencies = F); remotes::install_github("fgcz/exploreDE", dependencies = F)'

# Create Python virtual environment
RUN python3 -m venv /srv/shiny-server/venv
ENV PATH="/srv/shiny-server/venv/bin:$PATH"
# Install Python dependencies
RUN /srv/shiny-server/venv/bin/pip install --no-cache-dir wheel setuptools numpy cython

# Listening port
EXPOSE 3838

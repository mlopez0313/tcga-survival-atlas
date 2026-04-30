# =============================================================================
# Minimal Dockerfile for the TCGA Survival Atlas Shiny app.
#
# This image bundles only the app code. It expects the data cache and the
# cohort root to be mounted in at runtime (the data is too large to bake in
# and most of it lives on a research filesystem):
#
#   docker build -t tcga-survival-atlas .
#   docker run --rm -p 3838:3838 \
#       -v /path/to/data_cache:/srv/cache:ro \
#       -v /path/to/tcga_survival_cohorts:/srv/cohorts:ro \
#       -e TCGA_SURVIVAL_CACHE_DIR=/srv/cache \
#       -e TCGA_SURVIVAL_COHORT_ROOT=/srv/cohorts \
#       tcga-survival-atlas
#
# For an air-gapped image, run shiny_app/scripts/prepare_app_data.R against
# the source data first, then COPY the resulting data_cache/ here and skip
# the cache volume mount.
# =============================================================================

FROM rocker/shiny:4.4.1

RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
        libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libharfbuzz-dev \
        libfribidi-dev libgit2-dev pandoc && \
    rm -rf /var/lib/apt/lists/*

RUN R --no-save -e "install.packages(c( \
        'shiny','bslib','bsicons','htmltools','data.table','ggplot2', \
        'jsonlite','DT','plotly','shinyWidgets','base64enc','scales', \
        'ggrepel','RColorBrewer','viridis','pheatmap','patchwork','cowplot', \
        'matrixStats','cluster','factoextra','uwot','survival','survminer', \
        'thematic','dplyr','tidyr','stringr','glue','fs','arrow' \
    ), repos = 'https://packagemanager.posit.co/cran/__linux__/jammy/latest')"

WORKDIR /srv/shiny-server/tcga-survival-atlas
COPY . /srv/shiny-server/tcga-survival-atlas/

ENV TCGA_SURVIVAL_APP_ROOT=/srv/shiny-server/tcga-survival-atlas \
    TCGA_SURVIVAL_CACHE_DIR=/srv/cache \
    TCGA_SURVIVAL_COHORT_ROOT=/srv/cohorts \
    TCGA_SURVIVAL_LOG_LEVEL=info \
    PORT=3838 \
    HOST=0.0.0.0

EXPOSE 3838

# Healthcheck calls into the same validator the in-app Diagnostics card uses.
HEALTHCHECK --interval=30s --timeout=15s --start-period=20s --retries=3 \
    CMD Rscript /srv/shiny-server/tcga-survival-atlas/scripts/sanity_check.R || exit 1

CMD ["Rscript", "/srv/shiny-server/tcga-survival-atlas/scripts/run_app.R"]

#!/usr/bin/env bash

## Create environment w/ conda ## %%%%%%%%%%%%%%%%%%%%%%%%%%
# Check if environment exists in a conditinonal statement
if [[ -z $(conda env list | grep r_env) ]]; then
  mamba create -y -n r_env r-essentials r-base
fi
mamba activate r_env

## Set up renv ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rscript -e "
userdir <- unlist(strsplit(Sys.getenv('R_LIBS_USER'),':'))[1L]
repoloc <- 'https://cran.r-project.org/'
install.packages('renv', lib = userdir, repos = repoloc)
"
Rscript -e "renv::init()"

# Install specific packages ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Rscript -e "install.packages('NMF', repos='http://R-Forge.R-project.org')"
Rscript -e "
source('http://renozao.github.io/repotools/install.R')
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')
BiocManager::install('Biobase')
repotools::install.pkgs('NMF', force=FALSE, quiet=FALSE)
"

# visutils and Seurat is a dependency
Rscript -e "
if (requireNamespace('devtools', quietly = TRUE)) {
  userdir <- unlist(strsplit(Sys.getenv('R_LIBS_USER'),':'))[1L]
  repoloc <- 'https://cran.r-project.org/'
  install.packages('devtools', lib = userdir, repos = repoloc)
}
install.packages("fftwtools")
BiocManager::install("EBImage")
devtools::install_github('iaaka/visutils')
"

# Rscript -e "
# if (requireNamespace('remotes', quietly = TRUE)) {
#   userdir <- unlist(strsplit(Sys.getenv('R_LIBS_USER'),':'))[1L]
#   repoloc <- 'https://cran.r-project.org/'
#   install.packages('remotes', lib = userdir, repos = repoloc)
# }
# remotes::install_github('satijalab/seurat', 'seurat5', quiet = TRUE)
# "
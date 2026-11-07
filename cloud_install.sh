#!/usr/bin/env bash
set -euo pipefail

# Codex Cloud Linux environment install script for the Erlang studies.
if command -v sudo >/dev/null 2>&1; then
  sudo apt-get update
  sudo env DEBIAN_FRONTEND=noninteractive apt-get install -y \
    r-base r-base-dev build-essential libcurl4-openssl-dev libssl-dev \
    libxml2-dev libblas-dev liblapack-dev
else
  apt-get update
  DEBIAN_FRONTEND=noninteractive apt-get install -y \
    r-base r-base-dev build-essential libcurl4-openssl-dev libssl-dev \
    libxml2-dev libblas-dev liblapack-dev
fi

Rscript -e 'install.packages(c("foreach", "doParallel", "doSNOW", "RhpcBLASctl", "Rcpp", "RcppArmadillo", "distributions3", "extraDistr", "magrittr", "purrr", "rlist", "assertthat", "fields"), repos="https://cloud.r-project.org", Ncpus=2)'
Rscript -e 'for (p in c("foreach", "doParallel", "doSNOW", "RhpcBLASctl", "Rcpp", "RcppArmadillo", "distributions3", "extraDistr", "magrittr", "purrr", "rlist", "assertthat", "fields")) library(p, character.only=TRUE)'

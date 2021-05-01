#!/usr/bin/env bash
LOPASS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

## setup conda
# __conda_setup="$('/opt/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
# eval "$__conda_setup"
conda config --add channels bioconda
conda config --add channels conda-forge
sudo conda env create -n lopass --file conda-lopass.yml

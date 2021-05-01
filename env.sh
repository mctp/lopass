#!/usr/bin/env bash
LOPASS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export LOPASS_DIR

conda activate lopass
export PATH=$LOPASS_DIR/bin:$PATH

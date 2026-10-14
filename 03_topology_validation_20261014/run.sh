#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export JULIA_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
if [ "$#" -eq 0 ]; then
  echo 'Usage: ./run.sh inputs/NEW_CASE.toml [inputs/ANOTHER_NEW_CASE.toml ...]'
  echo 'Copy a saved TOML and change name first; existing run directories are refused.'
  exit 2
fi
exec julia --project=. simulate.jl "$@"

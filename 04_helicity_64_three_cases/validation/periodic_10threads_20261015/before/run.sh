#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [ "$#" -lt 1 ]; then
  echo 'Usage: ./run.sh {single_ring|hopf_link|trefoil} [--check|--resume]'
  exit 2
fi
case "$1" in single_ring|hopf_link|trefoil) case_name="$1";; *) echo 'Unknown case' >&2; exit 2;; esac
shift
case "${VORTEX_THREADS:-2}" in [1-4]) ;; *) echo 'VORTEX_THREADS must be 1, 2, 3, or 4' >&2; exit 2;; esac
export JULIA_NUM_THREADS="${VORTEX_THREADS:-2}" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
exec julia --startup-file=no --project=. "${case_name}.jl" "$@"

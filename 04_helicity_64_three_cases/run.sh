#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [ "$#" -lt 1 ]; then
  echo 'Usage: ./run.sh {single_ring|hopf_link|trefoil} [--check|--resume]'
  exit 2
fi
case "$1" in single_ring|hopf_link|trefoil) case_name="$1";; *) echo 'Unknown case' >&2; exit 2;; esac
shift
threads="${VORTEX_THREADS:-10}"
case "$threads" in ''|*[!0-9]*|0*) echo 'VORTEX_THREADS must be a positive integer without leading zeros' >&2; exit 2;; esac
available_cpus="$(getconf _NPROCESSORS_ONLN)"
if [ "$threads" -gt "$available_cpus" ]; then
  echo "VORTEX_THREADS=$threads exceeds $available_cpus available logical CPUs" >&2
  exit 2
fi
# Run cases sequentially; BLAS/OpenMP do not create nested worker pools.
export JULIA_NUM_THREADS="$threads" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
exec julia --startup-file=no --project=. "${case_name}.jl" "$@"

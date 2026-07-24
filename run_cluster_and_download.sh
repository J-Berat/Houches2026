#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIRECTORY="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CLUSTER_HOST="${CLUSTER_HOST:-jberat@sr650node230}"

printf '\nDYNAMO — REMOTE COMPUTATION\n'
printf 'Cluster: %s\n\n' "${CLUSTER_HOST}"

ssh -t "${CLUSTER_HOST}" \
    'cd "$HOME/Houches2026" && julia --threads=auto --startup-file=no --project=. run_figures.jl'

printf '\nRemote computation completed successfully.\n'
printf 'Downloading the generated figures to the laptop...\n'

CLUSTER_HOST="${CLUSTER_HOST}" \
    bash "${SCRIPT_DIRECTORY}/download_figures.sh"

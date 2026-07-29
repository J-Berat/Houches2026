#!/usr/bin/env bash
set -euo pipefail

repository_directory="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "${repository_directory}"

exec julia --threads=auto --startup-file=no --project=. \
    run_figures.jl --interactive "$@"

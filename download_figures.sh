#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIRECTORY="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

CLUSTER_HOST="${CLUSTER_HOST:-PSMN_sr650node230}"
REMOTE_FIGURE_DIRECTORY="${REMOTE_FIGURE_DIRECTORY:-~/Houches2026/figures}"
LOCAL_FIGURE_DIRECTORY="${LOCAL_FIGURE_DIRECTORY:-${SCRIPT_DIRECTORY}/figures}"

if ! command -v rsync >/dev/null 2>&1; then
    printf 'Error: rsync is not installed on this laptop.\n' >&2
    exit 1
fi

mkdir -p "${LOCAL_FIGURE_DIRECTORY}"

printf '\nDYNAMO — DOWNLOAD FIGURES\n'
printf 'Cluster          : %s\n' "${CLUSTER_HOST}"
printf 'Remote directory : %s\n' "${REMOTE_FIGURE_DIRECTORY}"
printf 'Local directory  : %s\n\n' "${LOCAL_FIGURE_DIRECTORY}"

rsync \
    --archive \
    --compress \
    --partial \
    --progress \
    "${CLUSTER_HOST}:${REMOTE_FIGURE_DIRECTORY%/}/" \
    "${LOCAL_FIGURE_DIRECTORY%/}/"

FIGURE_COUNT="$(
    find "${LOCAL_FIGURE_DIRECTORY}" -type f \
        \( -name '*.png' -o -name '*.pdf' \) |
        wc -l |
        tr -d ' '
)"

printf '\nDOWNLOAD COMPLETE\n'
printf 'Downloaded figures: %s\n' "${FIGURE_COUNT}"
printf 'Local directory   : %s\n' "${LOCAL_FIGURE_DIRECTORY}"

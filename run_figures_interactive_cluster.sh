#!/usr/bin/env bash
set -euo pipefail

script_directory="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cluster_host="${CLUSTER_HOST:-PSMN_sr650node230}"
remote_data_path="${DYNAMO_REMOTE_DATA_PATH:-}"

printf '\nDYNAMO — REMOTE INTERACTIVE COMPUTATION\n'
printf 'Cluster            : %s\n' "${cluster_host}"
printf 'Remote repository  : ~/Houches2026\n'
printf 'The following menu runs on the cluster and can access /Xnfs.\n\n'

if [[ -n "${remote_data_path}" ]]; then
    printf -v quoted_remote_data_path '%q' "${remote_data_path}"
    ssh -t "${cluster_host}" \
        "cd \"\$HOME/Houches2026\" && DYNAMO_COMPARISON_REPOSITORY=${quoted_remote_data_path} bash run_figures_interactive.sh"
else
    ssh -t "${cluster_host}" \
        'cd "$HOME/Houches2026" && bash run_figures_interactive.sh'
fi

printf '\nRemote computation completed successfully.\n'

download_after_run="${DOWNLOAD_AFTER_RUN:-true}"
case "${download_after_run}" in
    1|true|TRUE|yes|YES|y|Y)
        printf 'Downloading generated figures to the laptop...\n'
        CLUSTER_HOST="${cluster_host}" \
            bash "${script_directory}/download_figures.sh"
        ;;
    *)
        printf 'Figure download skipped (DOWNLOAD_AFTER_RUN=%s).\n' \
            "${download_after_run}"
        ;;
esac

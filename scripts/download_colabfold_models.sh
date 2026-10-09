#!/usr/bin/env bash
set -euo pipefail
# AlphaFold2 multimer-v3 parameters for ColabFold. This reproduces the layout of
# colabfold.download.download_alphafold_params (a params/ subdirectory plus the
# success marker colabfold_batch checks) using only curl and tar, so it works on
# GPU hosts without Docker, including the conda and micromamba profiles.
destination="${1:-models/colabfold}"
mkdir -p "$destination"
destination="$(cd "$destination" && pwd)"
params="$destination/params"
marker="$params/download_complexes_multimer_v3_finished.txt"
url='https://storage.googleapis.com/alphafold/alphafold_params_colab_2022-12-06.tar'

# ColabFold writes this marker empty, so test for existence rather than size.
if [[ -f "$marker" ]]; then
    echo "Using existing ColabFold parameters: $params"
    echo "Pass --colabfold_models $destination"
    exit 0
fi

mkdir -p "$params"
# Stream the archive straight into place; it is several GB and needs no second copy.
# The marker is written only after extraction succeeds, so an interrupted download
# is retried rather than silently reused.
# --no-same-owner: the archive records upstream uids, which root would try to restore
# and fail on many container filesystems.
curl --fail --location --retry 3 "$url" | tar -x --no-same-owner -C "$params"
test -s "$params/params_model_1_multimer_v3.npz"
: > "$marker"

echo "Downloaded AlphaFold2 multimer-v3 parameters to $params"
echo "Pass --colabfold_models $destination"

#!/usr/bin/env bash
set -euo pipefail
# Fetch the third-party model code BinderFlow drives. These checkouts are large and
# are not tracked in this repository; the revisions below match containers/*.Dockerfile
# so the conda profiles execute the same code as the container images.

RFDIFFUSION_REV="${RFDIFFUSION_REV:-86507b6538f51fce57b5a72477165f03999ed7ae}"
PROTEINMPNN_REV="${PROTEINMPNN_REV:-8907e6671bfbfc92303b5f79c4b5e6ce47cdef57}"

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
destination="${1:-$root}"
mkdir -p "$destination"
destination="$(cd "$destination" && pwd)"

fetch() {
    local name=$1 url=$2 rev=$3
    local path="$destination/$name"
    if [[ -e "$path/.git" ]] && [[ "$(git -C "$path" rev-parse HEAD 2>/dev/null || true)" == "$rev" ]]; then
        echo "Using existing $name at ${rev:0:7}"
        return 0
    fi
    if [[ -e "$path" ]]; then
        echo "$name already exists at $path but is not at ${rev:0:7}." >&2
        echo "Remove or move it and re-run to avoid discarding local changes." >&2
        return 1
    fi
    mkdir -p "$path"
    git -C "$path" init --quiet
    git -C "$path" remote add origin "$url"
    # Fetch only the pinned commit: full history of these repositories is not needed.
    git -C "$path" fetch --quiet --depth 1 origin "$rev"
    git -C "$path" checkout --quiet --detach FETCH_HEAD
    echo "Fetched $name at ${rev:0:7}"
}

fetch RFdiffusion https://github.com/RosettaCommons/RFdiffusion.git "$RFDIFFUSION_REV"
fetch ProteinMPNN https://github.com/dauparas/ProteinMPNN.git "$PROTEINMPNN_REV"

# Fail loudly now rather than inside a Nextflow task.
test -s "$destination/RFdiffusion/scripts/run_inference.py"
test -d "$destination/RFdiffusion/env/SE3Transformer/se3_transformer"
test -s "$destination/ProteinMPNN/protein_mpnn_run.py"
test -s "$destination/ProteinMPNN/vanilla_model_weights/v_48_020.pt"

cat <<EOF
Third-party code ready. ProteinMPNN ships its own weights; RFdiffusion weights come
from scripts/download_models.sh. Pass these to Nextflow:
  --rfdiffusion_script $destination/RFdiffusion/scripts/run_inference.py
  --proteinmpnn_script $destination/ProteinMPNN/protein_mpnn_run.py
EOF

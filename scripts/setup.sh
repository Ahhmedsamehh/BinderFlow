#!/usr/bin/env bash
set -euo pipefail
# One-shot setup for a fresh clone: fetches everything large that this repository
# deliberately does not track, then prints the command to run the workflow.
# Needs roughly 5 GB of disk and an internet connection. Safe to re-run: each step
# skips work that is already complete.
#
#   THIRD_PARTY_DIR  where RFdiffusion/ and ProteinMPNN/ are placed (default: repo root)
#   MODELS_DIR       RFdiffusion binder checkpoint   (default: <repo>/models/rfdiffusion)
#   COLABFOLD_DIR    AlphaFold2 multimer parameters (default: <repo>/models/colabfold)

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/.." && pwd)"
third_party="${THIRD_PARTY_DIR:-$root}"
models="${MODELS_DIR:-$root/models/rfdiffusion}"
colabfold="${COLABFOLD_DIR:-$root/models/colabfold}"

bash "$here/setup_third_party.sh" "$third_party"
bash "$here/download_models.sh" "$models"
bash "$here/download_colabfold_models.sh" "$colabfold"

third_party="$(cd "$third_party" && pwd)"
models="$(cd "$models" && pwd)"
colabfold="$(cd "$colabfold" && pwd)"

cat <<EOF

Setup complete. Run the workflow with:

  nextflow run main.nf -profile micromamba -params-file params.yaml \\
    --models $models \\
    --colabfold_models $colabfold \\
    --rfdiffusion_script $third_party/RFdiffusion/scripts/run_inference.py \\
    --proteinmpnn_script $third_party/ProteinMPNN/protein_mpnn_run.py

Use -profile conda or -profile mamba instead if that is the tool you have.
The first run solves the conda environments, which takes several minutes.
EOF

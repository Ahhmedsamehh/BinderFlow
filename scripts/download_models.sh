#!/usr/bin/env bash
set -euo pipefail
# BinderFlow's binder-only stage explicitly selects this upstream checkpoint.
destination="${1:-models}"
mkdir -p "$destination"
checkpoint="$destination/Complex_base_ckpt.pt"
if [[ -s "$checkpoint" ]]; then
    echo "Using existing checkpoint: $checkpoint"
    exit 0
fi
curl --fail --location --retry 3 \
    'https://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt' \
    --output "$checkpoint.part"
test -s "$checkpoint.part"
mv "$checkpoint.part" "$checkpoint"
echo "Downloaded $checkpoint. Each real run records its SHA-256 in run.json."

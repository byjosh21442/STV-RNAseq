#!/usr/bin/env bash
# Install the Nextflow launcher (official self-install) into the active stv-rnaseq env and pin its version.
# Usage: conda activate stv-rnaseq && bash scripts/bash/install_nextflow.sh
# Source: https://docs.seqera.io/nextflow/install and https://docs.seqera.io/nextflow/updating-nextflow
set -euo pipefail

NXF_PIN="26.10.3"            # latest stable release on 2026-10-25 (GitHub releases)
ENV_NAME="stv-rnaseq"

[ "${CONDA_DEFAULT_ENV:-}" = "$ENV_NAME" ] || { echo "ERROR: activate '$ENV_NAME' first" >&2; exit 1; }
java -version 2>&1 | head -n 1

if [ ! -x "$CONDA_PREFIX/bin/nextflow" ]; then
  tmp="$(mktemp -d)"
  ( cd "$tmp" && export NXF_VER="$NXF_PIN" && curl -fsSL https://get.nextflow.io | bash )
  chmod +x "$tmp/nextflow"
  # The env's bin/ is already on PATH whenever the env is active (docs: "a directory already in your PATH").
  mv "$tmp/nextflow" "$CONDA_PREFIX/bin/nextflow"
  rm -rf "$tmp"
fi

# Env-scoped variables (no change to ~/.bashrc or global conda config); nf-core recommends capping the JVM heap.
conda env config vars set -n "$ENV_NAME" NXF_VER="$NXF_PIN" NXF_OPTS="-Xms1g -Xmx4g" >/dev/null
echo "Re-activate the env to apply NXF_VER/NXF_OPTS, then run: nextflow -version"

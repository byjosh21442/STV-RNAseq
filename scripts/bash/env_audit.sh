#!/usr/bin/env bash
# Read-only environment audit. Prints system info and tool versions; installs/modifies nothing.
set -u

section() { printf '\n=== %s ===\n' "$1"; }

check_tool() {
  local name="$1"; shift
  local path
  path="$(command -v "$name" 2>/dev/null || true)"
  if [ -n "$path" ]; then
    printf '%-14s FOUND   %s | ' "$name" "$path"
    "$@" 2>&1 | head -n 1
  else
    printf '%-14s MISSING\n' "$name"
  fi
}

section "Date"
date -u +"%Y-%m-%dT%H:%M:%SZ"

section "Working directory"
pwd

section "OS"
grep -E '^(PRETTY_NAME|VERSION_ID|VERSION_CODENAME)=' /etc/os-release

section "Kernel / WSL"
uname -srm
cat /proc/version
echo "WSL_DISTRO_NAME=${WSL_DISTRO_NAME:-<unset>}"
command -v wslinfo >/dev/null 2>&1 && wslinfo --version 2>/dev/null | head -n 3

section "CPU"
lscpu | grep -E '^(Model name|Architecture|CPU\(s\)|Thread\(s\) per core|Core\(s\) per socket|Socket\(s\))'
echo "nproc: $(nproc)"

section "Memory"
free -h

section "Disk"
df -h "$PWD" "$HOME" / 2>/dev/null

section "Core tools"
check_tool git          git --version
check_tool conda        conda --version
check_tool mamba        mamba --version
check_tool micromamba   micromamba --version
check_tool java         java -version
check_tool nextflow     nextflow -version
check_tool R            R --version
check_tool Rscript      Rscript --version
check_tool python3      python3 --version
check_tool python       python --version
check_tool docker       docker --version
check_tool apptainer    apptainer --version
check_tool singularity  singularity --version

section "Bioinformatics / transfer tools"
check_tool curl         curl --version
check_tool wget         wget --version
check_tool pigz         pigz --version
check_tool prefetch     prefetch --version
check_tool fasterq-dump fasterq-dump --version
check_tool fastqc       fastqc --version
check_tool multiqc      multiqc --version
check_tool salmon       salmon --version
check_tool STAR         STAR --version
check_tool hisat2       hisat2 --version
check_tool featureCounts featureCounts -v

section "Conda installs not on PATH"
for d in "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/mambaforge" /opt/conda; do
  [ -x "$d/bin/conda" ] && echo "$d/bin/conda -> $("$d/bin/conda" --version 2>&1)"
done
echo "(end)"

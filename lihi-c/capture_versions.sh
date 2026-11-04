#!/usr/bin/env bash
# Run in the original hicup_clean environment from an authenticated Iris session.
# No Slurm job, common.sh, relative sourcing, or SSH is needed by this helper.
set -euo pipefail
out="${1:?Usage: bash capture_iris_versions.sh /absolute/new/version_output_directory}"
[[ "$out" = /* ]] || { echo "Use an absolute output directory" >&2; exit 1; }
mkdir -p "$out"
# Record executable locations and reported versions for this run.
{
  date -u
  command -v hicup || true
  hicup --version 2>&1 || true
  command -v juicer_tools || true
  command -v bowtie2 || true
  bowtie2 --version 2>&1 || true
  samtools --version 2>&1 || true
} > "$out/tool_versions.txt"
# Save the R environment and the HiC-DC+ functions used for feature construction/modeling.
Rscript - <<'RS' > "$out/R_versions_and_HiCDCPlus_implementation.txt" 2>&1
library(HiCDCPlus)
print(packageVersion("HiCDCPlus"))
print(packageDescription("HiCDCPlus"))
print(sessionInfo())
for (n in c("construct_features", "construct_features_chr", "expand_1D_features", "HiCDCPlus", "HiCDCPlus_chr")) {
 cat("\nFUNCTION:",n,"\n")
 print(get(n, envir=asNamespace("HiCDCPlus")))
}
RS
# Also save the environment package list when conda or micromamba is available.
if command -v conda >/dev/null 2>&1; then
  conda list --explicit > "$out/conda_explicit.txt"
  conda list --revisions > "$out/conda_revisions.txt"
elif command -v micromamba >/dev/null 2>&1; then
  micromamba list --explicit > "$out/micromamba_explicit.txt"
fi
echo "Captured CURRENT installation metadata in $out; compare against historical run dates."

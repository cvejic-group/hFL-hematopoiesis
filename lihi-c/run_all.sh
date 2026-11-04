#!/usr/bin/env bash
# Explicit absolute code path works even when Slurm spools this script.
#SBATCH --job-name=disomy_lihic
#SBATCH --partition=cpu
#SBATCH --time=3-00:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
# Stop if any step fails, rather than continuing with incomplete results.
set -euo pipefail
# Use the explicit code directory supplied at submission, not the Slurm spool path.
code="${1:?Usage: run_all.sh /absolute/code_directory /absolute/config.json /absolute/samples.tsv}"
config="${2:?Absolute config.json required}"
samples="${3:?Absolute samples.tsv required}"
for path in "$code" "$config" "$samples"; do
 [[ "$path" = /* ]] || { echo "All three paths must be absolute: $path" >&2; exit 1; }
done
# Read the output directory from the same configuration used by steps 01 and 02.
work=$(python3 - "$config" <<'PYTHON'
import json,sys
from pathlib import Path
print(Path(json.load(open(sys.argv[1]))['work_dir']).expanduser().resolve())
PYTHON
)
mkdir -p "$work"
# Save software versions before running the five analysis steps.
bash "$code/capture_versions.sh" "$work/software_versions"
python3 "$code/01_run_hicup.py" "$config" "$samples"
python3 "$code/02_make_contact_maps.py" "$config"
# Fit each biological replicate separately at 5, 10 and 25 kb.
for sample in 16216 17759; do
 Rscript "$code/03_call_interactions.R" "$work/maps/$sample/$sample.hic" "$work/calls/$sample" "$sample" "${SLURM_CPUS_PER_TASK:-4}"
done
# Summarize library quality, then compare the two 25-kb interaction tables.
Rscript "$code/04_library_qc.R" "$work/hicup" "$work/library_QC"
Rscript "$code/05_replicate_qc.R" \
 "$work/calls/16216/16216_hg38_25000_AAGCTT_combinedresult.txt.gz" \
 "$work/calls/17759/17759_hg38_25000_AAGCTT_combinedresult.txt.gz" "$work/replicate_QC"

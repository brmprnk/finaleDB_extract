#!/bin/bash
#SBATCH --job-name=finaledb_convert
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=32:00:00
#SBATCH --output=logs/%x_%j.out
#
# Convert the downloaded fragment files of one cohort to BAMs: sbatch slurm/convert.sbatch cris|jiang [IDS]
#   IDS defaults to lists/<cohort>_ids.txt; input <FRAGMENT_DIR>/<cohort>_frag/, output <PREPROCESSING_DIR>/<cohort>/
#   Records run 4 threads each, CPUs / 4 at a time (16 on 64 CPUs, about 15 GB memory and roughly 12 GB temporary disk each).
#   For cris, the per-record MAPQ filter of the original BAMs is taken from lists/cris_conversion_settings.csv; set
#   SETTINGS=none to convert every record with the default filter (MAPQ > 5).
# Submit from the repository root after loading slurm/site.env (see slurm/site.env.example).
cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
. slurm/_site.sh
COHORT="${1:?usage: sbatch slurm/convert.sbatch cris|jiang [IDS]}"
IDS="${2:-lists/${COHORT}_ids.txt}"
THREADS=4
PARALLEL=$(( ${SLURM_CPUS_PER_TASK:-4} / THREADS )); [ "${PARALLEL}" -ge 1 ] || PARALLEL=1

SETTINGS_ARG=()
if [ "${SETTINGS:-}" != none ] && [ -f "lists/${COHORT}_conversion_settings.csv" ]; then
	SETTINGS_ARG=(-s "lists/${COHORT}_conversion_settings.csv")
fi

${RUN} bash scripts/convert.sh -i "${IDS}" -f "${FRAGMENT_DIR}/${COHORT}_frag" -o "${PREPROCESSING_DIR}/${COHORT}" \
	-g "${GENOME}" "${SETTINGS_ARG[@]}" -j "${PARALLEL}" -t "${THREADS}" ${CONVERT_TMP:+--temp "${CONVERT_TMP}"} \
	--logs "logs/convert_${COHORT}"

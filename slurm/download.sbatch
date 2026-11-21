#!/bin/bash
#SBATCH --job-name=finaledb_download
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#
# Download the fragment files of one cohort: sbatch slurm/download.sbatch cris|jiang [IDS]
#   IDS defaults to lists/<cohort>_ids.txt; output <FRAGMENT_DIR>/<cohort>_frag/
# Submit from the repository root after loading slurm/site.env (see slurm/site.env.example).
cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
. slurm/_site.sh
COHORT="${1:?usage: sbatch slurm/download.sbatch cris|jiang [IDS]}"
IDS="${2:-lists/${COHORT}_ids.txt}"

${RUN} bash scripts/download.sh -i "${IDS}" -o "${FRAGMENT_DIR}/${COHORT}_frag" -j "${SLURM_CPUS_PER_TASK:-4}"

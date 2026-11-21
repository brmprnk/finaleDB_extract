#!/bin/bash
#SBATCH --job-name=finaledb_index
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.out
#
# Index the BAMs of one cohort that have no (or an outdated) index: sbatch slurm/index.sbatch cris|jiang
# Not needed after slurm/convert.sbatch, which indexes every BAM; this repairs a folder.
# Submit from the repository root after loading slurm/site.env (see slurm/site.env.example).
cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
. slurm/_site.sh
COHORT="${1:?usage: sbatch slurm/index.sbatch cris|jiang}"

${RUN} bash scripts/index.sh -d "${PREPROCESSING_DIR}/${COHORT}" -j "${SLURM_CPUS_PER_TASK:-4}"

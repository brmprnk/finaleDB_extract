#!/bin/bash
#SBATCH --job-name=finaledb_chr21_test
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#
# Build the chr21 test set (lists/chr21_test.csv): sbatch slurm/chr21_test.sbatch [OUTDIR]
#   OUTDIR defaults to <PREPROCESSING_DIR>/chr21_test; fragments are read from <FRAGMENT_DIR>/cris_frag/ when present,
#   otherwise only chr21 is fetched from FinaleDB with tabix.
# Submit from the repository root after loading slurm/site.env (see slurm/site.env.example).
cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
. slurm/_site.sh
OUTDIR="${1:-${PREPROCESSING_DIR}/chr21_test}"

${RUN} bash scripts/make_chr21_test.sh -g "${GENOME}" -o "${OUTDIR}" -f "${FRAGMENT_DIR}/cris_frag" -j 4 -t 2

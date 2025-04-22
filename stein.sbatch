#! /bin/sh
#SBATCH --partition=general --qos=medium
#SBATCH --cpus-per-task=64
#SBATCH --mem=240000
#SBATCH --time=31:59:00
#SBATCH --job-name=stein_p
#SBATCH --mail-user=ibpronk
#SBATCH --mail-type=END
#SBATCH --output=slurm_%j.out # Set name of output log. %j is the Slurm jobId
#SBATCH --error=slurm_%j.err # Set name of error log. %j is the Slurm jobId

ml use /opt/insy/modulefiles;

./stein_parallel.sh 64

#! /bin/sh
#SBATCH --partition=general --qos=medium
#SBATCH --cpus-per-task=2
#SBATCH --mem=25000
#SBATCH --time=31:59:00
#SBATCH --job-name=DL
#SBATCH --mail-user=ibpronk
#SBATCH --mail-type=END
#SBATCH --output=slurm_%j.out # Set name of output log. %j is the Slurm jobId
#SBATCH --error=slurm_%j.err # Set name of error log. %j is the Slurm jobId

ml use /opt/insy/modulefiles;

./download_ee.sh

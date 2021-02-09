#!/bin/bash
#SBATCH --qos=medium
#SBATCH --time=40:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

module load generic
module load singularity

srun singularity exec -u ~/r_latest.sif bash invoker.sh "$@"

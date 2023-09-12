#!/bin/bash
#SBATCH --time=0:20:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --array=1-50

module load julia
srun julia slurmjob_nmonitoring.jl $SLURM_ARRAY_TASK_ID

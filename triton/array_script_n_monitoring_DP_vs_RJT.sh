#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --array=1-50

module load julia
srun julia slurmjob_nmonitoring_CVaR.jl $SLURM_ARRAY_TASK_ID
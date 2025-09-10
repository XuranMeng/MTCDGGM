#!/bin/bash
#SBATCH --job-name=prepare
#SBATCH --output=logdir/%j-%x.log
#SBATCH --error=logdir/error_%j-%x.log
#SBATCH --time=80:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10GB
#SBATCH --partition=standard
#SBATCH --account=yili0
#SBATCH --mail-type=ALL

# Load MATLAB module (adjust according to your HPC environment)
module load matlab

# Run MATLAB script
n=800
matlab -nodisplay -nodesktop -r "n=$n; run('mainunknownprepare.m'); exit;"
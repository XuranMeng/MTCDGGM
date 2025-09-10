#!/bin/bash
#SBATCH --job-name=unbGamma
#SBATCH --error=logdir/error_%j-%x.log
#SBATCH --time=80:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=60GB
#SBATCH --partition=standard
#SBATCH --account=yili0
#SBATCH --mail-type=ALL

# Load MATLAB module (adjust according to your HPC environment)
module load matlab

# Run MATLAB script
n=200
matlab -nodisplay -nodesktop -r "n=$n; run('mainunknownstep1_bGamma.m'); exit;"
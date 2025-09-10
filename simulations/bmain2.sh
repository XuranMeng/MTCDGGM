#!/bin/bash
#SBATCH --job-name=main2
#SBATCH --output=logdir/%j-%x.log
#SBATCH --error=logdir/error_%j-%x.log
#SBATCH --time=40:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=11
#SBATCH --mem=60GB
#SBATCH --partition=standard
#SBATCH --account=yili0
#SBATCH --mail-type=ALL

# Load MATLAB module (adjust according to your HPC environment)
module load matlab

# Run MATLAB script
n=800
matlab -nodisplay -nodesktop -r "n=$n; run('mainstep2_optimizationcomparetpr.m'); exit;"
#!/bin/bash
#SBATCH --job-name=main1cv
#SBATCH --output=logdir/%j-%x.log
#SBATCH --error=logdir/error_%j-%x.log
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=9
#SBATCH --mem=180GB
#SBATCH --partition=standard
#SBATCH --account=yili0
#SBATCH --mail-type=ALL

# Load MATLAB module (adjust according to your HPC environment)
module load matlab

# Run MATLAB script
n=100
matlab -nodisplay -nodesktop -r "n=$n; run('mainstep1_regressioncv.m'); exit;"

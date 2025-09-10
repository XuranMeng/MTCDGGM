#!/bin/bash
#SBATCH --job-name=unknown3
#SBATCH --output=logdir/%j-%x.log
#SBATCH --error=logdir/error_%j-%x.log
#SBATCH --time=50:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=11
#SBATCH --mem=80GB
#SBATCH --partition=standard
#SBATCH --account=yili1
#SBATCH --mail-type=ALL

# Load MATLAB module (adjust according to your HPC environment)
module load matlab

# Run MATLAB script
n=800
matlab -nodisplay -nodesktop -r "n=$n; run('mainunknownstep3_optim.m'); exit;"
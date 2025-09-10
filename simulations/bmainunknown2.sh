#!/bin/bash
#SBATCH --job-name=unbGamma
#SBATCH --error=logdir/error_%j-%x.log
#SBATCH --time=80:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=13
#SBATCH --mem=180GB
#SBATCH --partition=standard
#SBATCH --account=yili1
#SBATCH --mail-type=ALL

# Load MATLAB module (adjust according to your HPC environment)
module load matlab

# Run MATLAB script
n=400
lbde=0.6
pool=12
matlab -nodisplay -nodesktop -r "n=$n;lbde=$lbde;pool=$pool; run('mainunknownstep2_regfixlbd.m'); exit;"
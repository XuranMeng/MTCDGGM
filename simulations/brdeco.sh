#!/bin/bash
#SBATCH --job-name=decorelation
#SBATCH --output=logdir/%j-%x.log
#SBATCH --error=logdir/error_%j-%x.log
#SBATCH --time=40:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80GB
#SBATCH --partition=standard
#SBATCH --account=yili0
#SBATCH --mail-type=ALL

# Load R module (adjust based on your HPC environment)
module load R

# Set parameter
n=200

# Run R script with n as argument
Rscript comparemethod/decorelation/test.R $n
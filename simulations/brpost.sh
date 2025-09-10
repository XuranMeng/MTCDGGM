#!/bin/bash
#SBATCH --job-name=postcompare
#SBATCH --output=logdir/%j-%x.log
#SBATCH --error=logdir/error_%j-%x.log
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20GB
#SBATCH --partition=standard
#SBATCH --account=yili0
#SBATCH --mail-type=ALL

# Load R module (adjust based on your HPC environment)
module load R

# Set parameter
n=800

# Run R script with n as argument
Rscript comparemethod/postselection/trial.R $n
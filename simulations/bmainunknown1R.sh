#!/bin/bash
#SBATCH --job-name=unbGamma
#SBATCH --error=logdir/error_%j-%x.log
#SBATCH --output=logdir/output_%j-%x.log
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60GB
#SBATCH --partition=standard
#SBATCH --account=yili1
#SBATCH --mail-type=ALL

module load R

n=400

Rscript mainunknownstep1_bGamma.R $n
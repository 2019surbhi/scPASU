#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=saf
#SBATCH -n 20
#SBATCH --mem-per-cpu=3900

#SBATCH -o saf%J.out
#SBATCH -e saf_%J.err


### Load modules ###

module load R/4.1.0

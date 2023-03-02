#!/bin/bash

#SBATCH -t 100:00:00
#SBATCH --job-name=MRC_paleo
#SBATCH -p normal
#SBATCH --export=ALL
#SBATCH --nodes=1               
#SBATCH --output=MRC_paleo.out
#SBATCH --output=MRC_paleo.err
#SBATCH --ntasks-per-node=16


module load R

Rscript config.simulations.basins.april.paleo.hopper.R


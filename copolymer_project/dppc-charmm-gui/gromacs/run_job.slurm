#!/bin/bash
#SBATCH -N 1                    # Run on a single CPU
#SBATCH --time=02:00:00               # Time limit hrs:min:sec
#SBATCH --partion=batch
#SBATCH --job-name=SMA
#SBATCH --output=gmx_SMA_%j.log
#SBATCH --error=gmx_SMA_%j.err

module load gromacs-2022.3
./minim_equil.csh

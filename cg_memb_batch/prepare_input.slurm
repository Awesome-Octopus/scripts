#!/bin/bash
#SBATCH --job-name=gmx    
#SBATCH --mail-type=END          
#SBATCH --mail-user=morri361@miamioh.edu 	
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=32gb                     # Job memory request
#SBATCH --time=20:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log


module load gromacs-2022.3
./prepare_input.sh

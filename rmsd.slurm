#!/bin/bash
#SBATCH --job-name=gmx    
#SBATCH --mail-type=END     
#SBATCH --mail-user=morri361@miamioh.edu 	
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=4gb                     # Job memory request
#SBATCH --time=01:00:00               # Time limit hrs:min:sec
#SBATCH --output=gmx_mpi_%j.log   # Standard output and error log


module load gromacs-2022.3

#echo 6| gmx_mpi trjconv -f recentered_traj.xtc -s step7_production.tpr -o just_protein.xtc
echo 0 0| gmx_mpi rms -f just_protein.xtc -s just_protein.pdb -o rmsd.xvg 2> err.log

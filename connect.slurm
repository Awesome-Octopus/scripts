#!/bin/bash
#SBATCH --job-name=gmx    
#SBATCH --mail-type=END     
#SBATCH --mail-user=morri361@miamioh.edu 	
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=4gb                     # Job memory request
#SBATCH --time=01:00:00               # Time limit hrs:min:sec
#SBATCH --output=gmx_mpi_%j.log   # Standard output and error log


module load gromacs-2022.3


echo 6 0 | gmx_mpi trjconv -f step7_production.gro -s step7_production.tpr -o recentered_traj.gro -pbc mol -center 2>error.log
echo 0 | gmx_mpi trjconv -f recentered_traj.gro -s step7_production.tpr -conect -o connected_traj.pdb
echo 6 0 | gmx_mpi trjconv -f step7_production.xtc -s step7_production.tpr -o recentered_traj.xtc -pbc mol -center 2> error.log
sed -i '/ENDMDL/d' connected_traj.pdb


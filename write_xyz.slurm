#!/bin/bash
#SBATCH --job-name=coords    
#SBATCH --mail-type=END,FAIL          
#SBATCH --mail-user=morri361@miamioh.edu 	
#SBATCH --ntasks=1                    # Run on 1 CPU
#SBATCH --mem=64gb                     # Job memory request
#SBATCH --time=02:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log


module load vmd

echo -ne "source ~/scripts/write_xyz.tcl\n quit\n" | vmd -dispdev text recentered_traj.xtc connected_traj.pdb

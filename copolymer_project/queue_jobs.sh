#!/bin/bash

for file in $(ls input_structs/*.pdb); do

	
	bn=$(basename "${file}" ".pdb")
	num="${bn//[!0-9]/}"
	#		                                        working directory                              path to dssp               path to params     
	#		                                     	|                                                     |                   |         
	sbatch -o "${num}_sbatch.log" -J "${num}" run_job.slurm $(pwd)"/replicates/""${bn}" ${file} $(pwd)"/""scripts/dssp-2.2.1/mkdssp" $(pwd)"/params/"  	

done

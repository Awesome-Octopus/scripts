#!/bin/bash

#SBATCH --job-name=mono    
#SBATCH --mail-user=morri361@miamioh.edu
#SBATCH --mail-type=END 	
#SBATCH --nodes=1                    # Run on a single CPU
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -p batch
#SBATCH --time=100:00:00               # Time limit hrs:min:sec
#SBATCH --output=parallel_batch.log   # Standard output and error


./batch-run.sh
#./batch-run_12-20.sh
#./batch-run_23-30.sh
#./batch-run_33-40.sh
#./batch-run_42-50.sh
#./batch-run_53-60.sh
#./batch-run_62-70.sh
#./batch-run_73-80.sh
#./batch-run_83-90.sh
#./batch-run_93-100.sh

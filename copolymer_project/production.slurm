#!/bin/bash
#SBATCH --job-name=min-eq
#SBATCH --mail-type=END
#SBATCH --mail-user=morri361@miamioh.edu
#SBATCH --ntasks=2                    # Run on a single CPU
#SBATCH --mem=16gb                     # Job memory request
#SBATCH --time=99:00:00               # Time limit hrs:min:sec
#SBATCH --output=gmx_mpi_%j.log   # Standard output and error log

wrkdir=${1}"/"
bn=${wrkdir}$(basename "${2}" ".pdb")
dssp=${3}
to_params=${4}
	
	
mkdir -p "${wrkdir}"

cd "${wrkdir}"

        # Production
       gmx_mpi grompp -f "${to_params}step7_production.mdp" -o step7_production.tpr -c step6.6_equilibration.gro -p system.top -n index.ndx
       gmx_mpi mdrun -deffnm step7_production

cd ../..

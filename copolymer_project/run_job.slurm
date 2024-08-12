#!/bin/bash
#SBATCH -N 1                    # Run on a single CPU
#SBATCH -n 24
#SBATCH --time=12:00:00               # Time limit hrs:min:sec

wrkdir=${1}"/"
bn=${wrkdir}$(basename "${2}" ".pdb")
dssp=${3}
to_params=${4}
	
	
mkdir -p "${wrkdir}"

cp ${2} ${wrkdir}$(basename ${2})
python2 martinize.py -f "${bn}.pdb" -x "${bn}_CG.gro" -o "${bn}_CG.top" -p backbone \
-dssp "${dssp}" -ff martini22p > "${wrkdir}martinize.log" 2>&1

	# see this page for a detailed explanation of insane.py: 
	# https://github.com/Tsjerk/Insane/blob/master/README.md
	 
python2 insane.py -f "${bn}_CG.gro" -o "${bn}_CG_memb.gro" -p "${bn}_CG_memb.top" -d 10 -sol PW \
-l POPC:3 -l POPG:1 -charge auto -salt 0.15 > "${wrkdir}insane.log" 2>&1

	# remove all the old include statements. they are incorrectly formatted.
grep -vE "^#" "${bn}_CG_memb.top" > "${bn}""_temp.top"

	
	# put the base filename of all the itp files created in martinize2 in a file
	# NOTE if you change to martinize2 it autonames them molecule_x. change the search pattern here if you use that

	# ~~~~~~~~~~ EDIT ME ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# !!! There needs to be one line of "Protein_A      1" in your .top file for each (homo)monomer in your protein.
	# for now you have to do this manually. IF THE NUMBER OF SUBUNITS DOESNT MATCH YOU WILL GET AN ERROR
	# remove all include statements and change "Protein" to "Protein_A" in the itp file
	
	grep -vE "^#" "${bn}""_temp.top" | sed -E 's/Protein\s+[0-9]/Protein_A        5/' > "${wrkdir}system.top"
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	# fix sodium and cloride which are messed up for some reason
sed -i 's/NA+/NA /g' "${wrkdir}system.top"
sed -i 's/NA+/NA /g' "${bn}_CG_memb.gro"
sed -i 's/CL-/CL /g' "${wrkdir}system.top"
sed -i 's/CL-/CL /g' "${bn}_CG_memb.gro"
	
mv Protein_A.itp "${wrkdir}Protein_A.itp"

	# add back the include statements in the correct order with Protein_A.itp at the bottom of the header
for file in $(ls ${wrkdir}/*Protein*.itp); do
	echo "#include \"${file}\"" | cat - "${wrkdir}system.top" > "${bn}_temp.top" && mv "${bn}_temp.top" "${wrkdir}system.top"
done

	# !!!! for use with martini 2
	# then add martini top files, with the general forcefield martini_v2.2P.itp at the very top
echo "#include \""${to_params}"martini_v2.0_ions.itp\"" | cat - "${wrkdir}system.top" > "${bn}""_tempFile" && mv "${bn}""_tempFile" "${wrkdir}system.top"
echo "#include \""${to_params}"martini_v2.0_lipids_all_201506.itp\"" | cat - "${wrkdir}system.top" > "${bn}""_tempFile" && mv ${bn}_tempFile "${wrkdir}system.top"
echo "#include \""${to_params}"martini_v2.2P.itp\"" | cat - "${wrkdir}system.top" > "${bn}""_tempFile" && mv "${bn}""_tempFile" "${wrkdir}system.top"


module load gromacs-2022.3

	# !!!!!!!!!!!!!! NOTE: for this part you have to go through manually once and see what numbers each group has
	# we want to select POPC and POPG as a group and name it membrane, rename Protein as protein 
	# and everything else is grouped together and named solute.

gmx_mpi make_ndx -f "${bn}_CG_memb.gro" -o "${wrkdir}index.ndx" << EOF
name 1 protein
13 | 14
name 19 membrane
!1 & !19
name 20 solute
q
EOF

cd ${wrkdir}
cp ../../*.ssd .
export GMX_MAXCONSTRWARN=-1

	# step6.0 - soft-core minimization
gmx_mpi grompp -f ${to_params}step6.0_minimization.mdp -o step6.0_minimization.tpr \
-c "${bn}_CG_memb.gro" -r "${bn}_CG_memb.gro" -p system.top -n index.ndx \
-maxwarn -1 
	
gmx_mpi mdrun -deffnm step6.0_minimization

        # step6.1
gmx_mpi grompp -f "${to_params}step6.1_minimization.mdp" -o step6.1_minimization.tpr \
-c step6.0_minimization.gro -r ${bn}_CG_memb.gro -p system.top -n index.ndx -maxwarn -1 

gmx_mpi mdrun -deffnm step6.1_minimization 

unset GMX_MAXCONSTRWARN

        # Equilibration
cnt=2
cntmax=6
while [ ${cnt} -le ${cntmax} ]; do
    pcnt=$((cnt - 1))
    if [ $cnt -eq 2 ]; then
	gmx_mpi grompp -f "${to_params}step6.${cnt}_equilibration.mdp" -o "step6.${cnt}_equilibration.tpr" \
	-c step6.${pcnt}_minimization.gro -r ${bn}_CG_memb.gro  -p system.top -n index.ndx -maxwarn -1
    else
       gmx_mpi grompp -f "${to_params}step6.${cnt}_equilibration.mdp" -o "step6.${cnt}_equilibration.tpr" \
       -c step6.${pcnt}_equilibration.gro -r ${bn}_CG_memb.gro -p system.top -n index.ndx -maxwarn -1
    fi
    gmx_mpi mdrun -deffnm step6.${cnt}_equilibration
    cnt=$((cnt +1))
done

#### THIS DELETES INFORMATION THAT MIGHT BE USEFUL, BUT GREATLY REDUCES footprint
rm *.trr
rm *.xtc

        # Production
       gmx_mpi grompp -f "${to_params}step7_production.mdp" -o step7_production.tpr -c step6.6_equilibration.gro -p system.top -n index.ndx
       gmx_mpi mdrun -deffnm step7_production -cpo

cd ../..

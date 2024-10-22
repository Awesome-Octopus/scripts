#!/bin/csh
#
# Generated by CHARMM-GUI (http://www.charmm-gui.org)
#
# 1) Use Gromacs 5.1 or newer to run these simulations
#

# Minimization
setenv GMX_MAXCONSTRWARN -1
# step6.0 - soft-core minimization
# If you encountered "There are 1 perturbed non-bonded pair interaction ......" error message, 
# please modify rvdw and rcoulomb values from 1.1 to 2.0 in the step6.0_minimization.mdp file
#gmx grompp -f step6.0_minimization.mdp -o step6.0_minimization.tpr -c combined_input.pdb -r combined_input.pdb -p system.top -n combined_ndx.ndx -maxwarn 1
#gmx mdrun -deffnm step6.0_minimization

# step6.1
gmx_mpi grompp -f step6.1_minimization.mdp -o step6.1_minimization.tpr -c step6.0_minimization.gro -r combined_input.pdb -p system.top -n combined_ndx.ndx -maxwarn 1
gmx_mpi mdrun -deffnm step6.1_minimization
unsetenv GMX_MAXCONSTRWARN

# Equilibration
set cnt    = 2
set cntmax = 6

while ( ${cnt} <= ${cntmax} )
    @ pcnt = ${cnt} - 1
    if ($cnt == 2) then
        gmx_mpi grompp -f step6.${cnt}_equilibration.mdp -o step6.${cnt}_equilibration.tpr -c step6.${pcnt}_minimization.gro -r combined_input.pdb -p system.top -n combined_ndx.ndx
    else
        gmx_mpi grompp -f step6.${cnt}_equilibration.mdp -o step6.${cnt}_equilibration.tpr -c step6.${pcnt}_equilibration.gro -r combined_input.pdb -p system.top -n combined_ndx.ndx
    endif
    gmx_mpi mdrun -deffnm step6.${cnt}_equilibration
    @ cnt += 1
end

# Production
#gmx grompp -f step7_production.mdp -o step7_production.tpr -c step6.6_equilibration.gro -p system.top -n combined_ndx.ndx
#gmx mdrun -deffnm step7_production

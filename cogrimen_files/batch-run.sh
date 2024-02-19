#!/bin/bash

# a master shell script for prepping input for each file in a numbered 
# series and queueing up jobs for them
# with coarse grained implicit solvent membrane system in NAMD

# if you want to change something like temp, step time, etc
# and have it be universal accross replicates, then change
# the template file (template-input.inp) in the base directory
# since that is where all the copies are made from

# this script is intended to be run in the base ~/md/cg-implicit-s2e folder

module load vmd
module load openmpi-4.1.2-gnu-6.3.0-slurm
pathtonamd="/home/morri361/cogrimen_namd/Linux-x86_64-g++"


#ln -s $pathtonamd/charmrun charmrun
#ln -s $pathtonamd/namd2 cogrimen
# edit the input parameter file to be specific to this structure
for ((i=72; i<=72; i++)); do
     
    base_name=randomized-S2E-oriented-$i
    
    mkdir replicates/$base_name/
    mv replicates/$base_name.pdb replicates/$base_name/
    cp template-input.inp replicates/$base_name/$base_name.inp
    cd replicates/$base_name

    # if you have a different structure for your directory, you need to change
    # that in the following lines. this script assumes that your martini-protein.top, 
    # your dssp executable, and all tcl scripts are in the folder 2 levels up (../../)
    
    # uncomment the pipe to grep for debugging if necessary
    
    input=$base_name\.inp
    path=$(pwd)
    sed -i "s/set outputname.*/set outputname cogrimen-s2e-$i/" $input #| grep "set outputname cogrimen-s2e-$i"
  
    sed -i "s/structure.*/structure     S2E-processed-cg-input-$i\.psf/" $input #| grep "structure     S2E-processed-cg-input-$i\.psf" 
   
    sed -i "s/coordinates.*/coordinates   S2E-processed-cg-input-$i\.pdb/" $input #| grep "coordinates   S2E-processed-cg-input-$i\.pdb" 
        
    sed -i "s/martini-params/\.\.\/\.\.\/martini-params/" $input #| grep "\.\.\/\.\.\/martini-params"
    sed -i -E "s/CGIMSolventFile\s+\S+/CGIMSolventFile      \.\.\/\.\.\/solvpar_martini_COGRIMEN\.inp/" $input
    cp ../../00-make-AA-psf.tcl .
    sed -E -i "s/mol load pdb \S+\.pdb/mol load pdb $base_name\.pdb/" 00-make-AA-psf.tcl #| grep "mol load pdb $base_name\.pdb"
    sed -E -i "s/topology\s+(\S+)/topology \.\.\/\.\.\/\1/" 00-make-AA-psf.tcl #| grep "topology \.\.\/\.\.\/"
    
    cp ../../01-coarse-grain.tcl .
    sed -E -i "s/read_db\s+(\S+)/read_db \.\.\/\.\.\/\1/" 01-coarse-grain.tcl #| grep "read_db \.\.\/\.\.\/"
    cp ../../02a-make-initial-CG-psf.tcl .
    sed -E -i "s/topology\s+(\S+)/topology \.\.\/\.\.\/\1/" 02a-make-initial-CG-psf.tcl # | grep "topology \.\.\/\.\.\/"
    cp ../../02b-correct-CG-psf.tcl .
    sed -E -i "s/fix_martini_psf\.tcl/\.\.\/\.\.\/fix_martini_psf\.tcl/" 02b-correct-CG-psf.tcl | grep "\.\.\/\.\.\/fix_martini_psf\.tcl"
    sed -E -i "s/\s([^\s\/]+)\.top/\.\.\/\.\.\/\1\.top/" 02b-correct-CG-psf.tcl #| grep "\.\.\/\.\.\/\S+\.top"
    sed -E -i "s/cg-\S+ ([0-9]+) \.\S*/S2E-processed-cg-input-$i \1 \.\.\/\.\.\//" 02b-correct-CG-psf.tcl
    echo source 00-make-AA-psf.tcl | vmd -dispdev text > vmd.log 
    echo source 01-coarse-grain.tcl | vmd -dispdev text >> vmd.log
    echo source 02a-make-initial-CG-psf.tcl | vmd -dispdev text >> vmd.log
    echo source 02b-correct-CG-psf.tcl | vmd -dispdev text >> vmd.log
    
    # remove all instances of "P4 " from the psf file since it is not parametrized in our
    # .par files replace it with something similar
    sed -i 's/P4 /N0H/' "S2E-processed-cg-input-$i.psf"
    charmrun namd2 $input > "s2e-cogrimen-$i.log"
    cd ../../
    
done



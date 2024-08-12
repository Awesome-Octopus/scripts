#!/bin/bash

# arguments are to be passed to this script verbatim as they would be passed to the 
# md_helix_orientations script except without the filename for the output


conda activate md

cd replicates
for dir in */; do
	if [ -d "$dir" ]; then 
		cd "$dir"
		num=$(echo "$dir" | sed -En 's/.+[-|_]([0-9]+)\/$/\1/p') #extract the number from the directory
		if [ -f "step7_production.xtc" ]; then
			python3 ~/scripts/md_helix_orientations.py "$@" >> ../../averaged_angles.dat  
		fi
		cd ../
	fi
done
cd ../



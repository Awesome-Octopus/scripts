#!/bin/bash

#you need to change 1rep_memb.pdb to your concatenated pdb file
#you also need to include the path if this is not in your folder
#for this script to produce usable output redirect to a text file 
# ie: ./renumber.sh > new.pdb

echo $(sed -n '1p' 1rep_memb.pdb)
# for every line in your pdb file excluding the first 
for ((c=2; c < $(cat 1rep_memb.pdb |wc -l); c++))
do
# replace the atom number (the first instance of 1 or more digits in the line)
# with the line number
 let g=$c-1
 atomnum=$(printf "%7d" $g)
 sed -r -n $c"s/\s+[0-9]+/$atomnum/p" "1rep_memb.pdb"
done
echo END

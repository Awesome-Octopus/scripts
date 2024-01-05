#!/bin/bash

#######################################################################
# A unified script that should take output trajectories from the
# membrane production and the production run of the peptide in solution
# and combine them then modify all the relevent files to produce
# proper input for the combined run
#######################################################################



###############################################################
# you need to change minput to be the pdb of the last frame of
# membrane only production run as produced by VMD
# change prinput to be the last frame of your protein only
# production run saved with only (an)ions and protein
# change n_counterions to the number of counterions needed
# to neutralize your peptide charge
# if it is an anionic peptide you will need to change CL to NA
# later in the script
###############################################################

###############################################################
#  EDIT THESE IF NEEDED
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
###############################################################
f2="combined_input.pdb"
minput="step5_charmm2gmx.pdb"
prinput="protein_input.pdb"
out="updated_indeces.ndx"
let n_counterions=7

###############################################################

# remove END and TER lines from membrane pdb
sed -E '/^CRYST1\s*|END|TER/d' $minput > membrane.pdb

################################################################################
# if they aren't already there, add the MEMB IONS and WAT group tags to each
# atom in membrane pdb and remove X chain ID
# if you use other lipids than POPC and POPG, or other ions than NA and CL
# you will need to modify the parts pointed to by the ^ in the comment
################################################################################

sed -i -E 's/X/ /' membrane.pdb


# reformat all waters and move them to separate pdb
sed -n -E 's/(.*(WM|WP|W).*)($)/\1 WAT/p' membrane.pdb > temp_water.pdb
sed -i -E '/.*(WM|WP|W).*/d' membrane.pdb

# add MEMB group tag to the end of the line
sed -i -E 's/(.*(POPC|POPG|OANT|POPE|POPS).*)($)/\1 MEMB/' membrane.pdb
#                ^     ^

# separate out all ions, reformat them, and move them to a new pdb
sed -n -E 's/(.*(NA|CL).*)($)/\1 IONS/p' membrane.pdb > temp_ions.pdb
#                ^  ^

# remove the ions in the membrane pdb
sed -i -E '/.*(NA|CL).*/d' membrane.pdb
#              ^  ^

sed -i -E 's/(.*NA\s+)(ION)(\s+.*)/\1NA\3/' temp_ions.pdb
#               ^                    ^   
sed -i -E 's/(.*CL\s+)(ION)(\s+.*)/\1CL\3/' temp_ions.pdb
#               ^                    ^


# count the number of lines with MEMB in them
let nmemb=$(grep -c "MEMB" membrane.pdb)

#count how many atoms from the membrane pdb are ions
let n_membrane_ions=$(grep -c "IONS" temp_ions.pdb)

#the number of waters
let n_water=$(cat temp_water.pdb |wc -l)

#count how many are ions or water (will be the group solute in index)
let nsolute=$(expr $n_water + $n_membrane_ions + $n_counterions)

# extract n chlorides from protein pdb and properly format the lines
# append them to ions.pdb and delete the chain ID
# change CL below to your desired counterion if needed
grep ' CL ' $prinput | head -n $n_counterions | sed -E 's/ION/CL/' | sed -E 's/X/ /' | sed -E 's/$/IONS/' >> temp_ions.pdb

# in the protein input pdb get every line that doesn't include the header or NA or CL 
# then for every line that isn't TER or END correctly format them with PROA at the end 
# of the line


#now remove any lines that aren't protein atoms and reformat the remainder, removing chain ID
grep -v -E "TER|END|CL|NA|CRYST1" $prinput | sed -E 's/X/ /' | sed -E 's/$/PROA/' > protein.pdb  
#                   ^  ^

let nprot=$(cat protein.pdb | wc -l)

###################################################################################
# renumber the residue number on ions 
# and on water as each should have its own numbering scheme
###################################################################################


# for ions
gawk '{$5 = NR}{print $0}' temp_ions.pdb > ions.pdb


# for water since it has 3 atoms, if you are on a line 
# that divides evenly into 3 the next line will be a new molecule
gawk 'BEGIN {n = 1}
{ if(NR % 3 == 0 )
{
	{$5 = n}
	{n++}
}
else
{
	{$5 = n}
} 
{print $0}
}' temp_water.pdb > water.pdb

# now glue the individual pdbs together
cat membrane.pdb water.pdb ions.pdb protein.pdb > temp_combined.pdb


#############################
# renumber the atoms        #
#############################

# In each line the second field (atomnumber) is replaced with
# the line number.
# It's then reprinted to exactly match the format of step5_charmm2gmx.pdb

#        column spacing guide: if you need it, uncomment it
#echo 0........1.........2.........3.........4.........5.........6.........7.........8
gawk '{$2 = NR} {printf "%-4s%7d  %-3s %-4s%6d %10.3f %7.3f %7.3f %5.2f %5.2f      %-4s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' temp_combined.pdb > $f2


#add back the header at the first line of the new file
sed -i "1s/^/$(grep 'CRYST' $minput)\n/" $f2


#now add in the TER and END lines with approprate indexes for TER
tail -n 1 $f2 | gawk '{$2 = $2 + 1}{printf "%3s %7d %8s %7d\n", "TER", $2, $4, $5}' >> $f2
echo END >> $f2
 

###############################################################
# update the index file for the combined pdb you just created #
# with the newly assigned atom number as indeces              #
###############################################################

echo "Making the updated index file. Be patient. This can take a long time"
# make [ MEMBRANE ] index group
printf "[ MEMBRANE ]\n" >> $out
for (( i=1; i < $(expr $nmemb + 1); i++)) 
# everything with atom number less 
# than the number of memb atoms + 1
do
    if [ $(expr $i % 15) -eq 0 ] || [ $i -eq $nmemb ]
    then
        printf "%d\n" $i >> $out
    else
        printf "%d " $i >> $out
    fi
done

# make [ SOLUTE ] index group
let z=$(expr $nmemb + $nsolute + 1)
printf "[ SOLUTE ]\n" >> $out
for (( i= $(expr $nmemb + 1) ; $i < $z; i++ ))
# everything with atom number less than
# membrane atoms + number of water + number of ion +1
do
    if [ $(expr $(expr $i - $nmemb) % 15) -eq 0 ] || [ $i -eq $(expr $z - 1) ]
    then
        printf "%d\n" $i >> $out
    else
        printf "%d " $i >> $out
    fi
done

# make [ PROTEIN ] index group
let x=$(expr $z + $(grep -c -E 'PROA' $f2))  
printf "[ PROTEIN ]\n" >> $out
for (( i=$z ; $i < $x; i++))
do
    if [ $(expr $(expr $i - $z + 1) % 15) -eq 0 ]
    then
        printf "%d\n" $i >> $out
    else
        printf "%d " $i >> $out
    fi
done
printf "\n" >> $out

#########   clean up temp files   ####
rm temp_combined.pdb
rm protein.pdb
rm membrane.pdb
rm ions.pdb
rm water.pdb
rm temp_water.pdb
rm temp_ions.pdb
######################################

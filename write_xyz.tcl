
# made 1/22/23 by Andrew Morris
# this script will write the xyz coordinates of given beads for each frame
# use in conjunction with matlab script to get distances


# set res to be the list of residues you want distances between

# set res { 34 50 } <------ deprecated

# for VSLD resid = the residue number - 417
# since vmd views the first residue as 1 rather than
# the full sequence numbering change seq offset to 0 if the
# first residue of your VMD sequence is the first residue
# in your numbering convention

set seqoffset 0

# note that if the unmutated residue is a glycine
# there is no SC1 and this won't work
# in which case you can still get coordinates for the
# BB bead
# for gromacs all atom models the beta carbon has the name CB
# the alpha carbon is CA 
# and the oxygen where the spin is is ON
 
# which portion of your frames do you want to average over
# to produce your distribution

set startframe 100
set endframe 4700

# the name of the text file this script will produce
set outfile [open label_coordinates.dat w]

# making the header row ###########################################

############## write header for backbones first####################
puts $outfile "Spatial coordinates for [pwd]/connected.pdb"
set sel [atomselect top "name BB or name CA" frame first]
set ndx [$sel list]

puts -nonewline $outfile "frame:\t"

set sel [atomselect top "name SC1 or name CB" frame first]
append ndx [$sel list]

foreach num $ndx {
	set v [atomselect top "index $num" frame first]
	set id [$v get resname]
	append id "_[$v get name]"
	append id "_[expr [$v get resid] + $seqoffset]"
	append id "_[$v get chain]"
	puts -nonewline $outfile "${id}_X\t${id}_Y\t${id}_Z\t"
}


puts $outfile ""
######################################################

# for each selected bead in each frame output its X Y and Z
for {set i $startframe} {$i < [expr $endframe + 1]} {incr i} {
	

	# print frame number
	puts -nonewline $outfile "$i\t"


	####  N O T E : 
	####  gromacs has extremely high precision, but it also leads to large files, and difficulty reading 
	####  printed output will be truncated to 3 decimal places
	####  there sould be no reason to need more precision for our purposes, but
	####  ---- IF YOU NEED MORE PRECISION ---
	####  delete the "format %.2f" in the lines below
 	
	foreach index $ndx {
		set sel [atomselect top "index $index" frame $i]
		 
		set t1 [format %.2f [$sel get x]]
		set t2 [format %.2f [$sel get y]]
		set t3 [format %.2f [$sel get z]]
	
		puts -nonewline $outfile "$t1\t$t2\t$t3\t"
	}
	puts $outfile ""
}

close $outfile

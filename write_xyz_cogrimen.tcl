
# made 1/05/24 by Andrew Morris
# derived from write_xyz.tcl
# intended for use with cogrimen NAMD
# this script will write the xyz coordinates of given beads for each frame
# This assumes that the only beads are protein, otherwise you will
# need to modify the atomselect on line 69


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
# if nothing is given average over all frames
if {$argc == 0} {
    set startframe 0
    set endframe [molinfo top get last]
    # if one arg is given that is the start frame number
} elseif {$argc == 1} {
    set startframe [lindex $argv 0]
    set endframe [molinfo top get last]
	# if its 2 is specifies the start and the end
} elseif {$argc == 2} {
    set startframe [lindex $argv 0]
    set endframe [lindex $argv 1]
}
# the name of the text file this script will produce
# the file name needs to be passed to the script as the 3rd argument
# otherwise take the name of the traj file and use it as a base name
if {$argc == 3} {
	set outfile [open [lindex $argv 2] w]
} else {
	set infile_name [molinfo top get name]
	set length [string length $infile_name]
	puts $infile_name
	set a [expr $length - 5]
	puts $a
	set bn [string range $infile_name 0 $a]
	puts $bn
	set fn [pwd]/$bn.dat
	puts $fn
	set outfile [open $fn w+]
	
}

# making the header row ###########################################

############## write header for backbones first####################
set name [molinfo top get name]
puts $outfile "Spatial coordinates for [pwd]/$name"
set sel [atomselect top "all" frame first]
set ndx [$sel list]

puts -nonewline $outfile "frame:\t"

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

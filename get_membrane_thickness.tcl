#######################################################################
# change the following to change how many residues there are in your protein
# and which frame number past which to average
#######################################################################
set startframe 2

set outfile [open avg_po4_dist.dat w]


# you need to write your own function just to take the mean of an array
proc average {list} {
    expr {[tcl::mathop::+ {*}$list 0.0] / max(1, [llength $list])}
}

# make a header for variable columns
puts -nonewline $outfile "frame"

puts $outfile "\tdiff_of_po4_avg"

	# get the number of frames
set nf [molinfo top get numframes]

	# loop over all frames from startframe to end
	
	
for { set i $startframe } { $i < $nf } { incr i } {

	
	set phos [atomselect top "name PO4" frame $i]
	puts -nonewline $outfile "$i"
	set phos_z_pos [$phos get z]

	# mean z of all phosphates
	set zavgpo4 [average $phos_z_pos]
	# make separate lists for upper and lower leaflets z positions
	set z_upper {}
	set z_lower {}
	
	# sort each z value of a phosphate into the upper or lower leaflet pile
	# based on whether it is above or below the mean, then subtract the average
	# po4 z position to give height reletive the center of the bilayer at 0
	foreach z $phos_z_pos {
		if {$z > $zavgpo4} {
	
			lappend z_upper [expr $z - $zavgpo4]
		} else {
			lappend z_lower [expr $z - $zavgpo4]
		}
	}
	
	# take their averages separately
	
	set avg_z_upper [average $z_upper]
	set avg_z_lower [average $z_lower]
	set diff_z [expr $avg_z_upper - $avg_z_lower]
	puts $outfile "\t$diff_z"
}
close $outfile


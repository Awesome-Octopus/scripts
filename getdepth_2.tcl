#################################
# change the following to change how many residues there are in your protein
# and which frame number past which to average
#######################################################################
set startframe 30
set aanum 30
set outfile [open depthplot.dat w]


# you need to write your own function just to take the mean of an array
proc average {list} {
    expr {[tcl::mathop::+ {*}$list 0.0] / max(1, [llength $list])}
}

# make a header for variable columns
puts -nonewline $outfile "frame"
for { set s 1 } { $s < [expr $aanum + 1 ] } { incr s } {
	puts -nonewline $outfile "\tres_$s"	
}
puts $outfile "\tupper_po4\tlower_po4"

atomselect macro prot { name BB or name SC1 or name SC2 or name SC3 or name SC4 or name SCN or name SCP }

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
	# po4 z position to give height relative the center of the bilayer at 0
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
		
		# for each amino acid find the z coordinate and subtract the z coordinate avg of the PO4
		# for both the upper and lower leaflets
################################################################		
#### for  loop over chain A, B , ... 
################################################################
	for { set l 1 } { $l < [expr $aanum + 1] } { incr l } {
		set calpha($l) [atomselect top "name BB and resid $l" frame $i]
		
		#set z_upper_dist [expr [$calpha($l) get z] - $avg_z_upper]
		#set z_lower_dist [expr [$calpha($l) get z] - $avg_z_lower]
		set aaposition [expr [$calpha($l) get z] - $zavgpo4]
		
		#puts -nonewline $outfile "\t$aaposition"
		
		# consider the distance for each point to be the lower 
		# of the two distances between that point and the average 
		# for the upper and lower leaflet
		# if { [expr abs ($z_upper_dist)] > [expr abs ($z_lower_dist)] } {
		#	puts -nonewline $outfile "\t$z_lower_dist"
		# } else {
		#	puts -nonewline $outfile "\t$z_upper_dist"
		# }
		
	}
	#puts -nonewline $outfile "\t$avg_z_upper"
	#puts $outfile "\t$avg_z_lower"
		
	
}
close $outfile


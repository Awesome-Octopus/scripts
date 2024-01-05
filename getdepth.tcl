set startframe 500
set aanum 75
set outfile [open depthplot.dat w]


# you need to write your own function just to take the mean of an array! -stupid
proc average {list} {
    expr {[tcl::mathop::+ {*}$list 0.0] / max(1, [llength $list])}
}

#make a header for variable columns
puts -nonewline $outfile "frame"
for { set s 1 } { $s < [expr $aanum + 1 ] } { incr s } {
	puts -nonewline $outfile "\tres_$s"	
}
puts $outfile ""

atomselect macro prot { name BB or name SC1 or name SC2 or name SC3 or name SC4 or name SCN or name SCP }

	# get the number of frames
set nf [molinfo top get numframes]

	# loop over all frames from startframe to end
for { set i $startframe } { $i < [expr $nf + 1] } { incr i } {

		#find all PO4 within 7 angstroms of a protein bead
	set phos [atomselect top "name PO4 and within 7 of prot" frame $i]
	puts -nonewline $outfile "$i"
	set zavgpo4 [average [$phos get z]]
	

		# for each amino acid find the z coordinate and subtract the z coordinate avg of the PO4
	for { set l 1 } { $l < [expr $aanum + 1] } { incr l } {
		set calpha($l) [atomselect top "name BB and resid $l" frame $i]
		puts -nonewline $outfile "\t[expr [$calpha($l) get z] - $zavgpo4]"
		
	}
	puts $outfile ""
}

close $outfile


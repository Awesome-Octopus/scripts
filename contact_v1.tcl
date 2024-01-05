


# frame after which you want want to average over
set startframe 1
# number of amino acids in the protein
set aanum 56

# the name of the bead you want to count within proximity
set moiety PO4

# the name of the sidechain moiety to check
set moiety2 SCP

# the distance in Angstroms for the cutoff to be considered in contact with
set threshold 3.5


set outfile [open KR_contacts.dat w]


#LIST OF indeces of positive side chains of R and K residues
set resnums { 5 6 16 24 34 36 51 54 55 }

for { set w 0 } { $w < [llength $resnums] } { incr w } {
	 set sel [atomselect top "name $moiety2 and resid [lindex $resnums $w]"]
	 set SC_ndx($w) [$sel list]

}

#make a header for variable columns

# uncomment for optional frame number column on the left
#puts -nonewline $outfile "frame\t"
foreach x $resnums {
	puts -nonewline $outfile "$x\t"	
}
puts $outfile ""

# get the number of frames
set nf [molinfo top get numframes]



# loop over all frames from startframe to end
for { set i $startframe } { $i < [expr $nf + 1] } { incr i } {
	
	# uncomment for optional frame number column on the left
	# puts -nonewline $outfile "$i"
	

	# for every positive group in K and R sidechains
	for { set y 0 } { $y < [array size SC_ndx] } { incr y } {
		
		
		# select however many of the target moiety are within the threshold of your sidechain bead at that frame
		set sel [atomselect top "name $moiety and within $threshold of index $SC_ndx($y)" frame $i] 
		
		# print how many beads were selected this way
		puts -nonewline $outfile "[$sel num]\t"
		
	}
	puts $outfile ""
}
close $outfile


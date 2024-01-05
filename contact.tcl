# -------------------12/14/2022--------------------------
#---------------D-E-S-C-R-I-P-T-I-O-N--------------------
# This is a VMD script that will produce a text table
# showing, per frame, how many of a target bead is within
# a radius threshold of another kind of bead
# e.g. how many phosphates are within 3 A of a charged 
# sidechain group of each lysine in each frame
# the result is a contact map
# -------------------------------------------------------

#########################################################
# ------------variables to alter------------------------#
#########################################################
# frame after which you want want to average over, 
# set to 1 to consider every frame
set startframe 1

# what unique identifier to give the file
# currently based on the name of the folder
set unique [lindex [split [pwd] /] end-1]

# how long is the peptide sequence

# the name of the bead you want to count
# within proximity of your sidechain moiety
set query_moiety NC3
set query_moiety2 PO4

# these are your sidechain moieties to check
set moiety SCN
set moiety2 SCP

# the distance in Angstroms for the cutoff to be considered in contact with
set threshold 4

# the name of the text file containing the result
set outfile [open $moiety-$query_moiety-$moiety2-$query_moiety2-contact-$threshold-A-cutoff-$unique.dat w]

# lists of the indeces to your sidechain moeities which you will
# iterate through with calls to check your query moieties
set sel [atomselect top "type $moiety" frame 1]
set m1_ndx [$sel list]
set sel [atomselect top "type $moiety2" frame 1]
set m2_ndx [$sel list]

#-------------------make a header for variable columns-------------------


puts -nonewline $outfile "frame\t"
foreach ndx $m1_ndx {
	set sel [atomselect top "index $ndx" frame 1]
	puts -nonewline $outfile "[$sel get resname][$sel get resid][$sel get name]\t" 	
}

foreach ndx $m2_ndx {
	set sel [atomselect top "index $ndx" frame 1]
	puts -nonewline $outfile "[$sel get resname][$sel get resid][$sel get name]\t" 	
}

#needed to start a new line
puts $outfile {}

# get the number of frames
set nf [molinfo top get numframes]



# loop over all frames from startframe to end
for { set i $startframe } { $i < [expr $nf + 1] } { incr i } {
	
	# uncomment for optional frame number column on the left
	puts -nonewline $outfile "$i\t"
	
	#for each moiety index in each frame print how many of the target moiety are
	# within threshold distance in that frame
	foreach ndx $m1_ndx {
		set sel [atomselect top "index $ndx and within $threshold of name $query_moiety" frame $i] 	
		puts -nonewline $outfile "[$sel num]\t"	
	}
	
	foreach ndx $m2_ndx {
		set sel [atomselect top "index $ndx and within $threshold of name $query_moiety2" frame $i] 	
		puts -nonewline $outfile "[$sel num]\t"	
	}
	
	
	puts $outfile {}
}
close $outfile


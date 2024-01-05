proc distance_z {seltext1 N_d f_r_out f_d_out} {
# note N_d is the number of binds for the distribution

set sel1 [atomselect top "$seltext1"]

set nf [molinfo top get numframes]

##################################################
# Loop over all frames.                          #
##################################################
set outfile [open $f_r_out w]
for {set i 0} {$i < $nf} {incr i} {

  puts "frame $i of $nf"
  $sel1 frame $i

set center [measure center $sel1 weight mass]
set z_distance [lindex $center 2]
puts $outfile "$i $z_distance"
}
close $outfile

}

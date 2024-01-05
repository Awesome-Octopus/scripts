set fd1 [open prot_pnts.txt "w"]
set fd2 [open memb_pnts.txt "w"]
set prot [atomselect top "name BB or name SC1 or name SC2 or name SC3 or name SC4 or name SCN or name SCP"]

set memb [atomselect top "name NC3 or name PO4"]

puts $fd1 [$prot get z]
puts $fd2 [$memb get z]

close $fd1
close $fd2


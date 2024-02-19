#####
# use fix_martini_psf to clean up protein segments

# load all-atom protein 
set aamol [mol new AA.psf]
mol addfile AA.pdb waitfor all

# source script
source fix_martini_psf.tcl

#Usage: fix_martini_psf {aamol cgtopfile cgpsffile cgpdbfile outprefix {usedssp 0} {dssppath ""} {Ncharge 1} {Ccharge -1}}

# example using STRIDE
#fix_martini_psf $aamol ../../04-cgc-top-par-files/martini-top/martini-protein.top cg-ubiquitin-init.psf cg-ubiquitin-init.pdb cg-ubiquitin-fixed 0 0 0 0

#example using DSSP
fix_martini_psf $aamol ../../martini-protein.top cg-init.psf cg-init.pdb cg-fixed 1 . 0 0

exit

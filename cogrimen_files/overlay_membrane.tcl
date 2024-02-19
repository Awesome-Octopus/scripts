# get the filename for the dcd and psf files
regexp {(S2E-processed-cg-input-[0-9]+\.psf)} [ls "*.psf"] psf_file
regexp {(cogrimen-s2e-60-nm-memb-[0-9]+\.dcd)} [ls "*.dcd"] dcd_file
mol new $dcd_file type dcd waitfor all
mol addfile $psf_file type psf waitfor all
display rendermode Acrobat3D
color Display Background white
axes location Off
display projection Orthographic
rotate x by -90
# show only BB atoms
mol modselect 0 0 name BAS
# 5 frame trajectory smoothing
mol smoothrep 0 0 5
#make the bb trace white
mol modcolor 0 0 ColorID 16
mol modstyle 0 0 Lines 4.000000

mol selection name BAS and resid < 39 and resid > 11
mol addrep 0
mol modcolor 1 0 ColorID 1
mol modstyle 1 0 Licorice 1.000000 10.000000 10.000000
mol smoothrep 0 1 5

mol selection name BAS and resid < 46 and resid > 39
mol addrep 0
mol modcolor 2 0 ColorID 7
mol modstyle 2 0 Licorice 1.000000 10.000000 10.000000
mol smoothrep 0 2 5

mol selection name BAS and resid < 64 and resid > 52
mol addrep 0
mol modcolor 3 0 ColorID 21
mol modstyle 3 0 Licorice 1.000000 10.000000 10.000000
mol smoothrep 0 3 5


animate speed 0.7

# this is a special vmd templated custom proc for making graphics objects
proc membrane {width} {
    set bottom {}
    lappend bottom -40 -40 [expr -1*$width/2]
    set top {} 
    lappend top -40 -40 [expr $width/2]
    draw materials on
    draw material Glass3
    draw color gray
    draw cylinder $bottom $top radius 80 filled yes 
}
membrane {40}

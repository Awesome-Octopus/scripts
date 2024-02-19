#martini-tutorial - script to fix martini psf - July 16 2011 
# files required aa.pdb aa.psf cg.pdb cg.psf 
#To use:
# load in vmd aa.pdb aa.psf
#source fix_martini_psf.tcl


proc fix_martini_psf {aamol cgtopfile cgpsffile cgpdbfile outprefix {usedssp 0} {dssppath ""} {Ncharge 1} {Ccharge -1}} {
# Generate a correct psf for a martini structure
# Apply the necessary patches for secondary structure in the process
# We assume that segment names in each file are identical

# If usedssp is nonzero, DSSP will be used instead of STRIDE
# In this case, DSSPPATH must be set if DSSP is not in your default path

# Since charge on C and N terminus should be treated independent of each other, we have the options separately
# If Ncharge is any number other than +1, N termini is uncharged 
# If Ccharge is any number other than -1, C termini is uncharged

    puts "\n\nApplying fix_martini_psf"
    puts "Usage: fix_martini_psf {aamol cgtopfile cgpsffile cgpdbfile outprefix {usedssp 0} {dssppath \"\"} {Ncharge 1} {Ccharge -1}}"
    puts "Command entered: fix_martini_psf $aamol $cgtopfile $cgpsffile $cgpdbfile $outprefix $usedssp $dssppath $Ncharge $Ccharge\n\n"

  set cgpatchlist [list]

  set aaprot [atomselect $aamol "protein and name CA"]
  set numprotres [$aaprot num]
puts "We have $numprotres residues."
  if {$numprotres == 0} {
    puts "No protein residues detected! Aborting script.\n\n"
    exit 0
  }

  set chainnames [lsort -unique [$aaprot get segname]]
  $aaprot delete

  foreach chain $chainnames {
    set sel [atomselect $aamol "segname $chain"]

    # DSSP seems to be confused by CTER residues, claiming the backbone is incomplete for that residue
    # This seems to be because there is no atom named O (instead there are OT1 and OT2)
    # Therefore we check to see if there are any atoms named OT1, and if so, rename them to O for tmp-martini.pdb
    # which will be fed into DSSP
    # There does not seem to be any problem with NTER residues.
    set cterchecksel [atomselect $aamol "segname $chain and name OT1"]
    $cterchecksel set name O
    $sel writepsf tmp-martini.psf
    $sel writepdb tmp-martini.pdb
    $sel delete
    set tmpmol [mol new tmp-martini.psf]
    mol addfile tmp-martini.pdb waitfor all

    set sel [atomselect $tmpmol "segname $chain and name CA"]
    set numres [$sel num]
    if {$numres == 1} {
    # check to see if this is a free residue
    # if so, don't call STRIDE, since it only returns an error
    # just apply the appropriate free residue patch based on the resname
      puts "This amino acid is a free residue."
      set ss ""
      set pss ""
      set nss ""
      set resn [$sel get resname]
      set resi [$sel get resid]
      set sspatch [get_patch_id $ss $pss $nss $resn $resi $chain]
      if {$sspatch != {}} {lappend cgpatchlist $sspatch}
      puts "\nReconstructing PSF file... \n"
      package require psfgen
      resetpsf

      # (this should be the Martini protein topology file)
      topology $cgtopfile
      readpsf $cgpsffile
      coordpdb $cgpdbfile

      # apply all the patches in the list
      foreach patch $cgpatchlist {
        patch [lindex $patch 0] [lindex $patch 1]
      }

      writepsf ${outprefix}.psf
      writepdb ${outprefix}.pdb

      exit 1
    }

    if {$usedssp == 0} {
      set ssarr [$sel get structure]
      puts "\nFinal contents of ssarr using STRIDE:"
      puts "$ssarr \n"
    } else {
       set sellen [$sel num]
      	if {$sellen == 1 } {
	set ssarr [$sel get structure]		
        } else {
        set ssarr [get_ss_from_dssp "$dssppath" "tmp-martini.pdb"]
        puts "\nFinal contents of ssarr using DSSP:"
        puts "$ssarr \n"
        #puts "**end**"
        }
    }
    set resnarr [$sel get resname]
    set resiarr [$sel get resid]

    for {set i 0} {$i < [llength $ssarr]} {incr i} {
      set ss [lindex $ssarr $i]

#     puts "************does this print on using dssp ?!****"
#      puts "this residue index no: "
#      puts $i
#      puts $ss
#      puts "***************"


      set resn [lindex $resnarr $i]
      set resi [lindex $resiarr $i]
#puts "ssarr is [llength $ssarr] long"
#puts "Looking at residue $i"
      # test for a free amino acid chain
      if {$i == 0 && $i == [expr [llength $ssarr] -1]} {

        #puts "This amino acid is a free residue."
        puts "You should never see this section, since we already checked for the free residue case outside of this loop."
        #set ss [lindex $ssarr $i]
	#set pss ""
	#set nss ""

      	#set sspatch [get_patch_id $ss $pss $nss $resn $resi $chain]
      	#if {$sspatch != {}} {lappend cgpatchlist $sspatch}
	#continue
      }

      # check to see if this is the first residue and apply patch
      if {$i == 0} {
        puts "This is the first residue of segment ${chain}."
        set pss ""
        # if so, append the correct n-terminal patch to the list
	if {$Ncharge == 1} { 
		lappend cgpatchlist [list NT${ss} ${chain}:$resi]
        	#puts "A charged N-terminus was requested."
        	puts "Secondary structure is ${ss}.  Applying patch [list NT${ss} ${chain}:$resi] to segname ${chain} resid ${resi}. (charged N-terminus)"
	} else {
               lappend cgpatchlist [list NU${ss} ${chain}:$resi]
        	#puts "An ucharged N-terminus was requested."
        	puts "Secondary structure is ${ss}.  Applying patch [list NU${ss} ${chain}:$resi] to segname ${chain} resid ${resi}. (uncharged N-terminus)"
	}
	continue
      } else {
        set pss [lindex $ssarr [expr $i - 1]]
      }

      # check to see if this is the last residue and appy patch
      if {$i == [expr [llength $ssarr] - 1] } {
        set nss ""
        puts "This is the last residue of segment ${chain}."
        # if so, append the correct c-terminal patch to the list
        if {$Ccharge == -1} {
       		lappend cgpatchlist [list CT${ss} ${chain}:$resi]
        	#puts "A charged C-terminus was requested."
        	puts "Secondary structure is ${ss}. Applying patch [list CT${ss} ${chain}:$resi] to segname ${chain} resid ${resi}. (charged C-terminus) \n"
	} else {
		lappend cgpatchlist [list CU${ss} ${chain}:$resi]
        	#puts "An uncharged C-terminus was requested."
        	puts "Secondary structure is ${ss}. Applying patch [list CU${ss} ${chain}:$resi] to segname ${chain} resid ${resi}. (uncharged C-terminus) \n"
	}
        continue
      } else {
        set nss [lindex $ssarr [expr $i + 1] ]
      }

      # if this is a residue somewhere else in the chain, call get_patch_id to figure out the right patch and append it to the list
      set sspatch [get_patch_id $ss $pss $nss $resn $resi $chain]
      puts "Secondary structure is ${ss}.  Applying patch $sspatch to segname ${chain} resid ${resi}."
      if {$sspatch != {}} {lappend cgpatchlist $sspatch}
    }
    $sel delete
    mol delete $tmpmol
  }


  #####
  ### Reconstruct the psf file with the secondary structure patches.
  #####
  puts "\nReconstructing PSF file... \n"
  package require psfgen
  resetpsf

  # (this should be the Martini protein topology file)
  topology $cgtopfile
  readpsf $cgpsffile
  coordpdb $cgpdbfile

  # apply all the patches in the list
  foreach patch $cgpatchlist {
    patch [lindex $patch 0] [lindex $patch 1]
  }

 # set execstring "autopsf -mol $cgmol -top $cgtopfile -prefix $outprefix [join $cgpatchlist]"
 # puts $execstring
 # eval $execstring

  writepsf ${outprefix}.psf
  writepdb ${outprefix}.pdb

  file delete -force -- tmp-martini.pdb
  file delete -force -- tmp-martini.psf
}

proc get_ss_from_dssp {dssppath pdbfile} {
  puts "\n--start of get_ss_from_dssp routine--"
  set dssplogfile [open "fix_martini_dssp_results.log" w+]
  puts "The full DSSP log will be written to fix_martini_dssp_results.log"
  set dsspstream [open "| [file join ${dssppath} dsspcmbi] -v $pdbfile"]
  set ssarr [list]
  set dsspRefarr [list]
  
  # readingSSdata is a flag indicating whether the current line contains data, vs. unwanted header information or an empty line
  set readingSSdata 0

  while {[gets $dsspstream line] >= 0} {
    puts $dssplogfile $line
#    set findout [string index $line 16]
#    puts $findout
    #puts $readingSSdata

#prints the 17 th character of $line
    #(the 17th character of $line, which should contain the secondary structure)

    if {$readingSSdata == 1} {
      if {$line == ""} {
        set readingSSdata 0
        continue
      }
      #puts "ssarr before dssp-> stride: $ssarr"
      lappend ssarr [conv_dssp_to_stride_ss [string index $line 16]]
      lappend dsspRefarr [string index $line 16]
      #puts "ssarr after dssp -> stride: $ssarr \n"
    } elseif {[string equal -length 10 "  #  RESIDUE AA" "[string range $line 0 14]"]} {
      # this marks the end of the header and the beginning of the secondary structure data
      set readingSSdata 1
    }
  }
  catch {close $dsspstream}
  catch {close $dssplogfile}
  puts "secondary structure from DSSP: $dsspRefarr"
  puts "after conversion to STRIDE   : $ssarr"
  puts "--end of get_ss_from_dssp routine--\n"

  return $ssarr
}





# Listed below are the label conventions for STRIDE and DSSP:

#  DSSP:                                     STRIDE:
#  H = alpha helix                           H = alpha helix
#  I = pi helix                              I = pi helix
#  G = 3-10 helix                            G = 3-10 helix
#  T = hydrogen-bonded turn                  T = turn
#  E = extended strand (beta sheet)          E = extended conformation/beta sheet
#  B = residue in isolated beta bridge       B = isolated bridge
#  S = bend                                  C = coil, i.e. none of the above
#" " = none of the above, or error


# convert a single secondary structure entry from dssp to stride convention
proc conv_dssp_to_stride_ss {dsspstruct} {
  switch -exact $dsspstruct {
    "S" {return "B"}
    "B" {return "E"}
    " " {return "C"}
    default {return $dsspstruct}
  }
}



# Return the patch statement needed to correct the BAS bead type in psfgen for this residue.
# pss and nss are the secondary structures of the previous and next residues.
proc get_patch_id {ss pss nss resname resid segname} {

  set firsttwo [get_firsttwo $resname]
  set helstart 0
  set helend 0
  # if the next residue is not helical, this might be the end of a helix
  if {$nss != "H" && $nss != "I" && $nss != "G"} {set helend 1}
  # if the previous residue was not helical, this might be the start of a helix
  if {$pss != "H" && $pss != "I" && $pss != "G"} {set helstart 1}
  # if $nss is empty, either this is the last residue, or this is a free residue
  if {$nss == {}} {set helend 2}
  # if $pss is empty, either this is the first residue, or this is a free residue
  if {$pss == {}} {set helstart 2}
  # if both flags are set to 2, this is a free residue, apply free residue patch
  if {$helstart > 1 && $helend > 1} {
    return [list ${firsttwo}F $segname:$resid]
  }

  switch -exact $ss {
    "C" { set lastone "C"}
    "B" { set lastone "B"}
    "E" { set lastone "E"}
    "T" { set lastone "T"}
    "H" -
    "G" -
    "I" { 
          set lastone "H"
          if {$helstart > 0} { set lastone "N" } 
          if {$helend > 0} { set lastone "I" } 
    }
  }

  return [list ${firsttwo}${lastone} $segname:$resid]

}



# return an appropriate beginning for a CG backbone patch, given the residue name
proc get_firsttwo {name} {
  switch -exact $name {
    "GLY" {return "GC"}
    "PRO" {return "PC"}
    "ALA" {return "AC"}
    default {return "OC"}
  }
}



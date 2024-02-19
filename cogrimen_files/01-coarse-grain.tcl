
# VMD won't let you redefine "protein"
# and it does not recognize CG residues as being protein residues, even though they have the appropriate resname.
# so we can create an atomselect macro called "cgprotein" to deal with this
atomselect macro cgprotein {resname ALA ARG ASN ASP CYS GLN GLU GLY HSD HSE HSP HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL}

package require cgtools

set aamol [mol new AA.psf]
mol addfile AA.pdb waitfor all

# read in .cgc files
::cgtools::read_db martini-protein.cgc

# coarse grain
::cgtools::apply_database $aamol cg.pdb cg.rcg

exit


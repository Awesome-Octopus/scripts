mol load pdb randomized_S2E_0.pdb
set protsel [atomselect top "protein"]
$protsel writepdb temp.pdb

package require psfgen
resetpsf
pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD

topology all27_prot_lipid_cmap.top

segment P {
  pdb temp.pdb
  first nter
  last cter
}
coordpdb temp.pdb P

guesscoord
writepsf AA.psf
writepdb AA.pdb

file delete -force -- temp.pdb

exit

# Define a procedure to find and select residues based on atom name and 
# distance from a point. The selects all molecules that contain at least
# one atom named "atom_name" within "distance" angrom radius from the
# atom with index value "ndx"

proc whole_mol_sel_by_radius {atom_name distance ndx} {
    # Select atoms of the given name within the distance threshold
    set selection [atomselect top "name $atom_name and within $distance of index $ndx"]

    # Get the list of unique residue IDs (resids)
    set resids [lsort -unique [$selection get resid]]

    # Create a selection string for these residues
    set selection_string ""
    foreach resid $resids {
        append selection_string "resid $resid or "
    }

    # Remove the trailing ' or '
    set selection_string [string trimright $selection_string " or "]

    # Return the selection string
    return $selection_string
}

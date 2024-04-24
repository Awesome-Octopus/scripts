#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 12:31:06 2024

@author: andrew
"""


def check_clash(cutoff_distance=0.36, angle_type='phi', anchor='N', **kwargs):
    """
    Checks if a rotation arround the phi or psi bond (default: phi) of the
    input residue has resulted in atoms closer than the cutoff radius (in nm)
    Return any atom serial numbers that clash. To be invoked after performing
    rotation. The cutoff distance is based on 2* the carbon vand der waals
    radius
    """

    # ---------------------------- CAVEATS ------------------------------------
    # kwargs should only point to a SINGLE residue. If kwargs match atoms on
    # multiple residues or the same residue on different chains, undefined
    # behavior will result. I COULD write it in to check that this is the case,
    # but it seems computationally expensive to check every time.
    # ------------- IT IS ON THE USER TO PASS KWARGS PROPERLY -----------------

    # extract the residue number from kwargs, if an atom name was specified
    # remove it
    res = kwargs.pop('res_num')
    if 'atom_name' in kwargs.keys():
        kwargs.pop('atom_name')

    # everything to the n-terminal side of the bond will not clash with others
    # which are also on that side of the bond. check each atom on the n
    # terminal side with every atom that isn't n-terminal to that bond
    # we take advantage of the fact that these should all be contiguous in
    # memory except hydrogen on the residue in question

    same_chain_atms = select(**kwargs)
    kwargs['res_num'] = res

    if angle_type == 'phi':

        kwargs['atom_name'] = 'N'   # get amide nitrogen
        amide_nitrogen = select(**kwargs)
        if len(amide_nitrogen) > 1:
            raise ValueError

        # split the chain onto the atoms on one side of the bond and those
        # on the other side
        n_term_frag = same_chain_atms[0:amide_nitrogen[0]+1]
        c_term_frag = same_chain_atms[amide_nitrogen[0]+1:]
        try:
            kwargs['atom_name'] = 'HN'
            amide_hydrogen = select(**kwargs)
            n_term_frag.append(amide_hydrogen[0])
            # add on the amide hydrogen, whose index is not contiguous with
            # the other n-terminal atoms, if it exists and remove it from the c
            # terminal list

            c_term_frag.remove(amide_hydrogen[0])

        except ValueError:
            pass

    elif angle_type == 'psi':

        # select the amide carbon and oxygen from the current residue and add
        # them to the list composed of all serial numbers from the next amide
        # nitrogen onwards
        kwargs['atom_name'] = 'C'
        amide_carbon = select(**kwargs)

        if len(amide_carbon) > 1:
            raise ValueError

        kwargs['atom_name'] = 'N'
        kwargs['res_num'] = res + 1
        next_amide_nitro = select(**kwargs)
        kwargs['res_num'] -= 1
        n_term_frag = same_chain_atms[0:next_amide_nitro[0]]
        c_term_frag = same_chain_atms[next_amide_nitro[0]:]

        # remove the amide carbon and oxygen from n terminal side and add to
        # c terminal side, again, because they are discontiguous
        try:

            kwargs['atom_name'] = 'C'
            amide_carbon = select(**kwargs)
            c_term_frag.append(amide_carbon[0])
            n_term_frag.remove(amide_carbon[0])
            kwargs['atom_name'] = 'O'
            amide_oxygen = select(**kwargs)
            c_term_frag.append(amide_oxygen[0])
            n_term_frag.remove(amide_oxygen[0])

        except ValueError:
            pass

    else:    # if input was invalid
        raise ValueError
        print(
            f'{angle_type} is not a valid input for parameter angle_type. '
            f'Enter either phi or psi')

    clash_atom_ser_nums = []
    search_list = []
    for r in info_table:
        search_list.append(r['ser_num'])

    if anchor == 'N':

        # remove the c-terminal fragment from the searchlist
        # so it doesn't search against itself
        for atom in c_term_frag:
            search_list.remove(atom)
        # then search each atom in the fragment against every atom not a part
        # of that fragment
        for atom in c_term_frag:
            for query in search_list:
                dist = coordinates[atom].distance_between(coordinates[query])
                if dist < cutoff_distance:
                    clash_atom_ser_nums.append([atom, query])

    elif anchor == 'C':

        # remove the c-terminal fragment from searchlist
        # so it doesn't search against itself
        for atom in n_term_frag:
            search_list.remove(atom)
        # then search each atom in the fragment against every atom not a part
        # of that fragment
        for atom in n_term_frag:
            for query in search_list:
                dist = coordinates[atom].distance_between(coordinates[query])
                if dist < cutoff_distance:
                    clash_atom_ser_nums.append([atom, query])

    else:
        raise ValueError
        print(
            f'{anchor} is not a valid input for parameter anchor. Enter '
            f'either N or C')
    return clash_atom_ser_nums

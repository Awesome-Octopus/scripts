#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 01:36:49 2023

@author: andrew
"""

# <codecell> Header
import random
import os
import matplotlib as mpl
from Bio.PDB.PDBParser import PDBParser
# how we will import the raw text
# from Bio.PDB import PDBIO

import numpy as np
import math
parser = PDBParser(PERMISSIVE=1)
# permissive=1 allows it to accept incomplete structures without error


# structure_id = "mode_7"
# whole_pdb = parser.get_structure("randomized_S2E",
#                                  "/home/andrew/test-cogrimen/randomized_S2E_0.pdb")
info_table = []


class vector (np.ndarray):
    def __new__(cls, vect):
        if type(vect) == list:
            vect = np.asarray(vect)
        obj = vect.view(cls)
        return obj

    def get_length(self):
        """
        Get the norm of the vector.

        Returns
        -------
        Float64
            DESCRIPTION.

        """
        return np.linalg.norm(self)

    def unitize(self):
        """
        Get the unit vector.

        Raises
        ------
        ZeroDivisionError
            If user tries to unitize zero vector.

        Returns
        -------
        vector
            A vector of length 1 colinear with the input vector.

        """
        length = self.get_length()
        if length == 0:
            raise ZeroDivisionError(
                'Error in vector.unitize(): Attempt to unitize the zero '
                'vector\n')
        else:
            return (self/length)

    def is_nonzero(self):
        if self.get_length() == 0:
            return False
        else:
            return True

    def dot(self, other_vector):
        return np.dot(self, other_vector)

    def distance_between(self, other_vector):
        d = self - other_vector
        return d.get_length()

    def angle_between(self, other_vector):

        try:
            a = self.unitize()
            b = other_vector.unitize()
            c = a.dot(b)
            return np.arccos(c)
        except ZeroDivisionError as err:
            print(
                'Error in vector.angle_between(): Attempt to take an angle '
                'with the zero vector\n',
                err)

    def is_orthogonal(self, other_vector):
        if self.angle_between(other_vector) == np.pi/2:
            return True
        else:
            return False

    def rotate_arround(self, angle, axis_vector=[]):

        # rotate a 2d or a 3d vector arround another 2d or 3d vector by the
        # given angle according to right hand rule. intentionally
        # not defined for vectors of rank 4 or greater

        if len(axis_vector) != 0:
            if axis_vector.is_nonzero() is True:
                axis_vector = axis_vector.unitize()
                if axis_vector.size == 3 and self.size == 3:
                    x = axis_vector[0]
                    y = axis_vector[1]
                    z = axis_vector[2]
                    a = np.cos(angle)
                    b = 1 - a
                    c = np.sin(angle)
                    m1 = a + (x**2)*b
                    m2 = x*y*b - z*c
                    m3 = x*z*b + y*c
                    m4 = y*x*b + z*c
                    m5 = a + (y**2)*b
                    m6 = y*z*b - x*c
                    m7 = z*x*b - y*c
                    m8 = z*y*b + x*c
                    m9 = a + (z**2)*b
                    # these terms describe a general rotation matrix for a 3d
                    # vector arround another arbitrary vector
                    rot_mat = np.asarray([[m1, m2, m3],
                                          [m4, m5, m6],
                                          [m7, m8, m9]])
                    return np.matmul(self, rot_mat)
                else:
                    raise ValueError(
                        'Error in vector.rotate_arround(): Dimensions of '
                        'vectors do not agree or vectors are neither rank 2 '
                        'nor rank 3\n')
            else:
                raise ValueError(
                    'Error in vector.rotate_arround(): Attempt to rotate '
                    'arround the zero vector\n')
        elif self.size == 2:
            rot_mat = np.asarray([[np.cos(angle), -1*np.sin(angle)],
                                  [np.sin(angle), np.cos(angle)]])
            return np.matmul(self, rot_mat)
        else:
            raise ValueError(
                'Error in vector.rotate_arround(): No axis of rotation '
                'provided for a rank 3 vector\n')

    def project_onto_normal_plane(self, other_vector):

        # given some vector and a second vector return a vector which is the
        # component of the first vector which is orthogonal to the second, in
        # other words, if the second vector defined the normal of a plane, give
        # the shadow of the first vector projected onto that plane

        a = self.dot(other_vector.unitize())
        shadow = self - a*other_vector
        return shadow


# <codecell> Functions


def select(**kwargs):

    # this function will return the indices in a list of all rows in the info
    # list who match the input arguments
    global info_table
    serial_numbers = []
    info = info_table
    for key, value in kwargs.items():
        if isinstance(value, list):
            info = [s for s in info if s[key] in value]
            # cull the list down to just items that match the criterion
            # then feed the culled list back through for the next non-'any'
            # criterion
        else:
            info = [s for s in info if s[key] == value]
    for i in info:
        serial_numbers.append(i['ser_num'])
        serial_numbers = list(set(serial_numbers))
    # return the serial numbers for everything that survived the culling
    # convert to a set to eliminate redunancy, then back to a list for indexing
    if not serial_numbers:
        raise ValueError(
            'Error in select: no atom(s) found that match criteria')
    return serial_numbers


def get_phi_psi(res_num=None):
    """
    Update the phi and psi values in the info_table for the specified residues.

    Parameters
    ----------
    res_num : TYPE, optional
        DESCRIPTION. The default is None

    Returns
    -------
    None.

    """
    selected_residues = []

    if res_num is None:
        res_num = []
        ser_nums = select(atom_name='CA')
        for num in ser_nums:
            res_num.append(info_table[num]['res_num'])
    else:
        if isinstance(res_num, int):
            res_num = [res_num]

    for n in res_num:
        n_ser_num = select(atom_name='N', res_num=n)
        ca_ser_num = select(res_num=n, atom_name='CA')
        c_ser_num = select(res_num=n, atom_name='C')
        ca_to_c_bond_vect = coordinates[c_ser_num] - coordinates[ca_ser_num]
        ca_to_c_bond_vect = ca_to_c_bond_vect[0]

        # for some reason subtracting vectors returns a list of
        # vectors with a single element instead of the vector itself making
        # the [0] indexing necessary

        ca_to_n_bond_vect = coordinates[n_ser_num] - coordinates[ca_ser_num]
        ca_to_n_bond_vect = ca_to_n_bond_vect[0]

        try:
            x = n - 1
            # if there is no residue that has a number higher than ours we are
            # at the c terminus, for which there is no psi value, otherwise you
            # need the next amide nitrogen
            prev_c_ser_num = select(res_num=x, atom_name='C')
            n_to_prev_c_bond_vect = coordinates[prev_c_ser_num] - \
                coordinates[n_ser_num]
            n_to_prev_c_bond_vect = n_to_prev_c_bond_vect[0]
            projection_a = n_to_prev_c_bond_vect.project_onto_normal_plane(
                ca_to_n_bond_vect)
            projection_b = ca_to_c_bond_vect.project_onto_normal_plane(
                -1*ca_to_n_bond_vect)
            phi = projection_a.angle_between(projection_b)
        except ValueError:
            phi = float('nan')
        try:
            x = n + 1
            # if there is no residue that has a number higher than ours we are
            # at the c terminus, for which there is no psi value, otherwise you
            # need the next amide nitrogen
            next_n_ser_num = select(res_num=x, atom_name='N')
            c_to_next_n_bond_vect = coordinates[next_n_ser_num] - \
                coordinates[c_ser_num]
            c_to_next_n_bond_vect = c_to_next_n_bond_vect[0]
            projection_a = ca_to_n_bond_vect.project_onto_normal_plane(
                -1*ca_to_c_bond_vect)
            projection_b = c_to_next_n_bond_vect.project_onto_normal_plane(
                ca_to_c_bond_vect)
            psi = projection_a.angle_between(projection_b)
        except ValueError:
            psi = float('nan')
        selected_residues = select(res_num=n)
        for i in selected_residues:
            info_table[i]['phi'] = phi
            info_table[i]['psi'] = psi


def set_phi_psi(angle, angle_type='phi', anchor='N', **kwargs):
    """
    Rotate everything in the coordinates matrix so that the residue given has a
    phi or a psi angle as specified, with the angle specified in degrees by the
    angle parameter. Which end will be held static while the other is rotated
    is specified by the anchor parameter
    """

    global info_table
    global coordinates

    if 'ser_num' in kwargs.keys():
        if isinstance(kwargs['ser_num'], list):
            raise ValueError(
                'you tried to set phi psi values by referencing more than one'
                ' atom serial number as keyword arguments')
        else:
            kwargs['res_num'] = info_table[kwargs['ser_num']]['res_num']
        del kwargs['ser_num']
    # if a serial number for an atom is passed, instead refer to the residue
    # number that that atom is on

    kwargs['atom_name'] = ['CA']
    # if an atom name was specified ignore it
    # and search for 'CA'
    ca_ser_nums = select(**kwargs)

    for ndx, ca in enumerate(ca_ser_nums):
        # for every ca on residues where we want to change the angle

        res = info_table[ca_ser_nums[ndx]]['res_num']
        mod = info_table[ca_ser_nums[ndx]]['model_name']
        sub = info_table[ca_ser_nums[ndx]]['submodel']
        chn = info_table[ca_ser_nums[ndx]]['chain']
        # get its info

        ser_nums_on_same_chain = select(
            model_name=mod, submodel=sub, chain=chn)
        # get all the atoms on the same chain as it

        n_ser_num = select(res_num=res, model_name=mod,
                           submodel=sub, chain=chn, atom_name='N')[0]
        # get the amide nitrogen on the same residue

        c_ser_num = select(res_num=res, model_name=mod,
                           submodel=sub, chain=chn, atom_name='C')[0]
        # and the carboxyl carbon on the same residue

        # save so we can add back later after rotation
        start_ca_position = coordinates[ca].copy()

    ############################     PHI    ###################################
        if angle_type == 'phi':
            current_phi = info_table[ca]['phi']

            if math.isnan(current_phi):

                get_phi_psi(res_num=res)
                current_phi = info_table[ca]['phi']
                # if the angle hasn't been calculated yet for this residue,
                # do so

            if info_table[ca]['res_name'] != 'PRO' and \
                    math.isnan(current_phi) is False:
                # we can't rotate the phi bond on a proline or for 1st residue,
                # so check to make sure it isn't one of these
                rotate_by = (angle - current_phi)

                # find the difference between the desired phi angle and current
                # phi angle in
                # degrees

                if anchor == 'N':

                    end_pt = ser_nums_on_same_chain[-1]

                    ###########################################################
                    # the last atom on our subset that is on the same chain.
                    # Note that this assumes that atoms are listed in order n
                    # to c terminus as is pdb convention we need all this to
                    # know what the range of serial numbers to apply the
                    # rotation to will be. Each time we do a rotation, we only
                    # want to rotate residues on the same model, submodel, and
                    # chain that are either n terminal (with c anchor, in which
                    # case all serial numbers up to ours) or c terminal
                    # (with n anchor, in which case starting with ours up to
                    # the last on this subset)
                    ###########################################################

                    coordinates[n_ser_num:end_pt + 1, :] -= start_ca_position
                    # center the origin on the alpha carbon only for the
                    # rotating subset. note that this assumes that atoms are
                    # listed in order n to c terminus as is pdb convention.

                    n_coord = coordinates[n_ser_num]
                    # the n coordinate on the same residue now defines the
                    # vector for the bond ca to n

                    NCA_bond_vect = n_coord.unitize()
                    # when the vector is reversed and unitized this defines our
                    # axis of rotation, which is ca to n

                    for n, vect in enumerate(coordinates[n_ser_num: end_pt + 1, :]):
                        coordinates[(n_ser_num + n), :] = vect.rotate_arround(
                            rotate_by, NCA_bond_vect)
                    for n, vect in enumerate(coordinates[n_ser_num:
                                                         end_pt+1, :]):
                        coordinates[(n_ser_num + n), :] = vect + \
                            start_ca_position

                elif anchor == 'C':

                    start_pt = ser_nums_on_same_chain[0]
                    # get the first atom serial number that is on the same
                    # model, submodel, and chain

                    coordinates[start_pt:ca + 1, :] -= start_ca_position
                    # center the origin on the alpha carbon only for the
                    # rotating subset. note that this assumes that atoms are
                    # listed in order n to c terminus as is pdb convention.
                    n_coord = coordinates[n_ser_num]
                    NCA_bond_vect = -1*n_coord.unitize()
                    # when the vector is reversed and unitized this defines our
                    # axis of rotation, which is n to ca

                    for n, vect in enumerate(coordinates[start_pt:ca, :]):
                        coordinates[(start_pt + n), :] = vect.rotate_arround(
                            rotate_by, NCA_bond_vect)
                    coordinates[start_pt:ca + 1] += start_ca_position
                else:
                    raise ValueError(
                        'the 3rd positional argument for set_phi_psi must be '
                        'either N or C, if none is given N is taken as the '
                        'default')
            else:
                raise ValueError('Warning: the phi value for residue',
                                 info_table[ca_ser_nums[ndx]]['res_num'],
                                 'can not be set because it is either a proline'
                                 ' or the n-terminus. The request was ignored.')

    ############################     PSI     ##################################

        elif angle_type == 'psi':
            current_psi = info_table[ca]['psi']
            if math.isnan(current_psi):
                current_psi = get_phi_psi(res_num=res)
                current_psi = info_table[ca]['psi']
                # if the angle hasn't been calculated yet for this residue,
                # do so

            if math.isnan(current_psi) is False:

                # we can't rotate the psi bond for the last residue, so check
                # to make sure it isn't

                rotate_by = (angle - current_psi)
                # find the difference between the desired phi angle and
                # current phi angle in degrees

                if anchor == 'N':

                    end_pt = ser_nums_on_same_chain[-1]
                    # the last atom on our subset that is on the same chain.
                    # note that this assumes that atoms are listed in order
                    # n to c terminus as is pdb convention

                    for n, vect in enumerate(coordinates[ca: end_pt + 1, :]):

                        coordinates[(ca + n), :] = vect - start_ca_position

                    # center the origin on the alpha carbon only for the
                    # rotating subset. note that this assumes that atoms are
                    # listed in order n to c terminus as is pdb convention.
                    c_coord = coordinates[c_ser_num]
                    # the c coordinate on the same residue now defines the
                    # vector for the bond ca to c
                    CAC_bond_vect = -1*c_coord.unitize()
                    # when the vector is r unitized this defines our axis of
                    # rotation, which is ca to c
                    for n, vect in enumerate(coordinates[c_ser_num: end_pt + 1, :]):
                        coordinates[(
                            c_ser_num + n), :] = vect.rotate_arround(rotate_by, CAC_bond_vect)
                    for n, vect in enumerate(coordinates[ca: end_pt + 1, :]):
                        coordinates[(ca + n), :] = vect + start_ca_position

                elif anchor == 'C':

                    start_pt = ser_nums_on_same_chain[0]
                    # get the first atom serial number that is on the same
                    # model, submodel, and chain

                    coordinates[start_pt:c_ser_num+1, :] -= start_ca_position
                    # center the origin on the alpha carbon only for the
                    # rotating subset. note that this assumes that atoms are
                    # listed in order n to c terminus as is pdb convention.
                    c_coord = coordinates[c_ser_num]
                    CAC_bond_vect = c_coord.unitize()
                    # when the vector is unitized this defines our axis of
                    # rotation, which is ca to n
                    for n, vect in enumerate(coordinates[start_pt: c_ser_num + 1, :]):
                        coordinates[(
                            start_pt + n), :] = vect.rotate_arround(rotate_by,
                                                                    CAC_bond_vect)
                    coordinates[start_pt:c_ser_num + 1] += start_ca_position
            else:
                raise ValueError('Warning: the psi value for residue',
                                 info_table[ca]['res_num'],
                                 'can not be set because it is at the c-terminus. '
                                 'The request was ignored.')
        else:
            raise ValueError(
                f'Warning: the dihedral angle type "{repr(angle_type)}" is '
                'not a valid selection. It must be either phi or psi. If '
                'none is given, phi is default')


def write_pdb(outfile=None):
    """
    Output a pdb text file of all currently loaded models with a default
    filename if none is given

    Returns
    -------
    None.

    """
    if outfile is None:
        model_name = info_table[0]['model_name']
    else:
        model_name = outfile

    # if the filename is already taken augment the number at the end of the
    # filename

    try:
        f = open(f'{model_name}', 'x')
    except FileExistsError:
        count = 0
        while True:
            try:
                f = open(f'{model_name}_{count}', 'x')
            except FileExistsError:
                count = count + 1
            else:
                break

    f.write(
        'MODEL        1                                                 '
        '               \n')
    num_termini = 0
    new_ser_num = 0
    model_count = 1
    # we need to keep track of how many times we hit the end of a chain because
    # the TER line gets its own serial number and we have to increase all
    # subsequent serial numbers by 1 for each TER we have encountered

    for n, a in enumerate(info_table):
        new_ser_num = new_ser_num + 1
        element = a['atom_name'][0]

        # the element is always the first letter in the name eg HG1 for gamma
        # hydrogen

        if element not in ['C', 'N', 'O', 'S', 'H']:

            # consider anything that is not one of these a heteroatom, note
            # this doesn't include phosphorous
            tag = 'HETATM'

        else:
            tag = 'ATOM  '

        if n != 0:
            if info_table[n-1]['chain'] != a['chain']:

                # if the chain for this is different from the last

                num_termini = num_termini + 1
                f.write(
                    f'TER   {n+num_termini:>5d}      '
                    f'{info_table[n-1]["res_name"]:<} '
                    f'{info_table[n-1]["chain"]}'
                    f'{info_table[n-1]["res_num"]: > 4d}'
                    f'{chr(10)}')
                new_ser_num = new_ser_num + 1
            if info_table[n-1]['submodel'] != a['submodel']:
                model_count = model_count + 1
                f.write(
                    'ENDMDL                                               '
                    '                           \n')
                f.write(
                    f'MODEL     {model_count:>4d}                           '
                    f'                                       {chr(10)}')
                num_termini = 0
                new_ser_num = 1
                # reset the numbering when you encounter a new model
        f.write(
            f'{tag}{new_ser_num:>5d} {a["atom_name"]:<4s} {a["res_name"]:<4s}'
            f'{a["chain"]}{a["res_num"]:>4d}{coordinates[n][0]:>12.3f}'
            f'{coordinates[n][1]:>8.3f}{coordinates[n][2]:>8.3f}  1.00  '
            f'0.00          {element:>2s}{chr(10)}')

    # after the very last atom print the final footer

    f.write(
        f'TER   {new_ser_num + 1:>5d}      {info_table[n-1]["res_name"]:<} '
        f'{info_table[n-1]["chain"]}{info_table[n-1]["res_num"]:>4d}         '
        f'                                               {chr(10)}')

    # this follows pdb guidelines for fixed width columns

    f.write(
        'ENDMDL                                                             '
        '             \n')
    f.write(
        'END                                                                '
        '             ')

    f.close()


def ramachandran():
    get_phi_psi()
    alpha_carbons = select(atom_name='CA')
    phis = []
    psis = []
    res_nums = []
    for n, res in enumerate(alpha_carbons):
        phis.append(info_table[res]['phi'])
        psis.append(info_table[res]['psi'])
        res_nums.append(info_table[res]['res_num'])
    return mpl.pyplot.scatter(phis, psis)


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


def random_backbone(n_structs, residue_list_filename, cutoff_distance=0.36,
                    anchor='N', max_tries=20):
    f = open(residue_list_filename, 'r')
    residue_list = f.readlines()

    # trim off the newline and convert to integer for the imported residue
    # numbers whose backbones you wish to randomize
    for ndx, line in enumerate(residue_list):
        residue_list[ndx] = int(line[0:-1])

    # loop the the generation for each of the n_structs you want to generate
    for ndx in range(0, n_structs):
        for num in residue_list:
            try:
                set_phi_psi(random.uniform(0, 2*np.pi),
                            'phi', 'N', res_num=num)

                # if your backbone adjustment resulted in a clash, retry with a
                # different random angle

                tries = 0
                while check_clash(angle_type='phi', anchor='N',
                                  res_num=num) != [] and tries < max_tries:

                    set_phi_psi(random.uniform(0, 2*np.pi), 'phi', 'N',
                                res_num=num)

                    tries += 1
                    if tries < max_tries:
                        print(f'This is the {tries}th try to unclash the phi '
                              f'angle for residue number {num} on model number'
                              f' {ndx}')
                    else:
                        print(f'The maximum number of attempts to unclash '
                              f'model number {ndx} was reached. The model has'
                              f' been discarded.')
            except ValueError:
                print(f'Cant set phi angle for N-terminus or proline at '
                      f'residue {num} of model number {ndx}. Skipping this'
                      f' residue')
                pass
            try:
                tries = 0
                set_phi_psi(random.uniform(0, 2*np.pi),
                            'psi', 'N', res_num=num)
                while check_clash(angle_type='psi', anchor='N',
                                  res_num=num) != [] and tries < 20:
                    set_phi_psi(random.uniform(0, 2*np.pi),
                                'psi', 'N', res_num=num)
                    tries += 1
                    if tries < max_tries:
                        print(f'This is the {tries}th try to unclash the psi '
                              f'angle for residue number {num} on model '
                              f'number {ndx}')
                    else:
                        print(f'The maximum number of attempts to unclash '
                              f'model number {ndx} was reached. The model has'
                              f' been discarded.')
            except ValueError:
                print(f'Cant set phi for C terminus at residue {num} '
                      f'of model number {ndx}. Skipping this residue')

                pass
        ramachandran()
        write_pdb()


def orient(center_sn, axis_sn):
    '''


    Parameters
    ----------
    center_sn : int
        serial number of the atom which is to be centered at [0 0 0]
    allignment_sn : int
        serial number of the atom which will define the new + z axis relative
        to the center atom

    Returns
    -------
    None.

    '''

    global coordinates
    global info_table

    if (len(center_sn) != 1) or (len(axis_sn) != 1) or (center_sn == axis_sn):
        raise ValueError

    else:

        ser_nums = range(0, len(coordinates))
        # the atom to be made [0,0,0]
        center = coordinates[center_sn[0]].copy()
        axis_vector = coordinates[axis_sn[0]]  # atom to be set to [0,0,z]
        #print(f'center is: {center}')
        #print(f'axis vector is: {axis_vector}')
        for n in ser_nums:
            coordinates[n] = coordinates[n] - center

        # find angle between desired axis and yz plane
        x_ax = vector([1, 0, 0])
        y_ax = vector([0, 1, 0])
        z_ax = vector([0, 0, 1])
        #print(f'after recentering, now the axis vector is: {axis_vector}')
        xy = axis_vector.project_onto_normal_plane(z_ax)
        #print(f'the vector projected onto the xy plane is: {xy}')
        yz_ang = y_ax.angle_between(xy)
        #print(f'the angle to the yz plane is: {yz_ang}')

        # rotate everything around the z axis by that angle
        for n in ser_nums:
            coordinates[n] = coordinates[n].rotate_arround(yz_ang, z_ax)

        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # ~~~~~~~~~~~~~~ IMPORTANT NOTE ~~~~~~~~~~~~~~~~~~~~~~
        # for some reason, the direction of the axis of rotation must be
        # reversed for some structures to get the alignment vector to
        # properly be in the yz plane, so we need to check whether the x
        # element after rotation is within floating point error of 0,
        # and rotate in the opposite direction if it is not. that is why
        # the following code block is needed
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if axis_vector[0] < 0.00001 and axis_vector[0] > -0.00001:
            pass
        else:
            # print('\n~~~~~~ERROR: THIS X COORDINATE OF THE ALIGNMENT VECTOR'
            #      ' SHOULD NOW BE 0 (within floating point error) BUT IT '
            #     'ISNT! ~~~~~~~~~~\n\nRetrying by reversing the '
            #    'direction of rotation ...')
            for n in ser_nums:
                coordinates[n] = coordinates[n].rotate_arround(-2*yz_ang, z_ax)
            if axis_vector[0] < 0.00001 and axis_vector[0] > -0.00001:
               # print('ERROR RESOLVED: the alignment vector is now: '
                #      f'{axis_vector}')
                pass
            else:
                # as a last resort (which should never be necessary), try
                # incremental rotations arround the z axis
                attempts = 1
                while axis_vector[0] > 0.00001 and axis_vector[0] < -0.00001 \
                        and attempts < 51:
                    for n in ser_nums:
                        coordinates[n] = coordinates[n].rotate_arround(
                            0.01, z_ax)
                if axis_vector[0] > 0.00001 and axis_vector[0] < -0.00001:
                    raise ValueError

        # now find angle between desired axis and current + z axis
        # and rotate everything by that
        z_ang = axis_vector.angle_between(z_ax)
        #print(f'the angle relative to the +z axis is now: {z_ang}')

        for n, vect in enumerate(coordinates):

            if n != center_sn:
                coordinates[n] = vect.rotate_arround(
                    -1*z_ang, x_ax)


def helix_vector(nterm_res, cterm_res):
    '''

    Parameters
    ----------
    nterm_res : int
        the residue number of the N-terminal most residue
        of the helix in question
    cterm_res : int
        the residue number of the C-terminal most residue
        of the helix in question

    Returns
    -------
        vector
        An R3 vector corresponding to the line of best fit along the CA
        atoms of the helix that indicates the direction the helix is pointing in
    '''
    # get ser nums of all alpha carbons in that range
    c_alpha_ser_nums = select(res_num=(list(range(nterm_res, cterm_res+1))),
                              atom_name='CA')
    n_atoms = len(c_alpha_ser_nums)
    c_alpha_coords = np.zeros(shape=(n_atoms, 3))
    for n, s in enumerate(c_alpha_ser_nums):
        c_alpha_coords[n] = coordinates[s]

    # find the mean of the position vectors of the CAs and recenter
    # them on this
    mean_coord = np.mean(c_alpha_coords, axis=0).view(np.ndarray)
    recentered_coords = (c_alpha_coords - mean_coord).view(np.ndarray)
    # this is mathematically equivalent to principle component analysis
    cov_matrix = np.cov(recentered_coords, rowvar=False, bias=True)
    eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
    max_eigval_ndx = np.argmax(eigenvalues)
    unit_vector = eigenvectors[:, max_eigval_ndx]
    #unit_vector /= np.linalg.norm(unit_vector)
    return unit_vector

# <codecell> Import atom data

# we are making an numpy array filled with the x y and z for each atom in each
# row doing it as a fixed array because this is much more memory efficient


# for w in range(0, 100):
w = 0

fname = f'randomized_S2E/randomized-S2E-oriented-{w}.pdb'
outfile_name = f'randomized_S2E/randomized-S2E-oriented-{w}.pdb'
structure_id = "mode_7"
whole_pdb = parser.get_structure(outfile_name, fname)
n = 0
info_table = []
for model in whole_pdb:
    for chain in model:
        for residue in chain:
            for atom in residue:
                a_tuple = atom.full_id
                a_dict = {
                    'ser_num': n, 'model_name': a_tuple[0],
                    'submodel': a_tuple[1], 'chain': a_tuple[2],
                    'res_num': a_tuple[3][1], 'res_name':
                    residue.get_resname(), 'atom_name': a_tuple[4][0],
                    'phi': float('nan'), 'psi': float('nan')
                }
                # the last 2 are empty fields where phi & psi are assigned

                info_table.append(a_dict)
                n = n + 1
coordinates = np.zeros(shape=(n, 3))
coordinates = coordinates.view(vector)
n = 0
for model in whole_pdb:
    for chain in model:
        for residue in chain:
            for atom in residue:
                coordinates[n] = atom.get_coord()
                n = n + 1

del (a_tuple, n, model, residue, structure_id, chain, atom, a_dict)
#     #print(f'~~~~~~~~~~~~Input file: {fname} ~~~~~~~~~~~~~~')
#     center = select(res_num=23,	atom_name='CA')
#     alignment_atom = select(res_num=13, atom_name='CA')
#     try:
#        orient(center, alignment_atom)
#     except ValueError:
#         print(f'!!!!!!!!!! ERROR: the input pdb file {fname} cannot be '
#               'properly oriented, and no oriented output structure was'
#               ' produced. This is likely due to formatting problems. '
#               'Inspect the input file !!!!!!!!!!!!!')
#         outfile_name = None
#     if outfile_name is not None:
#         #print('structure orientation sucessful. now writing output...')
#         # write_pdb(outfile_name)
#         #print(f'{outfile_name} written sucessfully.\n')
#         pass

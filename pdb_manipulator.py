#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 01:36:49 2023

@author: andrew
"""

# <codecell> Header
import random
import inspect
import os
import matplotlib as mpl
from Bio.PDB.PDBParser import PDBParser
# how we will import the raw text
# from Bio.PDB import PDBIO
import copy
import itertools
import numpy as np
import multiprocessing
import math
from scipy.special import comb


class vector (np.ndarray):
    def __new__(cls, vect):
        if type(vect) == list:
            vect = np.asarray(vect)
        obj = vect.view(cls)
        return obj

    def get_length(self):

        return np.linalg.norm(self)

    def unitize(self):

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
            # Define a tolerance for floating-point errors,
            # if the dot product is too close to -1 then arccos can't
            # be computed directly, but we can just tell it is arbitrary
            # close to pi
            tolerance = 1e-7

            if (c + 1) < tolerance:
                # Vectors are nearly anti-parallel, return pi
                return np.pi
            else:
                # Compute the angle using arccos
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
                    return np.matmul(rot_mat, self)
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


def convert_angle(angle):
    if angle >= 2*np.pi or angle <= -2*np.pi:
        v = angle % 2*np.pi
    else:
        v = angle
    if v > np.pi:
        new_angle = v - 2*np.pi
    elif v < -1*np.pi:
        new_angle = v + 2*np.pi
    else:
        new_angle = v
    return new_angle
# <codecell> Functions


def import_pdb(fname=None):
    # we are making an numpy array filled with the x y and z for each atom in each
    # row doing it as a fixed array because this is much more memory efficient
    # for w in range(0, 100):
    if fname is None:
        fname = 'simple.pdb'
    outfile_name = fname[0:-4]
    structure_id = "mode_7"
    parser = PDBParser(PERMISSIVE=1)
    # permissive=1 allows it to accept incomplete structures without error
    whole_pdb = parser.get_structure(outfile_name, fname)
    print(f'Coordinates for {fname} imported.')
    n = 0
    info_table = []
    coordinates = []
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
    coordinates = np.zeros((3, 3))
    coordinates = np.zeros(shape=(n, 3))
    coordinates = coordinates.view(vector)
    n = 0
    for model in whole_pdb:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coordinates[n] = atom.get_coord()
                    n = n + 1
    return coordinates, info_table


def select(info_table, **kwargs):

    # this function will return the indices in a list of all rows in the info
    # list who match the input arguments
    serial_numbers = []
    info = info_table
    db = kwargs.pop('debug_mode', False)
    # if db:
    #     print(f'debug_mode for calling select is {db}')
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

# !!! this function has grown way too complicated. it should now only allow
# !!! reference to one residue on one chain.
# !!! TO DO: REPLACE THIS WITH A SIMPLE FUNCTION that takes 4 vectors
# !!! and returns a dihedral angle, and another function that calls it
# !!! that determines which 4 vectors those are given a chain, model and
# !!! and residue, and angle type, and relies on a function is_terminus
# !!! to tell if a ser_number is on a residue that is on an n or c terminus.


def isolate_segment(info_table, atom_sn, anchor='N'):
    '''
        NOTE: this does not check that the atom a and B serial numbers correspond
        to directly bonded atoms on the same chain. that responsibility is on the
        calling function to make sure this is correct.

        Parameters
        ----------
        info_table : list of dictionaries
            DESCRIPTION.
        anchor : string
            which side of a chain to hold static, atoms on that side of the atom_sn
            will not have serial numbers returned in the list. Default is 'N'
        atom_sn : Int
            the serial number of an atom directly bonded to another with no 
            cyclization

        Returns
        -------
        segment_ser_num : list of int.

    '''

    # if the atom and anchor are the same then the side chain is included in
    # the segment to be used for rotation
    atom_name = info_table[atom_sn]['atom_name']
    if atom_name in ['N', 'C'] and anchor in ['N', 'C']:

        kwargs = {}
        kwargs['model_name'] = info_table[atom_sn]['model_name']
        kwargs['submodel'] = info_table[atom_sn]['submodel']
        kwargs['chain'] = info_table[atom_sn]['chain']

        # serial numbers on the same chain
        same_chain_atms = select(info_table, **kwargs)

        if atom_name == anchor:

            # this is a special case that demands extra list comprehension
            # due to them being put out of order info_table from what is
            # supposed to be the pdb format
            if atom_name == 'C':
                # get the index of the serial number that has the N atom
                kwargs['res_num'] = info_table[atom_sn]['res_num']
                kwargs['atom_name'] = 'N'
                n = select(info_table, **kwargs)
                kwargs.pop('atom_name')
                same_residue = select(info_table, **kwargs)
                sidechain = [info_table[s]['ser_num'] for s in same_residue
                             if info_table[s]['atom_name'] != 'C']

                # return everything up to that, but also the sidechain atoms
                # and alpha carbon but not the carbonyl
                if n != 0:
                    return same_chain_atms[0:n[0]] + sidechain
                else:
                    return sidechain

            else:
                # get the index of the serial number that has the C atom
                index = next((i for i, sn in enumerate(same_chain_atms)
                              if sn == atom_sn), None)
                # return everything after that

                if index is not None:

                    return same_chain_atms[index+1:]
                else:
                    raise ValueError

        else:
            # NOTE:
            # this method is simpler and saves computation but assumes that
            # atoms are indexed int he correct order by pdb convention
            # and that amide nitrogen is always 2 indexes away from carboxyl
            # in all cases, and the model complete otherwise your pdb must
            # be fixed
            if atom_name == 'C':

                index = next((i for i, sn in enumerate(same_chain_atms)
                              if sn == atom_sn), None)
                if index is not None:
                    return same_chain_atms[index+1:]
                else:
                    raise ValueError
            else:

                index = next((i for i, sn in enumerate(same_chain_atms)
                              if sn == atom_sn), None)
                if index is not None:
                    return same_chain_atms[0:index]
                else:
                    raise ValueError

    else:
        # in the future I will implement a way to directly manipulate side
        # chain rotamers but not now
        print(f'Attempt was made at segment selection on an atom named {atom_name} ,'
              f' with anchor = {anchor}. Currently, selection with atom names'
              'and/or anchors not either "C" or "N" is not implemented.')
        return None

    # depending on the anchor these may or may not be included


def measure_dihedral(point_a: vector, point_b: vector, point_c: vector,
                     point_d: vector):
    '''
    Parameters

    points must be connected in the arrangement seen below
      B---C
     /     \
    A       D

    Given points A, B, C, D as input, each vectors in R3, and for which no
    pair of these points are the same, and for which no combination of
    three of these points are co-linear to each other:
        (this conditional check will be handled by the calling function)

    The vector B to C (that is C - B) forms the axis of rotation.
    the plane P that bisects B to C is normal to it and is the plane of 
    rotation.
    The projection of the vector B to A (B - A) onto this plane forms the 
    reference vector for measuring angle, we will call it A'. The projection 
    of the vector C to D (D - C) onto P will be called D'. The dihedral angle 
    is the angle relative to the vector A' onto the normal plane P such that 
    the negative of this angle is used to compute a rotation matrix that 
    when applied to D results in the D' having an angle of 0 with A'



    Returns

       float64 the dihedral angle in radians between -pi and +pi

    '''
    b_to_a = point_a - point_b
    rotation_axis = point_c - point_b
    c_to_d = point_d - point_c
    a_prime = (b_to_a.project_onto_normal_plane(rotation_axis))
    d_prime = (c_to_d.project_onto_normal_plane(rotation_axis))
    theta = a_prime.angle_between(d_prime)
    # the way to determine whether this angle should be a positive or negative
    # is to see if the cross product of a and d prime is paralell or
    # anti-parallel to the axis of rotation as determined by the dot product
    if np.dot(np.cross(a_prime, d_prime), rotation_axis) < 0:
        theta = -theta

    return theta


def get_phi_psi(coordinates, info_table, **kwargs):
    """
    Update the phi and psi values in the info_table for the specified residues.

    Parameters
    ----------
    res_num : list, optional
        residue numbers for which you want the . The default is None
    debug_mode :
    Returns
    -------
    None.

    """

    if 'debug_mode' in kwargs:
        debug_mode = kwargs['debug_mode']
    else:
        debug_mode = False

    if debug_mode:
        print(f'the kwargs fed into get_phi_psi are {kwargs}')

    selected_residue_ser_nums = []
    if 'chain' in kwargs:
        selected_chains = kwargs.pop('chain', None)
    else:

        # if no chains were specified, make the selected chain the list
        # of all chains, we don't need to make this a kwarg because this is
        # default behaviour of select, but we do need a list we can iterate
        # through
        selected_chains = list({info_table[s]['chain'] for s in
                                range(len(info_table))})

    if debug_mode:
        print(f"the selected chains for get_phi_psi are {selected_chains}")

    # if atoms names were specified in kwargs ignore them and set them to
    # C alpha, otherwise specify c alpha
    kwargs['atom_name'] = 'CA'

    # if 'res_num' in kwargs:
    #     res_num = kwargs['res_num']
    #     if debug_mode:
    #         print('when modified to search for CA on a specific residue the'
    #               f' kwargs fed to select are:\n{kwargs}')
    #         print('get_phi_psi: the serial numbers that match these kwargs'
    #               f' are:\n{select(info_table, **kwargs)}')
    #     if isinstance(res_num, int):  # it must be a list
    #         res_num = [res_num]
    # else:
    #     res_num = list({info_table[s]['res_num'] for s in select(info_table,
    #                                                              **kwargs)})
    #     kwargs['res_num'] = res_num
    #     if debug_mode:
    #         print('when modified to search for CA on a specific residue the'
    #               f' kwargs fed to select are:\n{kwargs}')
    #         print('get_phi_psi: the serial numbers that match these kwargs'
    #               f' are:\n{select(info_table, **kwargs)}')
    #     ser_nums = select(info_table, **kwargs)
    #     for num in ser_nums:
    #         res_num.append(info_table[num]['res_num'])

    # if an explicit list of residues was passed, extract them
    if 'res_num' in kwargs:
        specified_residues = kwargs.pop('res_num', None)

        # if isinstance(specified_residues, int):
        #     specified_residues = [specified_residues]
        #     print(f'specified residues are {specified_residues}')
    else:
        specified_residues = None

    for ch in selected_chains:

        kwargs.pop('res_num', None)
        kwargs['chain'] = ch

        # if no residues were specified, make the selected residues be
        # the list of all residues ON THE CHAIN WE ARE ON IN THE ITERATION
        if specified_residues is None:

            residues_to_query = list({info_table[s]['res_num'] for s in
                                      select(info_table, **kwargs)})
            if debug_mode:
                print("get_phi_psi: no residues were explicitly passed, "
                      f"so the selected residues for chain {ch} are all"
                      f"of them:\n{residues_to_query}")

        # if residues were explicitly specified, iterate over those instead
        else:

            residues_to_query = specified_residues
            if isinstance(residues_to_query, int):
                residues_to_query = [residues_to_query]
            if debug_mode:
                print(
                    f'the residues being queried are:\n{residues_to_query}')

        for n in residues_to_query:

            kwargs['res_num'] = n
            kwargs['atom_name'] = 'N'
            n_ser_num = select(info_table, **kwargs)

            kwargs['atom_name'] = 'CA'
            ca_ser_num = select(info_table, **kwargs)

            kwargs['atom_name'] = 'C'
            c_ser_num = select(info_table, **kwargs)

            # if debug_mode:  # uncomment for extra nitty-gritty
            #     print(f"The current res num is {n}")
            #     print(f'the serial number of the N atom on this residue is'
            #           f" {n_ser_num}")
            #     print('the serial number of the alpha carbon on this'
            #           f"residue is {ca_ser_num}")
            #     print('the serial number of the carboxyl carbon on this'
            #           f"residue is {c_ser_num}")

            ca_to_c_bond_vect = coordinates[c_ser_num] - \
                coordinates[ca_ser_num]
            ca_to_c_bond_vect = ca_to_c_bond_vect.unitize()[0]

            # for some reason subtracting vectors returns a list of
            # vectors with a single element instead of the vector itself making
            # the [0] indexing necessary

            ca_to_n_bond_vect = coordinates[n_ser_num] - \
                coordinates[ca_ser_num]

            ca_to_n_bond_vect = ca_to_n_bond_vect.unitize()[0]

            try:
                kwargs['res_num'] = n - 1

                # if there is no residue that has a number lower than ours we
                # are at the N terminus, for which there is no phi value,
                # otherwise you need the previous carbonyl carbon

                kwargs['atom_name'] = 'C'
                prev_c_ser_num = select(info_table, **kwargs)
                n_to_prev_c_bond_vect = coordinates[prev_c_ser_num] - \
                    coordinates[n_ser_num]
                n_to_prev_c_bond_vect = n_to_prev_c_bond_vect.unitize()[0]

                projection_a = n_to_prev_c_bond_vect.project_onto_normal_plane(
                    ca_to_n_bond_vect)

                projection_b = ca_to_c_bond_vect.project_onto_normal_plane(
                    -1*ca_to_n_bond_vect)

                phi = projection_a.angle_between(projection_b)
                if debug_mode:
                    print(f'Phi angle for residue {n} on chain {ch}: {phi}')
            except ValueError:
                phi = float('nan')

                # if debug_mode:  # uncomment for extra nitty-gritty
                #     print('get_phi_psi: we called for the n-terminus and select'
                #           ' could not find an atom with res_num n - 1')

            try:

                # if there is no residue that has a number higher than ours we
                # are at the c terminus, for which there is no psi value,
                # otherwise you need the next amide nitrogen

                kwargs['res_num'] = n + 1
                kwargs['atom_name'] = 'N'
                next_n_ser_num = select(info_table, **kwargs)
                c_to_next_n_bond_vect = coordinates[next_n_ser_num] - \
                    coordinates[c_ser_num]
                c_to_next_n_bond_vect = c_to_next_n_bond_vect[0]
                projection_a = ca_to_n_bond_vect.project_onto_normal_plane(
                    ca_to_c_bond_vect)
                projection_b = c_to_next_n_bond_vect.project_onto_normal_plane(
                    -1*ca_to_c_bond_vect)
                psi = projection_a.angle_between(projection_b)
            except ValueError:
                psi = float('nan')

                # if debug_mode:  # uncomment for extra nitty-gritty
                #     print('get_phi_psi: we called for the c-terminus and select'
                #           ' could not find an atom with res_num n + 1')

            kwargs['res_num'] = n
            if 'atom_name' in kwargs:
                kwargs.pop('atom_name', None)

            selected_residue_ser_nums = select(info_table, **kwargs)
            # if debug_mode:
            # print(f"line 414: kwargs: {kwargs}")
            # g = {info_table[i]['res_num']
            #       for i in selected_residue_ser_nums}
            # c = {info_table[i]['chain'] for i in selected_residue_ser_nums}
            # print(f'serial numbers for get_phi_psi invoked on residues {g}'
            #       f' of chain(s) {c} are:\n{selected_residue_ser_nums}')

            for i in selected_residue_ser_nums:

                info_table[i]['phi'] = phi
                info_table[i]['psi'] = psi

                if debug_mode:

                    cur_res = str(info_table[i]['res_num']) + \
                        ' ' + info_table[i]['res_name'] + ' @ ' + \
                        info_table[i]['atom_name'] + \
                        ' on chain ' + info_table[i]['chain']
                    print(f'get_phi_psi: phi for {cur_res} is {phi}')
                    print(f'get_phi_psi: psi for {cur_res} is {psi}')


def set_phi_psi(coordinates, info_table, angle, angle_type='phi', anchor='N',
                **kwargs):
    """
    Rotate everything in the coordinates matrix so that the residue given has a
    phi or a psi angle as specified, with the angle specified in degrees by the
    angle parameter. Which end will be held static while the other is rotated
    is specified by the anchor parameter
    """

    # NOTE: this method only works if the coordinates for the amide nitrogen
    # the amide carbon, and the alpha carbon are defined for every residue
    # if they aren't the complete model needs to be built. However checking
    # for this every time would be a major drain on computation time.

    def same_residue(info_table, ser_num):
        # given a serial number, return all serial numbers on the same
        # residue and the same chain
        res = info_table[ser_num]['res_num']
        chain = info_table[ser_num]['chain']
        mod = info_table[ser_num]['model_name']
        smod = info_table[ser_num]['submodel']

        # this will take advantage of the fact that the serial numbers
        # should always be in the order N, CA, C, O, then the rest
        return(select(info_table, res_num=res, chain=chain, model_name=mod,
                      submodel=smod))

    def select_rotating_atoms(info_table, angle_type, anchor, ca_ser_num):
        # given the angle type, the anchor, and the alpha carbon sn, return
        # all the serial numbers that will be subject to rotation matrix
        # multiplication
        pass

    def find_rotation_angle(coordinates, info_table, angle_type, target_angle,
                            ser_nums):
        # given the residue pointed to by the ca_ser_num, and the desired
        # angle type, look up that angle, or call get_phi_psi if needed
        # to get it, and find the angle you must rotate by to end up
        # at the desired angle from the current one. THIS WILL ONLY work
        # properly if the kwargs point to only one residue on one chain

        # amide_nitrogen_sn is ser_nums [0]
        # c_alpha_sn  is ser_nums [1]
        # amide_carbon (carbonyl) is ser_nums[2]

        current_angle = info_table[ser_nums[1]][angle_type]
        if math.isnan(current_angle):

            get_phi_psi(coordinates, info_table, **kwargs)

        # if it is still undefined, you cannot take the angle
        if math.isnan(current_angle):
            return None
        else:
            return convert_angle((target_angle - current_angle))

    def rotate_by_bond(rotating_segment, rot_angl, rot_axis):
        rotated_coords = rotating_segment.rotate_arround(rot_angl, rot_axis)

    if 'debug_mode' in kwargs:
        debug_mode = kwargs['debug_mode']

    else:
        debug_mode = False

    if 'ser_num' in kwargs.keys():
        if isinstance(kwargs['ser_num'], list):
            if len(kwargs['ser_num']) > 1:
                raise ValueError(
                    'Can not set phi psi values by referencing more than'
                    ' one atom serial number as keyword arguments')
            else:
                kwargs['ser_num'] = kwargs['ser_num'][0]
        kwargs['res_num'] = info_table[kwargs['ser_num']]['res_num']
        del kwargs['ser_num']
    # if a serial number for an atom is passed, instead refer to the residue
    # number that that atom is on

    kwargs['atom_name'] = ['CA']
    # if an atom name was specified ignore it
    # and search for 'CA'
    ca_ser_nums = select(info_table, **kwargs)
    print(f'c alpha ser_nums are {ca_ser_nums}')

    for ndx, ca in enumerate(ca_ser_nums):

        print(f'for the alpha carbon serial number {ca}')
        # for every ca on residues where we want to change the angle
        # get all the atoms on the same chain as it
        kwargs.pop('res_num', None)
        kwargs.pop('atom_name', None)
        kwargs['model_name'] = info_table[ca_ser_nums[ndx]]['model_name']
        kwargs['submodel'] = info_table[ca_ser_nums[ndx]]['submodel']
        kwargs['chain'] = info_table[ca_ser_nums[ndx]]['chain']
        print(f'the kwargs fed into select are:\n{kwargs}')
        ser_nums_on_same_chain = select(info_table, **kwargs)
        print(f'the ser nums on the same chain are:\n{ser_nums_on_same_chain}')
        # then change kwargs to specify a specific residue on that chain
        kwargs['res_num'] = info_table[ca_ser_nums[ndx]]['res_num']

        # and then get the amide nitrogen on that same residue
        kwargs['atom_name'] = 'N'
        print(kwargs)
        n_ser_num = select(info_table, **kwargs)[0]

        # and the carboxyl carbon on the same residue
        kwargs['atom_name'] = 'C'
        c_ser_num = select(info_table, **kwargs)[0]

        # save so we can add back later after rotation
        start_ca_position = coordinates[ca].copy()
        print(
            f'for our selected residue, the starting '
            f'C alpha position vector is:\n{start_ca_position}')
        print(find_rotation_angle(coordinates,
              info_table, angle_type, angle, **kwargs))

    if False:
        # ############################     PHI    ###################################
        #     if angle_type == 'phi':
        #         # print('Entering get_phi_psi for phi\n\n')
        #         current_phi = info_table[ca]['phi']

        #         # if the angle hasn't been calculated yet for this residue,
        #         # (try) to do so
        #         if math.isnan(current_phi):

        #             # kwargs.pop('debug_mode', False)

        #             get_phi_psi(coordinates, info_table, **kwargs)
        #             current_phi = info_table[ca]['phi']
        #             if debug_mode:
        #                 print("current phi for residue "
        #                       f"{info_table[ca]['res_num']} of chain "
        #                       f"{info_table[ca]['chain']} at atom "
        #                       f"{info_table[ca]['atom_name']}"
        #                       f" is {info_table[ca]['phi']}")

        #         if info_table[ca]['res_name'] != 'PRO' and \
        #                 math.isnan(current_phi) is False:

        #             # we can't rotate the phi bond on a proline or for 1st residue,
        #             # so check to make sure it isn't one of these
        #             rotate_by = (angle - current_phi)

        #             # find the difference between the desired phi angle and current
        #             # phi angle in
        #             # degrees

        #             if anchor == 'N':

        #                 ###########################################################
        #                 # the last atom on our subset that is on the same chain.
        #                 # Note that this assumes that atoms are listed in order n
        #                 # to c terminus as is pdb convention we need all this to
        #                 # know what the range of serial numbers to apply the
        #                 # rotation to will be. Each time we do a rotation, we only
        #                 # want to rotate residues on the same model, submodel, and
        #                 # chain that are either n terminal (with c anchor, in which
        #                 # case all serial numbers up to ours) or c terminal
        #                 # (with n anchor, in which case starting with ours up to
        #                 # the last on this subset)
        #                 ###########################################################

        #                 end_pt = ser_nums_on_same_chain[-1]
        #                 # if debug_mode:

        #                 #     print(f"same chain atoms are {ser_nums_on_same_chain}")
        #                 #     print("the last atom on the same chain as ours is "
        #                 #           f"{info_table[end_pt]['atom_name']} on "
        #                 #           f"{info_table[end_pt]['res_num']} "
        #                 #           f"{info_table[end_pt]['res_name']} of chain "
        #                 #           f"{info_table[end_pt]['chain']}")
        #                 #     print(f"Before recentering, the coordinates are"
        #                 #           f":\n{coordinates[n_ser_num:end_pt + 1, :]}")

        #                 # center the origin on the alpha carbon only for the
        #                 # rotating subset. note that this assumes that atoms are
        #                 # listed in order n to c terminus as is pdb convention.
        #                 coordinates[n_ser_num:end_pt + 1, :] -= start_ca_position

        #                 # if debug_mode:
        #                 #     print(f"after recentering, the coordinates are:\n"
        #                 #           f"{coordinates[n_ser_num:end_pt + 1, :]}")

        #                 # the n coordinate on the same residue now defines the
        #                 # vector for the bond ca to n
        #                 n_coord = coordinates[n_ser_num]
        #                 CA_to_N_bond_vect = n_coord.unitize()

        #                 if debug_mode:
        #                     print(f"Applying rotation matrix: The rotation angle is"
        #                           f"{angle} - {current_phi} which is {rotate_by}"
        #                           f"\nThe axis of rotation is the normalized "
        #                           f"CA to N bond vector:\n{CA_to_N_bond_vect}")

        #                 # when the vector is unitized this defines our
        #                 # axis of rotation, which is ca to n

        #                 for n, vect in enumerate(coordinates[n_ser_num: end_pt + 1, :]):
        #                     coordinates[(n_ser_num + n), :] = vect.rotate_arround(
        #                         rotate_by, CA_to_N_bond_vect)

        #                 # if debug_mode:
        #                 #     print(f'the result of said rotation is:'
        #                 #           f'\n{coordinates[n_ser_num:(end_pt + 1)]}')

        #                 for n, vect in enumerate(coordinates[n_ser_num: end_pt+1, :]):
        #                     coordinates[(n_ser_num + n), :] = vect + \
        #                         start_ca_position
        #                 # print(
        #                 #     'Adding back the starting c alpha coordinates of the '
        #                 #     f'rotation center:\n{coordinates[n_ser_num:(end_pt +1)]}')

        #             elif anchor == 'C':

        #                 # get the first atom serial number that is on the same
        #                 # model, submodel, and chain
        #                 start_pt = ser_nums_on_same_chain[0]

        #                 # if debug_mode:
        #                 #     print(f"same chain atoms are {ser_nums_on_same_chain}")
        #                 #     print("the first atom on the same chain as ours is "
        #                 #           f"{info_table[start_pt]['atom_name']} on "
        #                 #           f"{info_table[start_pt]['res_num']} "
        #                 #           f"{info_table[start_pt]['res_name']} of chain "
        #                 #           f"{info_table[start_pt]['chain']}")

        #                 # center the origin on the alpha carbon only for the
        #                 # rotating subset. note that this assumes that atoms are
        #                 # listed in order n to c terminus as is pdb convention.
        #                 coordinates[start_pt:ca + 1, :] -= start_ca_position

        #                 # when the vector is reversed and unitized this defines our
        #                 # axis of rotation, which is n to ca
        #                 n_coord = coordinates[n_ser_num]
        #                 CA_to_N_bond_vect = n_coord.unitize()

        #                 for n, vect in enumerate(coordinates[start_pt:ca, :]):
        #                     coordinates[(start_pt + n), :] = vect.rotate_arround(
        #                         rotate_by, CA_to_N_bond_vect)

        #                 coordinates[start_pt:ca + 1] += start_ca_position

        #             else:
        #                 raise ValueError(
        #                     'the 5th positional argument for set_phi_psi must be '
        #                     'either N or C, if none is given N is taken as the '
        #                     'default')

        #             if debug_mode:
        #                 get_phi_psi(coordinates, info_table)
        #                 print('set_phi_psi: after performing rotation, phi for '
        #                       f"residue {info_table[ca]['res_num']} of chain "
        #                       f"{info_table[ca]['chain']} at atom "
        #                       f"{info_table[ca]['atom_name']}"
        #                       f" is {info_table[ca]['phi']*180/np.pi} degrees")

        #         else:
        #             rn = info_table[ca]['res_num']
        #             c = info_table[ca]['chain']
        #             raise ValueError(f'Warning: the phi value for residue {rn}'
        #                              f' on chain {c} can not be set because it is'
        #                              ' either a proline or the n-terminus.'
        #                              ' The request was ignored.')

        # ############################     PSI     ##################################

        #     elif angle_type == 'psi':

        #         # if the angle hasn't been calculated yet for this residue,
        #         # do so
        #         current_psi = info_table[ca]['psi']
        #         if math.isnan(current_psi):

        #             kwargs.pop('debug_mode', False)
        #             # print('Entering get_phi_psi for psi\n\n')
        #             get_phi_psi(coordinates, info_table, **kwargs)
        #             current_psi = info_table[ca]['psi']

        #             if debug_mode:

        #                 print("set_phi_psi: before rotation psi for residue "
        #                       f"{info_table[ca]['res_num']} of chain "
        #                       f"{info_table[ca]['chain']} at atom "
        #                       f"{info_table[ca]['atom_name']}"
        #                       f" is {info_table[ca]['psi']*180/np.pi} degrees")

        #         if math.isnan(current_psi) is False:

        #             # we can't rotate the psi bond for the last residue, so check
        #             # to make sure it isn't

        #             rotate_by = (angle - current_psi)
        #             # find the difference between the desired psi angle and
        #             # current psi angle in degrees

        #             if anchor == 'N':
        #                 end_pt = ser_nums_on_same_chain[-1]
        #                 # the last atom on our subset that is on the same chain.
        #                 # note that this assumes that atoms are listed in order
        #                 # n to c terminus as is pdb convention

        #                 # if debug_mode:
        #                 #     print(f"same chain atoms are {ser_nums_on_same_chain}")
        #                 #     print("the last atom on the same chain as ours is "
        #                 #           f"{info_table[end_pt]['atom_name']} on "
        #                 #           f"{info_table[end_pt]['res_num']} "
        #                 #           f"{info_table[end_pt]['res_name']} of chain "
        #                 #           f"{info_table[end_pt]['chain']}")
        #                 for n, vect in enumerate(coordinates[ca: end_pt + 1, :]):

        #                     coordinates[(ca + n), :] = vect - start_ca_position

        #                 # if debug_mode:
        #                 #     print(f"after recentering, the coordinates are:\n"
        #                 #           f"{coordinates[n_ser_num:end_pt + 1, :]}")

        #                 # center the origin on the alpha carbon only for the
        #                 # rotating subset. note that this assumes that atoms are
        #                 # listed in order n to c terminus as is pdb convention.
        #                 c_coord = coordinates[c_ser_num]
        #                 # the c coordinate on the same residue now defines the
        #                 # vector for the bond ca to c
        #                 CAC_bond_vect = -1*c_coord.unitize()
        #                 # if debug_mode:
        #                 #     print(f"Applying rotation matrix. the angle is"
        #                 #           f"{angle} - {current_psi} which is {rotate_by}"
        #                 #           f"\nThe axis of rotation is the normalized "
        #                 #           f"CAC bond vector:\n{CAC_bond_vect}")

        #                 # when the vector is r unitized this defines our axis of
        #                 # rotation, which is ca to c
        #                 s = coordinates[c_ser_num: end_pt + 1, :]
        #                 for n, vect in enumerate(s):
        #                     coordinates[(c_ser_num + n), :] = \
        #                         vect.rotate_arround(rotate_by, CAC_bond_vect)

        #                 for n, vect in enumerate(coordinates[ca: end_pt + 1, :]):
        #                     coordinates[(ca + n), :] = vect + start_ca_position

        #             elif anchor == 'C':

        #                 start_pt = ser_nums_on_same_chain[0]
        #                 # if debug_mode:
        #                 #     print(f"same chain atoms are {ser_nums_on_same_chain}")
        #                 #     print("the first atom on the same chain as ours is "
        #                 #           f"{info_table[start_pt]['atom_name']} on "
        #                 #           f"{info_table[start_pt]['res_num']} "
        #                 #           f"{info_table[start_pt]['res_name']} of chain "
        #                 #           f"{info_table[start_pt]['chain']}")

        #                 # get the first atom serial number that is on the same
        #                 # model, submodel, and chain

        #                 coordinates[start_pt:c_ser_num+1, :] -= start_ca_position
        #                 # center the origin on the alpha carbon only for the
        #                 # rotating subset. note that this assumes that atoms are
        #                 # listed in order n to c terminus as is pdb convention.
        #                 c_coord = coordinates[c_ser_num]
        #                 CAC_bond_vect = c_coord.unitize()
        #                 # when the vector is unitized this defines our axis of
        #                 # rotation, which is ca to n
        #                 for n, vect in enumerate(coordinates[start_pt:
        #                                                      c_ser_num + 1, :]):
        #                     coordinates[(start_pt + n), :] = \
        #                         vect.rotate_arround(rotate_by, CAC_bond_vect)
        #                 coordinates[start_pt:c_ser_num + 1] += start_ca_position

        #             else:
        #                 raise ValueError(
        #                     'the 5th positional argument for set_phi_psi must be '
        #                     'either N or C, if none is given N is taken as the '
        #                     'default')
        #             if debug_mode:
        #                 get_phi_psi(coordinates, info_table)
        #                 print('set_phi_psi: after performing rotation, psi for '
        #                       f"residue {info_table[ca]['res_num']} of chain "
        #                       f"{info_table[ca]['chain']} at atom "
        #                       f"{info_table[ca]['atom_name']}"
        #                       f" is {info_table[ca]['psi']*180/np.pi} degrees")

        #         else:
        #             rn = info_table[ca]['res_num']
        #             c = info_table[ca]['chain']
        #             raise ValueError(f'Warning: the psi value for residue {rn}'
        #                              f' on chain {c} can not be set because it is'
        #                              ' at the c-terminus. The request was'
        #                              ' ignored.')
        #     else:
        #         raise ValueError(
        #             f'Warning: the dihedral angle type "{repr(angle_type)}" is '
        #             'not a valid selection. It must be either phi or psi. If '
        #             'none is given, phi is default')
        # get_phi_psi(coordinates, info_table, **kwargs)
        pass


def write_pdb(coordinates, info_table, outfile=None):
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

    # if the filename is already taken, augment the number at the end of the
    # filename

    try:
        f = open(f'{model_name}.pdb', 'x')
    except FileExistsError:
        count = 0
        while True:
            try:
                f = open(f'{model_name}_{count}.pdb', 'x')
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


def ramachandran(info_table):
    get_phi_psi()
    alpha_carbons = select(info_table, atom_name='CA')
    phis = []
    psis = []
    res_nums = []
    for n, res in enumerate(alpha_carbons):
        phis.append(info_table[res]['phi'])
        psis.append(info_table[res]['psi'])
        res_nums.append(info_table[res]['res_num'])
    return mpl.pyplot.scatter(phis, psis)


def is_clashing(coordinates, search_set_indices, threshold):
    '''
    this will perform distance measurements for all coordinates contained in
    each of the lists inside search_set_indices against all coordinates in
    the other sublists (not those in the same sublist). If it encounters a
    distance less than the threshold, it will stop execution and return the
    threshold pair of indices.

    Parameters
    ----------
    coordinates : array of vectors
        the entire set of mulecular coordinates
    search_set_indices : list
        a list containing at least 2 lists  which each contain the indices of
        a set of coordinates to checked against all members of all other sets
        in the top list.


    Returns
    -------
    results : TYPE
        DESCRIPTION.

    '''

    # this inner function will divide up the pairwise search by a divide and
    # conquer multiprocessing algorithm to make it faster. if any of the
    # subprocesses halt, each is terminated.
    def parallel_search(ndx_list_a, ndx_list_b, threshold, terminate_event):
        results = multiprocessing.Queue()

        def pairwise_search():
            # print('top of pairwise search')
            for i in range(len(ndx_list_a)):
                for j in range(len(ndx_list_b)):
                    # print(
                    #     f'comparing against the {j} element of the'
                    #     'second group')
                    if terminate_event.is_set():
                        # print('aborting')
                        # if this or any other subprocess sent term signal
                        # just abort
                        return None
                    # print(ndx_list_a[i], ndx_list_b[j])
                    dist = coordinates[ndx_list_a[i]].distance_between(
                        coordinates[ndx_list_b[j]])
                    if dist < threshold:
                        # print('clash detected')
                        # Put the result into the queue created by the
                        # calling function
                        results.put([ndx_list_a[i], ndx_list_b[j]])
                        # then set the event variable
                        # to tell each subprocess to stop
                        terminate_event.set()
                        return None

        with multiprocessing.Pool() as pool:
            for i in range(len(ndx_list_a)):
                # print(
                #     f'in the {i}th element of the first coord set'
                #     ' starting async')
                pool.apply_async(
                    pairwise_search())
        pool.close()
        pool.join()
        return results

    if threshold is not None:

        terminate_event = multiprocessing.Event()

        # list all possible combinations of coordinate set pairs to be searched
        # against each other without redundancy

        # for the unique pairs search every element of the first against
        # every element of the second
        for pair in itertools.combinations(search_set_indices, 2):
            # print(pair)
            halt_indices = parallel_search(pair[0], pair[1], threshold,
                                           terminate_event)
            if not halt_indices.empty():
                return halt_indices.get()
        # if the function never picks up a distance below the threshold then
        # return None
        return None


def check_clash_all(coordinates, info_table, cutoff_distance=0.36):
    """
    checks every atom against every other atom not directly bonded to it for
    clashes. This is highly computationally expensive and to be avoided
    if possible.

    Parameters
    ----------
    cutoff_distance : TYPE, optional
        DESCRIPTION. The default is 0.36.

    Returns
    -------
    None.

    """

    pass


def check_internal_clash(coordinates, info_table, cutoff_distance=0.36,
                         angle_type='phi', anchor='N', **kwargs):
    """
    Check if a rotation arround the phi or psi bond (default: phi) of the
    input residue has resulted in atoms closer than the cutoff radius (in nm)
    Return any atom serial numbers that clash. Invoke after performing
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

    same_chain_atms = select(info_table, **kwargs)
    kwargs['res_num'] = res

    if angle_type == 'phi':

        kwargs['atom_name'] = 'N'   # get amide nitrogen
        amide_nitrogen = select(info_table, **kwargs)
        if len(amide_nitrogen) > 1:
            raise ValueError

        # split the chain onto the atoms on one side of the bond and those
        # on the other side
        n_term_frag = same_chain_atms[0: amide_nitrogen[0]+1]
        c_term_frag = same_chain_atms[amide_nitrogen[0]+1:]
        try:
            kwargs['atom_name'] = 'HN'
            amide_hydrogen = select(info_table, **kwargs)
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
        amide_carbon = select(info_table, **kwargs)

        if len(amide_carbon) > 1:
            raise ValueError

        kwargs['atom_name'] = 'N'
        kwargs['res_num'] = res + 1
        next_amide_nitro = select(info_table, **kwargs)
        kwargs['res_num'] -= 1
        n_term_frag = same_chain_atms[0: next_amide_nitro[0]]
        c_term_frag = same_chain_atms[next_amide_nitro[0]:]

        # remove the amide carbon and oxygen from n terminal side and add to
        # c terminal side, again, because they are discontiguous
        try:

            kwargs['atom_name'] = 'C'
            amide_carbon = select(info_table, **kwargs)
            c_term_frag.append(amide_carbon[0])
            n_term_frag.remove(amide_carbon[0])
            kwargs['atom_name'] = 'O'
            amide_oxygen = select(info_table, **kwargs)
            c_term_frag.append(amide_oxygen[0])
            n_term_frag.remove(amide_oxygen[0])

        except ValueError:
            pass

    else:    # if input was invalid
        raise ValueError
        print(
            f'{angle_type} is not a valid input for parameter angle_type. '
            f'Enter either phi or psi')

    clash_sn = is_clashing(
        coordinates, [n_term_frag, c_term_frag], cutoff_distance)
    if clash_sn is not None:
        return clash_sn
    else:
        return []
    # clash_atom_ser_nums = []
    # search_list = []
    # for r in info_table:
    #     search_list.append(r['ser_num'])

    # if anchor == 'N':

    #     # remove the c-terminal fragment from the searchlist
    #     # so it doesn't search against itself
    #     for atom in c_term_frag:
    #         search_list.remove(atom)
    #     # then search each atom in the fragment against every atom not a part
    #     # of that fragment
    #     for atom in c_term_frag:
    #         for query in search_list:
    #             dist = coordinates[atom].distance_between(coordinates[query])
    #             if dist < cutoff_distance:
    #                 clash_atom_ser_nums.append([atom, query])

    # elif anchor == 'C':

    #     # remove the c-terminal fragment from searchlist
    #     # so it doesn't search against itself
    #     for atom in n_term_frag:
    #         search_list.remove(atom)
    #     # then search each atom in the fragment against every atom not a part
    #     # of that fragment
    #     for atom in n_term_frag:
    #         for query in search_list:
    #             dist = coordinates[atom].distance_between(coordinates[query])
    #             if dist < cutoff_distance:
    #                 clash_atom_ser_nums.append([atom, query])

    # else:
    #     raise ValueError
    #     print(
    #         f'{anchor} is not a valid input for parameter anchor. Enter '
    #         f'either N or C')
    # return clash_atom_ser_nums


def check_external_clash(coordinates, info_table, ser_nums_a, ser_nums_b, cutoff_distance=0.36):

 # ######################## END GOAL #####################################
 # check first heuristically (since a positive result is definitive), then
 # pair-wise since that is definitive for a negative result in order to be garuanteed
 # to properly identify that this function catches all potential clashes for the chains specified
 # that occur between different chains
 # when invoked with check_internal_class it garuantees that in the case of
 # chains related by radial symmetry, that any one of those in a set checked this
 # way against any one other can be garuanteed to come from a set free of clashes between all members of the set
 # this means you only need to call this functions once to ensure the set is
 # clash-free. note this should change the output type from a list of clashing indeces
 # to a boolean statement as to whether any such classes exist.
 # look up a heuristic algorith such as Bounding Volume Hierarchy (BVH):
 # Distance Thresholding, Spatial Partitioning, and Hashing to implement this.
 # because pairwise comparison for all atoms between all other atoms for symmetrical chains
 # is otherwise what would have to be done
 #
 # For now:
 # just implement the definitive test, since this will work just as well, although less efficiently

    is_clashing = False

    def heuristic_test():
        pass

    def pairwise_test():
        for sa in ser_nums_a:
            for sb in ser_nums_b:
                if (coordinates[sa] - coordinates[sb]).get_length < cutoff_distance:
                    is_clashing = True
                    return is_clashing

    # if heuristic_test([chain_a, chain_b])

        pass


def random_backbone(coordinates, info_table, n_structs, residue_list_filename,
                    cutoff_distance=0.36, chainlist=None,
                    symmetry_groups='all', max_tries=20):

    f = open(residue_list_filename, 'r')
    residue_list = f.readlines()

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!! this processes correctly and sets up the anchored and free segments
# !!!! right but doesn't seem to change the phi and psi angles
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # !!! TO DO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # change it to accept input file format that include the chain along with
    # the residue designations for fixed and anchored residues, so that it can
    # be determined which ones are symmetrically related to which. This will
    # allow heterogenous set of chains while still allowing for assumptions
    # about symmetry to be used to simplify clash checking.
    # !!!

    # symmetry_groups will be a list of lists each containing Chain letters of
    # chains which are related to each other by symmetry this will be used to expedite
    # clash checking. if none are given, all chains will be treated as identical in
    # sequence qand related by symmetry together in one symmetry group, else if the value None is
    # passed, or an empty list is passed, they will all be treated as different

    # residue numbers after the header 'anchored:' in the residue file will
    # be put in a separate list, these residues cannot move at all
    anchor_header_line = None

    # residue numbers after the header 'fixed:' in the residue file will be in
    # a different group. these can move as a result of adjusting other
    # residues, but cannot have their own internal coordinates adjusted
    fixed_header_line = None

    # trim off the newline and convert to integer for the imported residue
    # numbers whose backbones you wish to randomize
    for ndx, line in enumerate(residue_list):
        residue_list[ndx] = line.strip()

    # get rid of any blank lines
    residue_list = [element for element in residue_list if element != '']
    for ndx, line in enumerate(residue_list):
        if not line.isdigit():
            if line == 'anchored:':
                if anchor_header_line is None:
                    anchor_header_line = ndx
                else:
                    raise ValueError('Only one contiguous group of residues'
                                     ' is permitted to be anchored')
            elif line == 'fixed:':
                fixed_header_line = ndx
            else:
                raise ValueError(f'Error: unrecognized header {line} in the'
                                 f'file {residue_list_filename}')
    for ndx, line in enumerate(residue_list):
        # print(line)
        if ndx != anchor_header_line and ndx != fixed_header_line:
            residue_list[ndx] = int(line)

    # conditional handling for formatting of the residue file
    if anchor_header_line is not None:
        if fixed_header_line is not None:
            if fixed_header_line != len(residue_list)-1:
                if anchor_header_line != len(residue_list)-1 and \
                        fixed_header_line > anchor_header_line:
                    anchored_residues = residue_list[anchor_header_line +
                                                     1:fixed_header_line]
                    fixed_residues = residue_list[fixed_header_line+1:]
                elif anchor_header_line != len(residue_list)-1 and \
                        fixed_header_line < anchor_header_line:
                    anchored_residues = residue_list[anchor_header_line+1:]
                    fixed_residues = []
            elif anchor_header_line != len(residue_list)-1:
                anchored_residues = residue_list[anchor_header_line+1:]
                fixed_residues = []
            else:
                anchored_residues = []
                fixed_residues = []
        elif anchor_header_line != len(residue_list)-1:
            anchored_residues = residue_list[anchor_header_line+1:]
            fixed_residues = []
        else:
            anchored_residues = []
            fixed_residues = []
    elif fixed_header_line is not None and fixed_header_line != \
            len(residue_list)-1:
        fixed_residues = residue_list[fixed_header_line+1:]
        anchored_residues = []
    else:
        fixed_residues = []
        anchored_residues = []

    # anchored residues only make sense if they are all contiguous in sequence
    # in the same chain. check that this is the case
    if len(anchored_residues) > 1:
        for ndx, res in enumerate(anchored_residues[:-1]):
            if res != anchored_residues[ndx + 1] - 1:
                raise ValueError('All residues in the anchored sequence must'
                                 ' be contiguous')
    print(f'anchored_residues are {anchored_residues}')
    print(f'fixed_residues are {fixed_residues}')

    # if no chains were passed, use all of them
    if chainlist is None:
        chainlist = set()
        for c in info_table:
            chainlist.add(c['chain'])
        chainlist = list(chainlist)
    else:
        # eliminate redundancy
        chainlist = set(chainlist)
        chainlist = list(chainlist)

    # a list of lists each containing the serial numbers from selected chains
    # in chain_list
    sn_by_chain = []
    for c in chainlist:
        sn = select(info_table, chain=c)
        sn_by_chain.append(sn)
        # print(f'sn_by_chain is {sn_by_chain}')

    free_residues = []
    for sn_set in sn_by_chain:
        free_residues.append({info_table[s]['res_num'] for s in sn_set})

    # whittle away the anchored residues from the list to iterate
    # over when randomizing angles
    free_residues = [[s for s in sn_set if s not in anchored_residues]
                     for sn_set in free_residues]

    print(f'after set comprehension, the free residues are\n{free_residues}')

    # if an anchor section exists, all residues before it are anchored at
    # C term and all residues after are anchored at N term
    n_anchor = None
    c_anchor = None

    # !!! NOTE: this currently will only work if the anchored residues are the
    # same for every chain. If anchors and fixed residues will be different
    # for different chains, this function should be called separately for
    # each instance passing only one chain each time
    # !!!

    if anchored_residues:
        # !!! in the future change this to a list containing minima for each
        # !!! chain segment of anchored residues
        n_anchor = min(anchored_residues)
        print(f'the n anchor is {n_anchor}')
        c_anchor = max(anchored_residues)
        print(f'the c anchor is {c_anchor}')
    # if an anchored segment exists, the free residues are split into those on
    # the n_terminal side of it (c_anchored) or the c terminal side
    # (n_anchored)

    c_anchored_segments = []
    n_anchored_segments = []
    print(f'there are {len(sn_by_chain)} chains being randomized')

    # for each chain that needs to be checked internally, add a list of
    # its n terminal and c terminal segments serial numbers as a
    # sublist to keep track of associated chains

    for i in range(len(sn_by_chain)):

        # if there are no anchors treat that whole chain of free residues as
        # as being N anchored at the first residue
        if n_anchor is None and c_anchor is None:
            n_anchored_segments.append(free_residues[i])
        else:
            if n_anchor is not None:

                seg = [res for res in free_residues[i] if res < n_anchor]

                # then put it in reverse order so that rotation iterates from C
                # towards n terminus

                print(f' the {i}th set of c_anchored segments is {seg}')
                seg.sort(reverse=True)
                c_anchored_segments.append(seg)
            if c_anchor is not None:
                # repeat this for n anchored segments
                seg = [res for res in free_residues[i] if res > c_anchor]
                print(seg)
                seg.sort(reverse=False)
                n_anchored_segments.append(seg)
        # print(f'the nterm seg is\n{c_anchored_segments}\nthe c term seg '
             # f'is\n{n_anchored_segments}')
    print(f'there are {len(c_anchored_segments)} c anchored segments')
    print(f'there are {len(c_anchored_segments)} n anchored segments')
    print(f'there are {len(chainlist)} chain segments')
    if len(c_anchored_segments) != len(n_anchored_segments) and \
            len(free_residues) != len(chainlist):
        raise ValueError('the number of lists in c_anchored_segments, '
                         'n_anchored_segments and the number of items in '
                         'chainlist must all be equal')

    print(f'the chain in this selection are\n{chainlist}\n')
    print(f'the C-anchored segments for each chain are {c_anchored_segments}')
    print(f'the N-anchored segments for each chain are {n_anchored_segments}')
    # loop the generation for each of the n_structs you want to generate
    for ndx in range(0, n_structs):
        for i in range(len(c_anchored_segments)):

            if c_anchored_segments[i] != []:

                for num in c_anchored_segments[i]:
                    try:
                        set_phi_psi(coordinates, info_table,
                                    random.uniform(0, 2*np.pi),
                                    angle_type='phi', anchor='C', res_num=num,
                                    chain=chainlist[i])

                        # if your backbone adjustment resulted in a clash,
                        # retry with a different random angle

                        tries = 0
                        while (check_internal_clash(
                                coordinates, info_table, cutoff_distance,
                                angle_type='phi', anchor='C', res_num=num,
                                chain=chainlist[i]) != []
                               and tries < max_tries):

                            set_phi_psi(coordinates, info_table,
                                        random.uniform(0, 2*np.pi),
                                        angle_type='phi', anchor='C',
                                        res_num=num, chain=chainlist[i])

                            tries += 1
                            if tries < max_tries:
                                print(f'This is the {tries}th try to unclash '
                                      f'the phi angle for residue number {num}'
                                      f' on model number {ndx}')
                            else:
                                print(f'The maximum number of attempts to '
                                      f'unclash model number {ndx} was '
                                      'reached. The model has been discarded.')
                    except ValueError:
                        print(f'Cant set phi angle for N-terminus or proline'
                              f' at residue {num} of model number {ndx}. '
                              'Skipping this residue')
                        pass
                    try:
                        tries = 0
                        set_phi_psi(coordinates, info_table,
                                    random.uniform(0, 2*np.pi),
                                    angle_type='psi', anchor='C',
                                    res_num=num, chain=chainlist[i])
                        while (check_internal_clash(coordinates, info_table,
                                                    angle_type='psi',
                                                    anchor='C',
                                                    res_num=num,
                                                    chain=chainlist[i]) != []
                                and tries < 20):
                            set_phi_psi(coordinates, info_table,
                                        random.uniform(0, 2*np.pi),
                                        angle_type='psi', anchor='C',
                                        res_num=num, chain=chainlist[i])
                            tries += 1
                            if tries < max_tries:
                                print(f'This is the {tries}th try to unclash '
                                      f'the psi angle for residue number '
                                      f'{num} on model number {ndx}')
                            else:
                                print(f'The maximum number of attempts to '
                                      f'unclash model number {ndx} was '
                                      'reached. The model has been discarded.')
                    except ValueError:
                        print(f'Cant set psi for C terminus at residue {num} '
                              f'of model number {ndx}. Skipping this residue.')
                        pass
            if n_anchored_segments[i] != []:
                for num in n_anchored_segments:
                    try:
                        set_phi_psi(coordinates, info_table,
                                    random.uniform(0, 2*np.pi),
                                    angle_type='phi', anchor='N',
                                    res_num=num, chain=chainlist[i])

                    # if your backbone adjustment resulted in a clash, retry
                    # with a different random angle

                        tries = 0
                        while (check_internal_clash(coordinates, info_table,
                                                    angle_type='phi',
                                                    anchor='N',
                                                    res_num=num,
                                                    chain=chainlist[i]) != []
                               and tries < 20):

                            set_phi_psi(coordinates, info_table,
                                        random.uniform(0, 2*np.pi),
                                        angle_type='phi', anchor='N',
                                        res_num=num, chain=chainlist[i])

                            tries += 1
                            if tries < max_tries:
                                print(f'This is the {tries}th try to unclash'
                                      f'the phi angle for residue number {num}'
                                      f' on model number {ndx}')
                            else:
                                print(f'The maximum number of attempts to '
                                      f'unclash model number {ndx} was '
                                      'reached. The model has been discarded.')
                    except ValueError:
                        print(f'Cant set phi angle for N-terminus or proline '
                              f'at residue {num} of model number {ndx}. '
                              'Skipping this residue.')
                        pass
                    try:
                        tries = 0
                        set_phi_psi(coordinates, info_table,
                                    random.uniform(0, 2*np.pi),
                                    angle_type='psi', anchor='N',
                                    res_num=num, chain=chainlist[i])
                        while (check_internal_clash(coordinates, info_table,
                                                    angle_type='psi',
                                                    anchor='N',
                                                    res_num=num,
                                                    chain=chainlist[i]) != []
                               and tries < 20):
                            set_phi_psi(coordinates, info_table,
                                        random.uniform(0, 2*np.pi),
                                        angle_type='psi', anchor='N',
                                        res_num=num, chain=chainlist[i])
                            tries += 1
                            if tries < max_tries:
                                print(f'This is the {tries}th try to unclash'
                                      f' the psi angle for residue number '
                                      f'{num} on model number {ndx}')
                            else:
                                print(f'The maximum number of attempts to '
                                      f'unclash model number {ndx} was '
                                      'reached. The model has been discarded.')
                    except ValueError:
                        print(f'Cant set phi for C terminus at residue {num} '
                              f'of model number {ndx}. Skipping this residue')
                        pass
            if free_residues != []:
                pass
        # # # ramachandran()
        # # write_pdb()


def orient(coordinates, info_table, center_sn, axis_sn):
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

    if (len(center_sn) != 1) or (len(axis_sn) != 1) or (center_sn == axis_sn):
        raise ValueError

    else:

        ser_nums = range(0, len(coordinates))
        # the atom to be made [0,0,0]
        center = coordinates[center_sn[0]].copy()
        axis_vector = coordinates[axis_sn[0]]  # atom to be set to [0,0,z]
        for n in ser_nums:
            coordinates[n] = coordinates[n] - center

        # find angle between desired axis and yz plane
        x_ax = vector([1, 0, 0])
        y_ax = vector([0, 1, 0])
        z_ax = vector([0, 0, 1])
        xy = axis_vector.project_onto_normal_plane(z_ax)
        yz_ang = y_ax.angle_between(xy)

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
        # print(f'the angle relative to the +z axis is now: {z_ang}')

        for n, vect in enumerate(coordinates):

            if n != center_sn:
                coordinates[n] = vect.rotate_arround(
                    -1*z_ang, x_ax)


def helix_vector(coordinates, info_table, nterm_res, cterm_res):
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
    c_alpha_ser_nums = select(info_table, res_num=(list(range(nterm_res, cterm_res+1))),
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
    # unit_vector /= np.linalg.norm(unit_vector)
    return unit_vector


def import_from_vmd(filename=None):

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #  UNDER CONSTRUCTION
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if filename is None:
        pass
    else:
        try:
            w = 1

            filename = f'cg-trajectories/S2E-processed-cg-input-{w}.pdb'
            # outfile_name = f'randomized_S2E/randomized-S2E-oriented-{w}.pdb'
            whole_pdb = PDBParser.get_structure(w, filename)
            n = 0
            global info_table
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

        except:
            pass


def clone_chain(coordinates, info_table, chain_list=None):
    '''


    Parameters
    ----------
    coordinates : Array of float 64
        the atomic coordinates of all atoms
    info_table : list
        list of dictionaries containing data for each coordinate
    chain_list : list of strings, optional
        Which chains should be cloned. If none given all will be selected.

    Returns
    -------
    coordinates : Array of float 64
        the updated atomic coordinates of all atoms after cloning chains.
    info_table : list
        the updated list of dictionaries containing data for each coordinate
        after cloning chains
    new_chain_list : list
        list of the chain letters created by this function

    '''
    # getting started, we need to know what chains we want copied, and which
    # ones already exist
    existing_chains = list({entry['chain'] for entry in info_table})
    if chain_list is None:
        # if no chain was specified, select all chains
        chain_list = existing_chains
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    availible_letters = [
        letter for letter in alphabet if letter not in existing_chains]
    new_chain_list = []
    for chain in chain_list:

        # assign a new letter for the new chain and remove tha letter from the
        # list of availible ones
        new_chain = min(availible_letters)
        new_chain_list.append(new_chain)
        availible_letters = [
            letter for letter in availible_letters if letter is not new_chain]
        mother_chain_sn = select(info_table, chain=chain)
        smallest_sn = min(mother_chain_sn)
        last_ndx = len(info_table)
        daughter_chain_sn = [x - smallest_sn +
                             last_ndx for x in mother_chain_sn]
        # make temp array that is slice of coordinates array where chain is stored
        # then concetenate the two
        new_coords = coordinates[mother_chain_sn[0]: mother_chain_sn[-1]+1]
        coordinates = np.concatenate((coordinates, new_coords)).view(vector)
        # print(coordinates)
        # print(info_table[0])
        info_table.extend(copy.deepcopy(
            info_table[mother_chain_sn[0]: mother_chain_sn[-1]+1]))
        for i, sn in enumerate(daughter_chain_sn):
            info_table[sn]['ser_num'] = sn
            info_table[sn]['chain'] = new_chain
    return coordinates, info_table, new_chain_list


def axial_symmetry(coordinates, info_table, chains, multiplicity, cofr=vector([0.0, 0.0, 0.0]), rot_axis=vector([0.0, 0.0, 1.0])):
    '''


        Parameters
        ----------
        chains : LIST
            list containing strings indicating chains that will be multimerized.
        multiplicity : INT16
            positive integer describing how many subunits in the arrangement
        cofr : VECTOR, optional
            the center of rotation which the transformation will be applied arround.
            The default is vector([0.0, 0.0, 0.0]).
        rot_axis : VECTOR, optional
            a vector that serves as axis of rotation after recentering on cofr.
            The default is vector([0.0, 0.0, 1.0]).

        Returns
        -------
        None.

        '''
    try:
        rot_axis = rot_axis.unitize()
    except ZeroDivisionError:
        print('The axis of rotation can not be the zero vector')

    if multiplicity < 2 or type(multiplicity) is not int:
        print('multiplicity for axial symmetry must be an integer 2 or greater')
        raise ValueError

    # a list containing lists of chains which are identical
    chain_groups = []
    if chains:
        for c in chains:
            group = []
            group.append(c)
            try:
                for x in range(0, (multiplicity - 1)):
                    coordinates, info_table, new_chains = clone_chain(
                        coordinates, info_table, (c)[0])
                    group.append(new_chains[0])
                chain_groups.append(group)
                print(f'the group is {group}')
            except ValueError:
                print(
                    f'chain {c} does not exist or can not be copied. skipping.')
        for group in chain_groups:
            for n, chain_id in enumerate(group):
                angle = 2*np.pi/multiplicity*n
                current_chain = select(info_table, chain=chain_id)
                for sn in current_chain:
                    coordinates[sn] = coordinates[sn].rotate_arround(
                        angle, vector([0, 0, 1]))

    else:
        print('No chains were provided')
        raise ValueError
    return coordinates, info_table

######################## for aligning multiple struture files en batch #######
#  center = select(info_table, res_num=23,	atom_name='CA')
#   alignment_atom = select(info_table, res_num=13, atom_name='CA')
#    try:
#         orient(center, alignment_atom)
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


# axial_symmetry(['A'], 5, vector([-5.2, -9, 0])) this produces excelent overlay
# over pdb 7k3g 1st submodel when applied to randomized-S2E-oriented-0.pdb

coordinates, info_table = import_pdb('simple.pdb')

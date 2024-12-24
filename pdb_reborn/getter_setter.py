#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 00:52:03 2024

@author: andrew
"""

import random
# how we will import the raw text
# from Bio.PDB import PDBIO
import copy
import itertools
import numpy as np
import math
from io_funcs import import_pdb, plot_model
from coord_funcs import cart2cylind, convert_angle
from vector_class import vector

def select(info_table, **kwargs):

    # this function will return the indices in a list of all rows in the info
    # list who match the input arguments
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


def get_dihedral_backbones(info_table, ca_ser_nums):

    # given a LIST of alpha carbon serial numbers, return sets of the serial
    # numbers for surrounding backbone atoms needed to get/set dihedrals

    if isinstance(ca_ser_nums, int):
        ca_ser_nums = [ca_ser_nums]
    bb_ser_nums = []
    kwargs = {}  # what we feed into select
    for ndx, sn in enumerate(ca_ser_nums):
        kwargs['res_num'] = info_table[sn]['res_num'] - 1
        kwargs['chain'] = info_table[sn]['chain']
        kwargs['model_name'] = info_table[sn]['model_name']
        try:
            # try to assign the previous carboxyl
            kwargs['atom_name'] = 'C'
            prev_c = select(info_table, **kwargs)
            prev_c = prev_c[0]
        except ValueError:
            prev_c = None
        amide = sn - 1  # amide should always be behind the C alpha
        kwargs['res_num'] += 1
        carboxyl = select(info_table, **kwargs)
        carboxyl = carboxyl[0]
        kwargs['atom_name'] = 'N'
        kwargs['res_num'] += 1
        try:
            # try to assign the next amide
            next_amide = select(info_table, **kwargs)
            next_amide = next_amide[0]
        except ValueError:
            next_amide = None
        bb_ser_nums.append([prev_c, amide, sn, carboxyl, next_amide])
    return bb_ser_nums


def isolate_segment(info_table, atom_sn, anchor='N'):
    '''
        DESCRIPTION.
        NOTE: this does not check that the atom a and B serial numbers correspond
        to directly bonded atoms on the same chain. that responsibility is on the
        calling function to make sure this is correct.

        Parameters
        ----------
        info_table : list of dictionaries

        anchor : string
            which side of a chain to hold static, atoms on that side of the
            atom_sn will not have serial numbers returned in the list.
            Valid input is 'N' or 'C'. Default is 'N'.

        atom_sn : Int
            the serial number of an atom directly bonded to another with no
            cyclization i.e. it is not a member of a ring

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
        kwargs['res_num'] = info_table[atom_sn]['res_num']
        k = kwargs['res_num']
        cur_res = [sn for sn in same_chain_atms if info_table[sn]['res_num'] == k]
        same_chain_cur_res_ndx = [same_chain_atms.index(
            sn) for sn in same_chain_atms if info_table[sn]['res_num'] == k]

        if atom_name == anchor:

            # this is a special case that demands extra list comprehension
            # due to them being put out of order info_table from what is
            # supposed to be the pdb format
            if atom_name == 'C':
                # get the index of the serial number that has the N atom

                kwargs['atom_name'] = 'N'
                n = select(info_table, **kwargs)
                kwargs.pop('atom_name')

                sidechain = [info_table[s]['ser_num'] for s in cur_res
                             if info_table[s]['atom_name'] not in ['C', 'O']]
                # return everything up to that, but also the sidechain atoms
                # and alpha carbon but not the carbonyl or oxygen

                return same_chain_atms[0:n[0]] + sidechain

            else:
                # get the index of the serial number that has the next amide
                next_amide = max(cur_res) + 1
                # number right after last in curr res should be amide Nitro

                # get everything on current res except amide n and hydro
                amide_excluded = [sn for sn in cur_res if
                                  info_table[sn]['atom_name'] not in ['N', 'H']]

                if next_amide is not None:
                    return same_chain_atms[next_amide:] + amide_excluded
                else:
                    return amide_excluded
                    print('reached the n terminus -----\n\n\n')

        else:

            if atom_name == 'C':  # grab everthing after carbonyl including oxygen

                next_res = kwargs['res_num'] + 1
                # if carboxyl O is there, it will right after the C.
                # include it.
                if info_table[atom_sn+1]['atom_name'] == 'O':
                    oxy = [atom_sn+1]
                next_amide = [sn for sn in same_chain_atms if (info_table[sn]
                              ['res_num'] == next_res and
                              info_table[sn]['atom_name'] == 'N')]
                if next_amide != []:
                    try:
                        return same_chain_atms[next_amide[0]:] + oxy
                    except:
                        return same_chain_atms[next_amide[0]:]
                else:
                    raise ValueError

            else:  # grab every residue before + the amide hydrogen

                cur_n = same_chain_atms.index(min(cur_res))
                cur_h = [same_chain_atms.index(
                    sn) for sn in cur_res if info_table[sn]['atom_name'] == 'H']

                if cur_h != []:
                    return same_chain_atms[0:cur_n] + cur_h
                else:
                    return same_chain_atms[0:cur_n]
    else:
        # in the future I will implement a way to directly manipulate side
        # chain rotamers but not now
        print(f'Attempt was made at segment selection on an atom named {atom_name} ,'
              f' with anchor = {anchor}. Currently, selection with atom names'
              'and/or anchors not either "C" or "N" is not implemented.')
        return None

    # depending on the anchor these may or may not be included

def recursive_search(iterable, target):
    for item in iterable:
        if isinstance(item, (list, tuple)):
            # If the item is iterable, recursively search within it
            found = recursive_search(item, target)
            if found:
                return True
        elif item == target:
            # If the item matches the target, return True
            return True
    # If the target is not found in the current iterable or its sub-iterables, return False
    return False


def measure_dihedral(set_of_points, **kwargs):
    '''
    Given points A, B, C, D as input, each vectors in R3, and for which no
    pair of these points are the same, and for which no combination of
    three of these points are co-linear to each other:
        (this conditional check will be handled by the calling function)

    The vector B to C forms the axis of rotation.
    Rturn the angle between plane ABC and plane BCD relative to the axis of
    rotation. The dihedral angle is the arccosine of the dot products of the
    normal vectors of this plane. the sign of the cross product bweteen the
    normal vectors determines if the angle is positive or negative relative
    to the axis of rotation

    Parameters
    ----------
    set_of_points : list of exactly four R3 vectors which must be connected
    in the arrangement seen below
      B ____ C
     /********\
    A**********D

    **kwargs : list of dictionaries. - optional - Only the keyword 'debug_mode'
    is used to provide verbose output. Default is None.

    Returns

       float64 the dihedral angle in radians between -pi and +pi

    '''
    a_to_b = set_of_points[1] - set_of_points[0]
    a_to_c = set_of_points[2] - set_of_points[0]
    b_to_c = set_of_points[2] - set_of_points[1]
    b_to_d = set_of_points[3] - set_of_points[1]
    rotation_axis = -b_to_c
    c_to_d = set_of_points[3] - set_of_points[2]

    normal_plane1 = vector(np.cross(a_to_b, a_to_c))
    normal_plane2 = vector(np.cross(b_to_c, b_to_d))

    crossp = np.cross(normal_plane1, normal_plane2)
    theta = normal_plane1.angle_between(normal_plane2)

    # the way to determine whether this angle should be a positive or negative
    # is to see if the cross product of the plane normals is parallel or
    # anti-parallel to the axis of rotation as determined by the dot product
    if np.dot(np.cross(normal_plane1, normal_plane2), rotation_axis) > 0:
        theta = -theta

    return theta


def get_phi_psi(coordinates, info_table, **kwargs):

    # !!!!!! re-write this to only take serial numbers this is getting too
    # complicated to handle this input
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

    # -------------- mode 1, get values by serial numbers
    if 'ser_num' in kwargs:
        kwargs.pop('debug_mode', False)
        if len(kwargs.items()) > 1:
            raise ValueError('If serial numbers are specified as the selection'
                             ' criterion, no other criteria can be given')
        else:
            ser_nums = kwargs.pop('ser_num')
            if isinstance(ser_nums, int):
                ser_nums = [ser_nums]
            alpha_sns = []
            for sn in ser_nums:

                # if it is not an alpha carbon, find the CA of the residue it
                # is on
                k = info_table[sn]
                if k['atom_name'] != 'CA':
                    kwargs['submodel'] = k['submodel']
                    kwargs['model_name'] = k['model_name']
                    kwargs['res_num'] = k['res_num']
                    kwargs['chain'] = k['chain']
                    kwargs['atom_name'] = 'CA'
                    ca_sn = select(info_table, **kwargs)[0]
                else:
                    ca_sn = sn

                alpha_sns.append(ca_sn)
            bb_ser_num_sets = get_dihedral_backbones(info_table, alpha_sns)

            # check if phi and psi are definable for this residue, then
            # pass backbone atom ser_nums to measure_dihedral to get them

            for group in bb_ser_num_sets:

                if group[0] is not None:  # if there was a prev residue
                    vects = [coordinates[sn] for sn in group[0:4]]
                    # pass first 4
                    phi = measure_dihedral(vects, debug_mode=True)
                    for sn in group[1:4]:
                        info_table[sn]['phi'] = phi
                else:  # else phi is NaN
                    for sn in group[1:4]:
                        info_table[sn]['phi'] = float('nan')
                if group[-1] is not None:  # if there is a next residue
                    # pass last 4, reversed
                    vects = [coordinates[sn] for sn in group[1:]]
                    psi = measure_dihedral(vects)
                    for sn in group[1:-1]:
                        info_table[sn]['psi'] = psi
                else:  # else psi is NaN
                    for sn in group[1:-1]:
                        info_table[sn]['psi'] = float('nan')

    # ------------- mode 2, parse kwargs and find phi/psi for all matches
    else:
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
            pass
            print(f"the selected chains for get_phi_psi are {selected_chains}")

            # if atoms names were specified in kwargs ignore them and set them to
            # C alpha, otherwise specify c alpha
        kwargs['atom_name'] = 'CA'

        # if an explicit list of residues was passed, extract them
        if 'res_num' in kwargs:
            specified_residues = kwargs.pop('res_num', None)

            if isinstance(specified_residues, int):
                specified_residues = [specified_residues]

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
                          f"so the selected residues for chain {ch} are all "
                          f"of them:\n{residues_to_query}")

            # if residues were explicitly specified, iterate over those instead
            else:

                residues_to_query = specified_residues
                if isinstance(residues_to_query, int):
                    residues_to_query = [residues_to_query]
                if debug_mode:

                    print(
                        f'the residues being queried are:\n{residues_to_query}')

                    # for n in residues_to_query:

                    # get the ser nums for all C alphas on residues across all
                    # chains we specified
            kwargs['res_num'] = residues_to_query
            kwargs['atom_name'] = 'CA'
            ca_ser_nums = select(info_table, **kwargs)

            bb_ser_num_sets = get_dihedral_backbones(info_table, ca_ser_nums)

            for group in bb_ser_num_sets:

                if group[0] is not None:  # if there was a prev residue
                    vects = [coordinates[sn]
                             for sn in group[0:4]]  # pass first 4
                    phi = measure_dihedral(vects, debug_mode=True)
                    for sn in group[1:4]:

                        info_table[sn]['phi'] = phi
                else:  # else phi is NaN
                    for sn in group[1:4]:

                        info_table[sn]['phi'] = float('nan')
                if group[-1] is not None:  # if there is a next residue
                    # pass last 4, reversed
                    vects = [coordinates[sn] for sn in group[1:]]
                    psi = measure_dihedral(vects)
                    for sn in group[1:-1]:

                        info_table[sn]['psi'] = psi
                else:  # else psi is NaN
                    for sn in group[1:-1]:
                        info_table[sn]['psi'] = float('nan')


def set_phi_psi(coordinates, info_table, angle, angle_type='phi', anchor='N',
                **kwargs):
    """
    Rotate everything in the coordinates matrix so that the residue given has a
    phi or a psi angle as specified, with the angle specified in radians by the
    angle parameter. Which end will be held static while the other is rotated
    is specified by the anchor parameter
    """

    # NOTE: this method only works if the coordinates for the amide nitrogen
    # the amide carbon, and the alpha carbon are defined for every residue
    # if they aren't the complete model needs to be built. However checking
    # for this every time would be a major drain on computation time.

    # !!!!!!!!!!!!!!!  CURRENTLY UNUSED !!!!!!!!!!!!!!!!!!!!!!!
    def same_residue(info_table, ser_num):
        # given a serial number, return all serial numbers on the same
        # residue and the same chain
        res = info_table[ser_num]['res_num']
        chain = info_table[ser_num]['chain']
        mod = info_table[ser_num]['model_name']
        smod = info_table[ser_num]['submodel']

        # this will take advantage of the fact that the serial numbers
        # should always be in the order N, CA, C, O, then the rest
        same_res = select(info_table, res_num=res, chain=chain, model_name=mod,
                          submodel=smod)

        n = (info_table[sn]['ser_num'] for sn in same_res if
             info_table[sn]['atom_name'] == 'N')
        ca = (info_table[sn]['ser_num'] for sn in same_res if
              info_table[sn]['atom_name'] == 'CA')
        c = (info_table[sn]['ser_num'] for sn in same_res if
             info_table[sn]['atom_name'] == 'C')

        return None

    def find_rotation_angle(target_angle, ser_num):
        # given the residue pointed to by the ca_ser_num, and the desired
        # angle type, look up that angle, or call get_phi_psi if needed
        # to get it, and find the angle you must rotate by to end up
        # at the desired angle from the current one. THIS WILL ONLY work
        # properly if the kwargs point to only one residue on one chain

        # amide_nitrogen_sn is ser_nums [0]
        # c_alpha_sn  is ser_nums [1]
        # amide_carbon (carbonyl) is ser_nums[2]

        current_angle = info_table[ser_num][angle_type]
        if math.isnan(current_angle):
            kwargs.pop('debug_mode', None)

            get_phi_psi(coordinates, info_table, ser_num=ser_num)

            current_angle = info_table[ser_num][angle_type]
        # if it is still undefined, you cannot take the angle
        if math.isnan(current_angle):
            return None
        else:
            return (target_angle - current_angle)

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


    kwargs['atom_name'] = 'CA'
    # if an atom name was specified ignore it
    # and search for 'CA'
    ca_ser_nums = select(info_table, **kwargs)

    for ndx, ca in enumerate(ca_ser_nums):

        backbone_sn_groups = get_dihedral_backbones(info_table, [ca])

        for bb_sn in backbone_sn_groups:

            # bb_sn is a gorup of backbone serial numbers for a given residue
            # there are 5. in the following order [prev C, amide, alpha, C, next amide]
            # the first and last are none if there is no residue before or after

            # if phi was there was a prev residue and phi was
            # the selected angle

            if bb_sn[0] is not None and angle_type == 'phi':

                if info_table[bb_sn[2]]['res_name'] != 'PRO':
                    center = coordinates[bb_sn[2]].copy()
                    # after centering on alpha carbon this now defines the axis
                    # of rotation
                    amide = coordinates[bb_sn[1]] - center

                    # we select the atoms that need to have the transformation
                    # applied to them
                    rotation_segment = isolate_segment(
                        info_table, bb_sn[1], anchor)

                    rot_ang = find_rotation_angle(angle, bb_sn[2])
                    for ser_num in rotation_segment:
                        coordinates[ser_num] -= center

                        coordinates[ser_num] = coordinates[ser_num].rotate_arround(
                            rot_ang, -amide)

                        coordinates[ser_num] += center

            if bb_sn[4] is not None and angle_type == 'psi':
                # if we specified psi and there an amide at N +1 (ie not c term)
                center = coordinates[bb_sn[2]].copy()

                # after centering on alpha carbon this now defines the axis
                # of rotation
                carboxyl = coordinates[bb_sn[3]] - center

                # we select the atoms that need to have the transformation
                # applied to them, based on carboxyl
                rotation_segment = isolate_segment(
                    info_table, bb_sn[3], anchor)

                rot_ang = find_rotation_angle(angle, bb_sn[2])
                for ser_num in rotation_segment:

                    coordinates[ser_num] -= center

                    # for some reason, even though the axis was defined
                    # correctly, you need to reverse the sign on the rotation
                    # for psi
                    coordinates[ser_num] = coordinates[ser_num].rotate_arround(
                        -rot_ang, -carboxyl)
                    coordinates[ser_num] += center

            get_phi_psi(coordinates, info_table, ser_num=bb_sn[2])
            # --------------------- Omega ----------------

            if bb_sn[4] is not None and angle_type == 'omega':

                # if we specified omega and we are not at c term

                # !!!---------Note: this is not elegant -----------------
                # rewrite later to generalize this by extending bb_sn to 6
                # elements that includes the next CA instead of reselcting
                # again inside this code block

                # get all the info for the last element (the next amine)
                next_res = info_table[bb_sn[4]]['res_num']
                c = info_table[bb_sn[2]]['chain']
                smod = info_table[bb_sn[2]]['submodel']
                mod = info_table[bb_sn[2]]['model_name']
                # use it to find the next c_alpha

                next_ca = select(info_table, res_num=next_res,
                                 chain=c, submodel=smod, model_name=mod,
                                 atom_name='CA')
                vects = [coordinates[bb_sn[2]], coordinates[bb_sn[3]],
                         coordinates[bb_sn[4]], coordinates[next_ca[0]]]
                rot_ang = angle - measure_dihedral(vects)

                # for omega the carboyl is the center
                # we select the atoms that need to have the transformation
                # applied to them, based on carboxyl to N bond
                rotation_segment = isolate_segment(
                    info_table, bb_sn[4], anchor)

                # !!! also update find rotation angle to handle arbitrary
                # dihedrals instead of manually recalculating here
                carboxyl = coordinates[bb_sn[3]].copy()
                next_n = coordinates[bb_sn[4]].copy() - carboxyl

                for ser_num in rotation_segment:

                    coordinates[ser_num] -= carboxyl

                    # for some reason, even though the axis was defined
                    # correctly, you need to reverse the sign on the rotation
                    coordinates[ser_num] = coordinates[ser_num].rotate_arround(
                        rot_ang, -next_n)
                    coordinates[ser_num] += carboxyl
                # NOTE this assumes that there was no problem setting this
                # angle, and doesn't actually measure and check
                for s in bb_sn[1:4]:
                    info_table[s]['omega'] = measure_dihedral(vects)



def is_clashing(coordinates, search_set_indices, threshold):
    '''
    given 2 lists of serial numbers, return true if any pairwise combination
    of points from coordinates associated with a serial number in A is closer
    than threshold from any coordinate indicated by a serial number in list b

    Parameters
    ----------
    coordinates : vector
        all coordinates
    search_set_indices : list of lists of integers
        list containing exactly 2 lists each of which contains serial numbers.
    threshold : float64
        value below which to halt the process and report a clash.

    Returns
    -------
    bool
        if a clash was detected.

    if debug mode is true distances is where the measured distances
    will be recorded. Otherwise, only return true false


    '''

    # !!! in the future try to convert this to a heuristic kdtree algorithm
    # which is much faster but difficult to understand and implement.
    # right now it is brute force and ~100x slower than it could be.

    # unpack the input lists
    set_a, set_b = search_set_indices

    pairs = [[a, b] for a in set_a for b in set_b]
    for pair in pairs:

        if (coordinates[pair[1]] -
                coordinates[pair[0]]).get_length() < threshold:
            return True
    return False


# !!! in the future, I would like to change this to just take an atom
# serial number to allow it to be more generalizable, to check all serial
# number on one side of the atom against all others on the other side,but
# for now it needs to fit into existing code.

def check_internal_clash(coordinates, info_table, cutoff_distance, angle_type,
                         anchor, **kwargs):

    if angle_type == 'phi':
        # get the index of the nitro of the current res
        kwargs['atom_name'] = 'N'
        atom_sn = select(info_table, **kwargs)
        if len(atom_sn) > 1:
            raise ValueError('no more than 1 atom can be specified for segment'
                             'isolation')
    elif angle_type == 'psi':
        # get the index of the carbonyl of the current res
        kwargs['atom_name'] = 'C'
        atom_sn = select(info_table, **kwargs)

    # get the sets of serial numbers on the same chain on each side of
    # a given atom. Note that the atom must be a carbonyl or amide nitrogen
    set_a = isolate_segment(info_table, atom_sn[0])

    set_b = isolate_segment(info_table, atom_sn[0], anchor='C')

    if is_clashing(coordinates, [set_a, set_b], cutoff_distance):
        return True
    return False


def rot_sym_clash_check(coordinates, **kwargs):
    '''
    Checks if a given set of coordinates can produce radial symmetry arround
    the center vector, when multiplicity copies are arranged along the rot_ax.

    Parameters
    ----------
    coordinates : Array of vector
        a SUBSET of the total set of coordinates, which is the group
        of molecules to be considered for symmetry testing. coordinates should
        be in the standard cartesian basis.
    multiplicity : Int
        An integer 2 or greater that indicates the number of symmetrical units.
    center_v : vector
        center of symmetry for the multimer to be tested. This should
        generally be outside the bounding volume of the molecule, and in
        standard cartesian basis.
    radial_v :vector
        radial axis. should be orthogonal to rot_ax, and in
        standard cartesian basis.
    rot_ax : vector
        The axis arround which rotational symmetry of the moleule will be
        tested.
    cutoff_distance : float
        the distance threshold below which is considered a clash.

    Returns
    -------
    clash : boolean
        True if arranging n copies symmetrically arround the center of rotation
        would result in a clash, otherwise false

    '''
    multiplicity = kwargs.get('multiplicity', None)
    radial_v= kwargs.get('radial_v', None)
    rot_ax = kwargs.get('axis_vector', None)
    radius = kwargs.get('radius', None)
    n_slices = kwargs.get('n_slices', 20)
    n_rings = kwargs.get('n_rings', 20)
    cutoff_distance = kwargs.get('cutoff_distance', 0.36)
    center_v = kwargs.get('center_v', None)

    def onion_method():
        # print('--------------------------------------')
        # get the highest and lowest values of theta
        # an easy hueristic is that if the angle needed to contain
        # all the point of a monomer in an arc arround the center
        # is less than 2 pi / n, then that monomer can have n fold
        # radial symmetry without clashes between subunits
        a = np.argmin(cylind_coords[:, 1])
        s = cylind_coords[a, 1]
        b = np.argmax(cylind_coords[:, 1])
        v = cylind_coords[b, 1]
        threshold = 2*np.pi/multiplicity
        if abs(v - s) < threshold:

            return False

        else:
            # sort cylindrical coordinates by z and bin them
            # find the range for each bin
            ndx_by_height = np.argsort(cylind_coords[:, 0])
            min_z = cylind_coords[ndx_by_height[0]][0]
            max_z = cylind_coords[ndx_by_height[-1]][0]

            dz = (max_z - min_z)/n_slices

            height_ndx_bins = []
            numel = len(ndx_by_height)
            k = 0
            for i in range(n_slices):
                height_ndx_bins.append([])
                low_limit = min_z + i*dz
                high_limit = low_limit + dz
                while k < numel and cylind_coords[ndx_by_height[k]][0] \
                    >= low_limit and cylind_coords[ndx_by_height[k]][0] \
                        < high_limit:

                    height_ndx_bins[i].append(ndx_by_height[k])
                    k += 1
            # because it is strictly less than the high limit,
            # the last ordered index must be tacked on to the last bin

            # if height_ndx_bins != []:
            height_ndx_bins[-1].append(ndx_by_height[-1])

            # now go divide your onion slice into onion "rings"
            angle_diffs = []
            for i in range(n_rings):

                # the indices in current slice
                curr_slice = height_ndx_bins[i].copy()
                # if one or fewer points in a slice, you cant
                # find an angle, but we know it will not
                # cause a clash, so skip it.
                if len(curr_slice) > 1:
                    # the actual coordinates of those
                    slice_coords = [cylind_coords[q] for
                                    q in curr_slice]
                    slice_coords = np.asarray(slice_coords)

                    # the indices of the indices in the slice sorted by radius
                    ndx_by_rad = np.argsort(slice_coords[:, 2])

                    numel = len(ndx_by_rad)

                    min_r = slice_coords[ndx_by_rad[0]][2]
                    max_r = slice_coords[ndx_by_rad[-1]][2]

                    dr = (max_r - min_r)/n_rings

                    rad_ndx_bins = []
                    k = 0
                    for j in range(n_rings):
                        rad_ndx_bins.append([])
                        low_limit = min_r + j*dr
                        high_limit = low_limit + dr

                        while k < numel and slice_coords[ndx_by_rad[k]][2] \
                            >= low_limit and slice_coords[ndx_by_rad[k]][2] \
                                < high_limit:
                            rad_ndx_bins[j].append(ndx_by_rad[k])
                            k += 1
                    rad_ndx_bins[-1].append(ndx_by_rad[-1])

                    # now go through each onion ring and if there is more than 1
                    # index, check if their angles indicate they would clash with
                    # another multimer

                    for j in range(n_rings):
                        # if there are 2 or more get all pairs of indices
                        if len(rad_ndx_bins[j]) > 1:
                            index_pairs = itertools.combinations(
                                rad_ndx_bins[j], 2)
                            # if the difference in angle between any 2 points
                            # on the same ring is over 2pi/multiplicity,
                            # there is a clash
                            for ndx_a, ndx_b in index_pairs:
                                ang_a = slice_coords[ndx_a][1]
                                ang_b = slice_coords[ndx_b][1]
                                diff = convert_angle(ang_a - ang_b)

                                angle_diffs.append(diff)

                                if diff <= -np.pi/multiplicity/2 or diff \
                                        >= np.pi/multiplicity/2:

                                    return True

                    height_ndx_bins[i] = rad_ndx_bins

        # if the function got this far there is no clash
        return False

    # copy subset of coordinates and recenter on center_v
    cylind_coords = coordinates.copy()
    cylind_coords -= center_v

    if not rot_ax.is_nonzero():
        raise ValueError('The axis of rotation cannot be the zero vector')
    if not radial_v.is_nonzero():
        raise ValueError('The radial axis cannot be the zero vector')

    # convert to cylindrical coordinates, and change basis so that rot_ax is
    # directly along the cylindrical axis, and all angles are relative to
    # radial_v
    cylind_coords = cart2cylind(cylind_coords, center_v, rot_ax, radial_v)

    will_clash = onion_method()

    return will_clash


def random_backbone(coordinates, info_table, n_structs, residue_list_filename,
                    **kwargs):

    f = open(residue_list_filename, 'r')
    residue_list = f.readlines()
    axis_vector = kwargs.get('axis_vector', None)
    radial_v = kwargs.get('radial_v', None)
    cutoff_distance = kwargs.get('cutoff_distance', 0.36)
    chainlist = kwargs.get('chainlist', None)
    max_tries = kwargs.get('max_tries', 20)
    center_sn = kwargs.get('center_sn', None)
    radial_sn = kwargs.get('radial_sn', None)

    symmetry_groups = kwargs.get('symmetry_groups', 'all')
    radius = kwargs.get('radius', None)
    multiplicity = kwargs.get('multiplicity', None)

    if not (center_sn and radial_sn):
        raise ValueError('rotational symmetry clash check requires the '
                         'additional parameters center_sn and radial_sn as kwargs')

    # !!! TO DO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # change it to accept input file format that include the chain along with
    # the residue designations for fixed and anchored residues, so that it can
    # be determined which ones are symmetrically related to which. This will
    # allow heterogenous set of chains while still allowing for assumptions
    # about symmetry to be used to simplify clash checking.
    # !!!

    # symmetry_groups will be a list of lists each containing Chain letters of
    # chains which are related to each other by symmetry this will be used to
    # expedite clash checking. if none are given, all chains will be treated
    # as identical in sequence qand related by symmetry together in one
    # symmetry group, else if the value None is passed, or an empty list is
    # passed, they will all be treated as different

    # residue numbers after the header 'anchored:' in the residue file will
    # be put in a separate list, these residues cannot move at all
    anchor_header_line = None

    # residue numbers after the header 'fixed:' in the residue file will be in
    # a different group. these can move as a result of adjusting other
    # residues, but cannot have their own internal coordinates adjusted
    fixed_header_line = None

    # trim off the newline and convert to integer for the imported residue
    # numbers whose backbones you wish to randomize

    # if your randomized backbone needs to consider self-symmetry in its clash
    # evaluation additional parameters are needed for the hueristic clash check
    # which should be passed as a list.

    # print(f'the kwargs are:\n\n\n\{kwargs}\n\n\n\n\n')
    if center_sn and radial_sn:

        kwargs['center_v'] = coordinates[center_sn]

        def symmetry_clash(coordinates):
            return rot_sym_clash_check(coordinates, **kwargs)
    else:
        # if symmetry is not to be considered just return false

        def symmetry_clash(coordinates):
            return False

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

    free_residues = []
    for sn_set in sn_by_chain:
        free_residues.append({info_table[s]['res_num'] for s in sn_set})

    # whittle away the anchored residues from the list to iterate
    # over when randomizing angles
    free_residues = [[s for s in sn_set if s not in anchored_residues and s
                      not in fixed_residues] for sn_set in free_residues]

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

        c_anchor = max(anchored_residues)

    # if an anchored segment exists, the free residues are split into those on
    # the n_terminal side of it (c_anchored) or the c terminal side
    # (n_anchored)

    c_anchored_segments = []
    n_anchored_segments = []


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

                seg.sort(reverse=True)
                c_anchored_segments.append(seg)
            if c_anchor is not None:
                # repeat this for n anchored segments
                seg = [res for res in free_residues[i] if res > c_anchor]
                seg.sort(reverse=False)
                n_anchored_segments.append(seg)

    if len(c_anchored_segments) != len(n_anchored_segments) and \
            len(free_residues) != len(chainlist):
        raise ValueError('the number of lists in c_anchored_segments, '
                         'n_anchored_segments and the number of items in '
                         'chainlist must all be equal')

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
                        print(f'\n---model: {ndx} residue: {num} ', end='')
                        tries = 0
                        while ( #symmetry_clash(coordinates)
                               #or
                               check_internal_clash(coordinates,
                                                       info_table,
                                                       cutoff_distance,
                                                       angle_type='phi',
                                                       anchor='C',
                                                       res_num=num,
                                                       chain=chainlist[i])) \
                                and tries < max_tries:

                            set_phi_psi(coordinates, info_table,
                                        random.uniform(0, 2*np.pi),
                                        angle_type='phi', anchor='C',
                                        res_num=num, chain=chainlist[i])

                            tries += 1
                            if tries < max_tries:
                                print('.', end='')
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
                        print(f'\n---model: {ndx} residue: {num} ', end='')
                        tries = 0
                        set_phi_psi(coordinates, info_table,
                                    random.uniform(0, 2*np.pi),
                                    angle_type='psi', anchor='C',
                                    res_num=num, chain=chainlist[i])
                        while (
                                # symmetry_clash(coordinates) or
                                check_internal_clash(coordinates,
                                                       info_table,
                                                       cutoff_distance,
                                                       angle_type='psi',
                                                       anchor='C',
                                                       res_num=num,
                                                       chain=chainlist[i])) \
                                and tries < max_tries:
                            set_phi_psi(coordinates, info_table,
                                        random.uniform(0, 2*np.pi),
                                        angle_type='psi', anchor='C',
                                        res_num=num, chain=chainlist[i])
                            tries += 1
                            if tries < max_tries:
                                print('.', end='')
                            else:
                                print(f'The maximum number of attempts to '
                                      f'unclash model number {ndx} was '
                                      'reached. The model has been discarded.')
                    except ValueError:
                        print(f'Cant set psi for C terminus at residue {num} '
                              f'of model number {ndx}. Skipping this residue.')
                        pass
            if n_anchored_segments[i] != []:
                for num in n_anchored_segments[i]:
                    try:
                        set_phi_psi(coordinates, info_table,
                                    random.uniform(0, 2*np.pi),
                                    angle_type='phi', anchor='N',
                                    res_num=num, chain=chainlist[i])
                        print(f'\n---model: {ndx} residue: {num} ', end='')
                        # if your backbone adjustment resulted in a clash, retry
                        # with a different random angle

                        tries = 0
                        while (#symmetry_clash(coordinates)
                              # or
                               check_internal_clash(coordinates,
                                                       info_table,
                                                       cutoff_distance,
                                                       angle_type='phi',
                                                       anchor='N',
                                                       res_num=num,
                                                       chain=chainlist[i])) \
                                and tries < max_tries:

                            set_phi_psi(coordinates, info_table,
                                        random.uniform(0, 2*np.pi),
                                        angle_type='phi', anchor='N',
                                        res_num=num, chain=chainlist[i])

                            tries += 1
                            if tries < max_tries:
                                print('.', end='')
                            else:
                                print(f'\nThe maximum number of attempts to '
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
                        while (#symmetry_clash(coordinates)
                               #or
                               check_internal_clash(coordinates,
                                                       info_table,
                                                       cutoff_distance,
                                                       angle_type='psi',
                                                       anchor='N',
                                                       res_num=num,
                                                       chain=chainlist[i])) \
                                and tries < max_tries:
                            set_phi_psi(coordinates, info_table,
                                        random.uniform(0, 2*np.pi),
                                        angle_type='psi', anchor='N',
                                        res_num=num, chain=chainlist[i])
                            tries += 1
                            if tries < max_tries:
                                print('.', end='')
                            else:
                                print(f'\nThe maximum number of attempts to '
                                      f'unclash model number {ndx} was '
                                      'reached. The model has been discarded.')
                    except ValueError:
                        print(f'Cant set phi for C terminus at residue {num} '
                              f'of model number {ndx}. Skipping this residue')
                        pass
            if free_residues != []:
                pass
    return coordinates, info_table

def get_transform_mat(source_coords, target_coords, source_center_sn, 
                      source_axis_sn, source_radial_sn, **kwargs):

    """
    Compute the 4x4 homogeneous transformation matrix that aligns the source coordinates to the target coordinates.

    The transformation matrix will rotate and translate the source points to match the target points.
    By default, both rotation and translation are applied. If only rotation is required, set the `translate` keyword argument to `False`.

    Parameters:
    ----------
    source_coords : numpy.ndarray
        The coordinates of the source object as a 2D array (n x 3).

    center_sn : int
        The index of the center point in the source coordinates.

    axis_sn : int
        The index of the axis point in the source coordinates.

    radial_sn : int
        The index of the radial point in the source coordinates.

    target_coords : numpy.ndarray
        The coordinates of the target object as a 2D array (n x 3).

    **kwargs : dict, optional
        A dictionary of keyword arguments:
        - `translate` (bool): If `True`, both rotation and translation are applied (default behavior).
        - If `False`, only rotation is applied, and the translation is ignored.

        - `target_center_sn` (int): If given, the index for the atom to be used as
            center coordinate, otherwise assigned the same as source_center_sn

        - `target_radial_sn' (int): If given, the index for the atom to be used as
            the radial axis coordinate, otherwise assigned the same as source_radial_sn

        - `target_axis_sn' (int): If given, the index for the atom to be used as
            the principle axis coordinate, otherwise assigned the same as source_axis_sn

            NOTE: either all three target coordinates keywords must be supplied, or none.
            incomplete set of coordinates raises an error

    Returns:
    -------
    transformation_matrix : numpy.ndarray
        A 4x4 homogeneous transformation matrix that aligns the source coordinates to the target coordinates.
        The matrix includes both rotation and, if requested, translation.

    Raises:
    ------
    ValueError
        If the number of arguments in `args` is not exactly 3 or if the input data is invalid.
    """

    translate = kwargs.get('translate', True)
    target_radial_sn = kwargs.get('target_radial_sn', None)
    target_center_sn = kwargs.get('target_center_sn', None)
    target_axis_sn = kwargs.get('target_axis_sn', None)

    v = [tf is None for tf in (target_axis_sn, target_center_sn, target_radial_sn)]
    if all(v):
        target_radial_sn = source_radial_sn
        target_center_sn = source_center_sn
        target_axis_sn = source_axis_sn
    elif any (v) and not all(v):
        raise ValueError ('missing one or more indices to define the target axis, if none are given the same set from the source structure is used by default')

    if translate == False:
        translation_vect = np.asarray([0, 0, 0]).T
        # if translation to superimpose the two structures is not desired
        # change the last column to zeros except the last element
        # [[ .  .  .  0 ]
        #  [ .  .  .  0 ]
        #  [ .  .  .  0 ]
        #  [ 0  0  0  1 ]]
    else:

        # translation vect is the difference vector between their center_atoms
        # as a column vector
        translation_vect = source_coords[source_center_sn].copy().reshape(-1,1) \
            - target_coords[target_center_sn].reshape(-1,1)
    print(translation_vect)
    # define internals for the source structure ----------------
    print(source_axis_sn, source_center_sn, source_radial_sn, target_axis_sn, target_center_sn, target_radial_sn)
    # the priciple axis that defines "up"
    print(source_coords[source_axis_sn], source_coords[source_center_sn], 
          source_coords[source_radial_sn])
    print(target_coords[target_axis_sn], target_coords[target_center_sn], 
          target_coords[target_radial_sn])
    source_principle_axis = vector(source_coords[source_axis_sn] -
                                   source_coords[source_center_sn])
    source_principle_axis = source_principle_axis.unitize()
    source_principle_axis = np.squeeze(source_principle_axis)
    print(f'source principle {source_principle_axis}')

    source_radial_axis = vector(source_coords[source_radial_sn] -
                                source_coords[source_center_sn])
    source_radial_axis = np.squeeze(source_radial_axis)
    source_radial_axis = source_radial_axis.project_onto_normal_plane(source_principle_axis)
    source_radial_axis = source_radial_axis.unitize()

    source_norm_axis = vector(np.cross(source_principle_axis,
                                       source_radial_axis))

    source_mat = np.vstack((source_principle_axis, source_radial_axis,
                            source_norm_axis))

    print('source_mat\n', source_mat, source_mat.shape)


    # define internals for the target structure ----------------

    target_principle_axis = vector(target_coords[target_axis_sn] -
                                   target_coords[target_center_sn])
    target_principle_axis = target_principle_axis.unitize()
    target_principle_axis = np.squeeze(target_principle_axis)
    print(f'target principle {target_principle_axis}')

    target_radial_axis = vector(target_coords[target_radial_sn] -
                                target_coords[target_center_sn])
    target_radial_axis = np.squeeze(target_radial_axis)
    target_radial_axis = target_radial_axis.project_onto_normal_plane(target_principle_axis)
    target_radial_axis = target_radial_axis.unitize()

    target_norm_axis = vector(np.cross(target_principle_axis,
                                       target_radial_axis))
    target_mat = np.vstack((target_principle_axis, target_radial_axis,
                            target_norm_axis))

    print('target_mat\n', target_mat, target_mat.shape)

    rot_mat = target_mat @ source_mat.T
    print('rotation matrix\n', rot_mat, rot_mat.shape)
    aug_rot_mat = np.hstack((rot_mat, translation_vect))
    # print('\naugmented rotation matrix\n', aug_rot_mat, aug_rot_mat.shape)

    B = np.asarray([0, 0, 0, 1])
    # print('\nB\n', B, B.shape)
    # print(transformation_matrix, transformation_matrix.shape)
    transformation_matrix = np.vstack((aug_rot_mat, B))
    # print(f'transformation matrix: {transformation_matrix}')
    return transformation_matrix

def orient(coordinates, target_coords, center_sn, axis_sn, radial_sn,
           **kwargs):

    # pass along kwargs received directly. to generate transformation matrix according
    # to kwargs passed to this function.
    transformation_matrix = get_transform_mat(coordinates, target_coords, 
                                              center_sn, axis_sn, radial_sn, 
                                              **kwargs)

    # pad on extra 1s for homogeneous coordinates
    padding = np.ones((len(coordinates), 1), dtype=float)
    padded_coordinates = np.hstack((coordinates, padding))

    padded_coordinates = np.vstack((padded_coordinates, vector([0, 0, 0, 1])))
    # print(f'size of coordinates {padded_coordinates.shape}')

    transformed_coordinates = padded_coordinates@transformation_matrix
    transformed_coordinates = transformed_coordinates[:-1, :3]

    return transformed_coordinates


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
        atoms of the helix that indicates the direction the helix is pointing
        in
    '''
    # get ser nums of all alpha carbons in that range
    c_alpha_ser_nums = select(info_table, res_num=(list(range(nterm_res,
                                                              cterm_res+1))),
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
    elif type(chain_list) is str:
        chain_list = [chain_list]
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    availible_letters = [
        letter for letter in alphabet if letter not in existing_chains]
    new_chain_list = []
    for chain in chain_list:

        # assign a new letter for the new chain and remove that letter from the
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
        # make temp array that is slice of coordinates array where chain
        # is stored then concetenate the two
        new_coords = coordinates[mother_chain_sn[0]: mother_chain_sn[-1]+1]
        coordinates = np.concatenate((coordinates, new_coords)).view(vector)
        info_table.extend(copy.deepcopy(
            info_table[mother_chain_sn[0]: mother_chain_sn[-1]+1]))
        for i, sn in enumerate(daughter_chain_sn):
            info_table[sn]['ser_num'] = sn
            info_table[sn]['chain'] = new_chain
    return coordinates, info_table, new_chain_list


def axial_symmetry(coordinates, info_table, multiplicity, radius,
                   center_sn, axis_sn, radial_sn, **kwargs):
    '''
        Parameters
        coordinates : array of vectors
            the molecular coordinates
        info_table : list of dictionaries
            list containing info on each coordinate
        multiplicity : INT16
            how many fold symmetry we want
        radius : float64
            the length of the radius from the center atom of our chain to the
            center of rotation
        center_sn : INT16
            serial number specifying the atom which will be used as the basis
            of the coordinates for the chain which we will recenter on
        axis_sn : INT16
            serial number specifying the atom which will be used to form the
            longitudinal axis by the vector axis atom - center atom
        radial_sn : INT16
            serial number specifying the atom which will be used to form the
            radial axis by the vector radial atom - center atom
        kwargs ~~~~~~~~
        chains : LIST
            list containing strings indicating chains that will be
            multimerized.
        multiplicity : INT16
            positive integer describing how many subunits in the arrangement
        cofr : VECTOR, optional
            the center of rotation which the transformation will be applied
            arround. The default is vector([0.0, 0.0, 0.0]).
        rot_axis : VECTOR, optional
            a vector that serves as axis of rotation after recentering on cofr.
            The default is vector([0.0, 0.0, 1.0]).
        check_method : string, optional
            either 'heuristic' which indicates that the onion method will be
            used for computational speed, or 'definitive' which indicates that
            exhaustive pairwise clash checking will be performed. 'definitive'
            is the default.
        threshold : float, optional
            the distance in angstroms below which a clash should be reported.
            propogated to the clash check functions. default is 0.36
        radial_angle : float
            the angle in radians between the vector from canter atom to
            rotation center and the radial axis

        Returns
        -------
        Coordinates, info_table.

        '''

    # parse kwargs
    n_slices = kwargs.get('n_slices', 20)
    n_rings = kwargs.get('n_rings', 20)
    cofr = kwargs.get('cofr', vector([0, 0, 0]))
    radial_angle = kwargs.get('radial_angle', 0)
    rot_axis = kwargs.get('rot_axis', vector([0, 0, 1]))
    check_method = kwargs.get('check_method', 'definitive')
    threshold = kwargs.get('threshold', 0.36)

    try:
        rot_axis = rot_axis.unitize()
    except ZeroDivisionError:
        print('The axis of rotation can not be the zero vector')

    if multiplicity < 2 or type(multiplicity) is not int:
        print('multiplicity for axial symmetry must be an integer 2 or'
              ' greater')
        raise ValueError

    # first orient the chain specified by your reference vector points

    coordinates, radial_axis = orient(coordinates, info_table, center_sn,
                                      axis_sn, radial_sn, radial_angle, rot_axis, recenter=True)

    # the center of rotation is an eigenvector of the radial axis
    symmetry_center = radial_axis*radius
    # duplicate the chain, and keep a list of chains related by symmetry

    # call rot_sym_clash_check to ensure that the structure that we have for
    # our monomer would not result in a clash if n fold symmetry was applied
    # rot_sym_clash_check(coordinates, info_table, multiplicity,
    #                     radial_axis, rot_axis, 10.3, n_slices, n_rings)

    chain_group = []
    chain = info_table[center_sn]['chain']
    chain_group.append(chain[0])

    # get ser nums of everything on your chain of interest, then isolate it
    # as another array
    isolated_chain_ser_nums = select(info_table, chain=chain)

    isolated_chain = coordinates[isolated_chain_ser_nums]
    # if they pass the rotational symmetry check you are good to apply copying
    # and rotation, otherwise return None

    def definitive_test():
        angle = 2*np.pi/multiplicity
        # this copies the entire set of coordinates and info_table with the
        # new chain added, but is simpler given how clone chain was implemented

        tmp = coordinates.copy()
        tmp_table = info_table.copy()
        clone_coords, clone_info, new_chain = clone_chain(tmp, tmp_table)

        # make a single copy, and apply 1 step of n fold symmetry, so that
        # the new subunit is neighboring the old one. If these 2 subunits don't
        # clash, the nature of rotational symmetry ensures that no other ones
        # do either
        clone_sns = select(clone_info, chain=new_chain)

        for sn in clone_sns:
            clone_coords[sn] -= symmetry_center

        # check every point in the start chain against every point in its
        return is_clashing(clone_coords, [isolated_chain_ser_nums, clone_sns],
                           threshold)

    if check_method == 'heuristic':
        if rot_sym_clash_check(coordinates, multiplicity, radial_axis,
                               rot_axis, radius, center_sn, int(
                                   n_slices), int(n_rings),
                               center_v=symmetry_center):
            print('clash in symmetry detected')
            return None
        else:
            try:
                for x in range(0, (multiplicity - 1)):
                    coordinates, info_table, new_chain = clone_chain(
                        coordinates, info_table, chain)
                    chain_group.append(new_chain[0])
            except ValueError:
                print(f'chain {chain} does not exist or can not be copied.'
                      ' skipping.')
            for n, chain_id in enumerate(chain_group):
                angle = 2*np.pi/multiplicity*n
                current_chain = select(info_table, chain=chain_id)
                for sn in current_chain:
                    coordinates[sn] -= symmetry_center
                    coordinates[sn] = coordinates[sn].rotate_arround(
                        angle, rot_axis)

    elif check_method == 'definitive':
        if definitive_test():
            print('clash in symmetry detected')

            ###### for testing only ####
            # try:
            #     for x in range(0, (multiplicity - 1)):
            #         coordinates, info_table, new_chain = clone_chain(
            #             coordinates, info_table, chain)
            #         chain_group.append(new_chain[0])
            # except ValueError:
            #     print(f'chain {chain} does not exist or can not be copied.'
            #           ' skipping.')
            # for n, chain_id in enumerate(chain_group):
            #     # print(f'n : {n}')
            #     angle = 2*np.pi/multiplicity*n
            #     current_chain = select(info_table, chain=chain_id)
            #     for sn in current_chain:
            #         coordinates[sn] -= symmetry_center
            #         coordinates[sn] = coordinates[sn].rotate_arround(
            #             angle, rot_axis)
            # plot_model(coordinates)
            #####################################
            return None
        else:
            try:
                for x in range(0, (multiplicity - 1)):
                    coordinates, info_table, new_chain = clone_chain(
                        coordinates, info_table, chain)
                    chain_group.append(new_chain[0])
            except ValueError:
                print(f'chain {chain} does not exist or can not be copied.'
                      ' skipping.')
            for n, chain_id in enumerate(chain_group):
                angle = 2*np.pi/multiplicity*n
                current_chain = select(info_table, chain=chain_id)
                for sn in current_chain:
                    coordinates[sn] -= symmetry_center
                    coordinates[sn] = coordinates[sn].rotate_arround(
                        angle, rot_axis)
    else:
        raise ValueError("options for check_method argument are"
                         "'definitive' or 'heuristic'")
    return coordinates, info_table

def ref_point_input_check (coordinates, info_table, *ser_nums, **kwargs):

    check_chain = kwargs.get('check_chain', True)
    check_model = kwargs.get('check_model', True)
    check_submodel = kwargs.get('check_submodel', True)

    # get chain for the inputs, raise a warning if points aren't unique and
    # on the same chain
    if check_chain:
        chains = [info_table[sn]['chain'] for sn in ser_nums]
        if not all( ch == chains[0] for ch in chains):
            raise ValueError('The supplied reference points for basis coordinates'
                             'are not from the same chain. If this is okay, '
                             'pass the check_chain=False kwarg')
    if check_model:
        models = [info_table[sn]['model'] for sn in ser_nums]
        if not all( m == models[0] for m in models):
            raise ValueError('The supplied reference points for basis coordinates'
                             'are not from the same model. If this is okay, '
                             'pass the check_model=False kwarg')
    if check_submodel:
        subs = [info_table[sn]['submodel'] for sn in ser_nums]
        if not all( s == subs[0] for s in subs):
            raise ValueError('The supplied reference points for basis coordinates'
                             'are not from the same submodel. If this is okay, '
                             'pass the check_submodel=False kwarg')

    return True

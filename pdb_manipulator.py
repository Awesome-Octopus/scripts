#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 01:36:49 2023

@author: andrew
"""

# %% Header
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
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt


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
        if self.get_length() < 1e-7:
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
                # print(
                #     f'the vector is now {self.get_length()} long\n the other vector is {other_vector.get_length()} long')
                return np.arccos(c)
        except ZeroDivisionError as err:
            print(
                'Error in vector.angle_between(): Attempt to take an angle '
                'with the zero vector\n',
                err)

    def is_orthogonal(self, other_vector):
        v = np.dot(self, other_vector)
        f_point_tollerance = 1e-7
        if v < f_point_tollerance and \
           v > -f_point_tollerance:
            return True
        else:
            return False

    def rotate_arround(self, angle, axis_vector=[]):

        # rotate a 2d or a 3d vector arround another 2d or 3d vector by the
        # given angle according to right hand rule. intentionally
        # not defined for vectors of rank 4 or greater
        # print(f'the starting vector is {self.get_length()} long')
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
                    a = np.matmul(rot_mat, self)
                    # print(f'afterwards it is now {a.get_length()} long')
                    return a
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
            return np.matmul(rot_mat, self)
        else:
            raise ValueError(
                'Error in vector.rotate_arround(): No axis of rotation '
                'provided for a rank 3 vector\n')

    def project_onto_normal_plane(self, other_vector):

        # given some vector and a second vector return a vector which is the
        # component of the first vector which is orthogonal to the second, in
        # other words, if the second vector defined the normal of a plane, give
        # the shadow of the first vector projected onto that plane
        # print(f'before projection, length is {self.get_length()}')
        a = self.dot(other_vector.unitize())
        shadow = self - a*other_vector
        # print(
        #     f'after projection its shadow is {shadow.get_length()} long, while it is hopefully still only {self.get_length()} long')
        return shadow
# %% Coordinate Functions


def change_basis(coordinates, center_v, x_prime, y_prime,
                 z_prime):
    '''
    given a serial number for a central atom and a serial number for an atom to
    define a long axis, orient everything on that chain so that the vector from
    the center atom to the axis atom alligns with the reference vector.
    Parameters
    ----------
    coordinates : Array of float 64
        the atomic coordinates of all atoms
    info_table : list
        list of dictionaries containing data for each coordinate
    center_sn : int
        serial number of the atom which is to be centered at [0 0 0]
    allignment_sn : int
        serial number of the atom which will define the new + z axis relative
        to the center atom
    reference_vect : vector
        the axis relative to the standard basis set that you wish to align
        the long axis to. If none given, [0 , 0, 1].


    Returns
    -------
    coordinates : Array of float 64
        The new coordinates after allignment

    '''

    # input testing
    ref_points = [center_v, x_prime, y_prime, z_prime]
    new_basis = [x_prime, y_prime, z_prime]
    if not all(a.is_orthogonal(b) for a in new_basis for b in
               new_basis if a is not b):
        raise ValueError('new basis vectors must be orthogonal')

    if any(not a.is_nonzero() for a in [x_prime, y_prime, z_prime]):
        raise ValueError('basis vectors can not be zero')

    x_ax = vector([1, 0, 0])
    y_ax = vector([0, 1, 0])
    z_ax = vector([0, 0, 1])

    x_scale = x_prime.get_length()
    y_scale = y_prime.get_length()
    z_scale = z_prime.get_length()

    copy_coords = coordinates.copy()

    # print(f'change basis center_v: {center_v}')
    # for n in range(len(copy_coords)):
    #     if np.array_equal(center_v, copy_coords[n]):
    #         c = copy_coords[n].copy()
    copy_coords -= center_v

    # if bases within fp error of one another, just return
    # the input
    if np.allclose(new_basis, [x_ax, y_ax, z_ax]):
        return copy_coords
    else:
        x_prime = x_prime.unitize()
        y_prime = y_prime.unitize()
        z_prime = z_prime.unitize()

    # if the dot products of the unitized vectors and their coresponding basis
    # vectors are all either 1 or -1 , we dont need to perform rotation, only
    # scalar multiplication along an axis
    d_prods = [np.dot(a, b) for a, b in zip([x_prime, y_prime, z_prime],
                                            [x_ax, y_ax, z_ax])]
    # print(f'dot products {d_prods}')
    if all(np.allclose(a, 1) or np.allclose(a, -1) for a in d_prods):

        # print('reflection')
        # perform a reflection
        pass

    # if the z axis isn't already alligned, do that first
    elif not np.allclose(z_prime, z_ax):

        xy = z_prime.project_onto_normal_plane(z_ax)
        # print(f'xy projection {xy}')

        if not xy.is_nonzero():
            yz_ang = None

        # if the bases aren't antiparalell or paralell
        elif not np.allclose(-z_prime, z_ax):
            yz_ang = y_ax.angle_between(xy)
            crossp = vector(np.cross(xy, y_ax))
            # if the Z component of the cross product is negative, reverse
            # the direction
            if crossp[2] < 0:
                yz_ang = -yz_ang

            # print(f'yz_ang {yz_ang}')
            # rotate everything around the z axis by that angle
            for n in range(len(copy_coords)):
                if copy_coords[n].is_nonzero():
                    copy_coords[n] = copy_coords[n].rotate_arround(
                        yz_ang, z_ax)
                x_prime = x_prime.rotate_arround(yz_ang, z_ax)
                y_prime = y_prime.rotate_arround(yz_ang, z_ax)
                z_prime = z_prime.rotate_arround(yz_ang, z_ax)
                # print(
                #     f'x prime: {x_prime}, y prime: {y_prime}, z prime: {z_prime}')
            z_ang = z_prime.angle_between(z_ax)
            crossp = vector(np.cross(z_prime, z_ax))

            # if cross product is on +x rotation is positive
            if crossp[0] < 0:
                z_ang = -z_ang
            for n in range(len(copy_coords)):
                if copy_coords[n].is_nonzero():
                    copy_coords[n] = copy_coords[n].rotate_arround(
                        z_ang, x_ax)
                x_prime = x_prime.rotate_arround(z_ang, x_ax)
                y_prime = y_prime.rotate_arround(z_ang, x_ax)
                z_prime = z_prime.rotate_arround(z_ang, x_ax)

    # z axes are now alligned, we now align x axes by rotation arround z
    if not np.allclose(x_prime, x_ax):
        x_ang = x_prime.angle_between(x_ax)
        crossp = vector(np.cross(x_prime, x_ax))
        if np.allclose(crossp, vector([0, 0, 0])):
            x_ang = np.pi
        elif crossp[2] < 0:
            x_ang = -x_ang
        for n in range(len(copy_coords)):
            if copy_coords[n].is_nonzero():
                copy_coords[n] = copy_coords[n].rotate_arround(
                    x_ang, z_ax)
            x_prime = x_prime.rotate_arround(x_ang, x_ax)
            y_prime = y_prime.rotate_arround(x_ang, x_ax)
            z_prime = z_prime.rotate_arround(x_ang, x_ax)

    # if all the dot products were 1, we are done, otherwise we have switched
    # to a left handed coordinate system and need to invert Y
    d_prods = [np.dot(a, b) for a, b in zip([x_prime, y_prime, z_prime],
                                            [x_ax, y_ax, z_ax])]
    if not all(np.allclose(a, 1) for a in d_prods):
        for n in range(len(copy_coords)):
            copy_coords[n][1] = -copy_coords[n][1]

    # lasty divide by the length of the new basis vectors if they weren't unit
    # vectors

    if not np.allclose(x_scale, 1):
        for n in range(len(copy_coords)):
            copy_coords[n][0] = copy_coords[n][0]/x_scale
    if not np.allclose(y_scale, 1):
        for n in range(len(copy_coords)):
            copy_coords[n][1] = copy_coords[n][1]/y_scale
    if not np.allclose(z_scale, 1):
        for n in range(len(copy_coords)):
            copy_coords[n][2] = copy_coords[n][2]/z_scale

    return copy_coords


def cylind2cart(coordinates, i_hat, j_hat, k_hat):
    """


    Parameters
    ----------
    coordinates : array of vectors
        an n by 3 array of 3d vectors given in cylindrical coordinates.
        Follows the convention z(height), theta (angle), r(radius)
    i_hat : vector
        Cartesian unit basis vector. conventionally the x axis.
    j_hat : vector
        Cartesian unit basis vector. conventionally the y axis.
    k_hat : TYPE
        Cartesian unit basis vector. conventionally the z axis.

    Returns
    -------
    cartesian : array of vectors
        returns the cartesian version of the inputted cylindrical coordinates.

    """
    cartesian = coordinates.copy()
    for i in range(len(coordinates)):
        z = coordinates[i][0]
        theta = coordinates[i][1]
        r = coordinates[i][2]
        cartesian[i][0] = r*np.cos(theta)
        cartesian[i][1] = r*np.sin(theta)
        cartesian[i][2] = z
    cartesian = change_basis(cartesian, vector([0, 0, 0]), i_hat, j_hat, k_hat)
    return cartesian


def cart2cylind(coordinates, center_vector, axial_vector, radial_vector):
    """
    converts standard cartesian coordinates to cylindrical coordinates
    (z, theta, r) recentered on center_vector with angles relative to
    radial_vector and height along axial_vector.

    Parameters
    ----------
    coordinates : TYPE
        DESCRIPTION.
    center_vector : vector
        cartesian spatial vector that will be subtracted from points
        before conversion to cylindrical to recenter on this point.
    axial_vector : vector
        vector in cartesian coordinates that will form the long axis
    radial_vector : vector
        vector in cartesian coordinates that will form the radial axis

    Returns
    -------
    coordinates.

    """
    # print(f'cart2cylind : axial_vector {axial_vector}')
    # print(f'cart2cylind : radial_vector {radial_vector}')
    if not axial_vector.is_nonzero() or not radial_vector.is_nonzero():
        raise ValueError('Basis vectors must be non-zero')
    if not axial_vector.is_orthogonal(radial_vector):
        raise ValueError('Basis vectors must be orthogonal')
    axial_vector = axial_vector.unitize()
    radial_vector = radial_vector.unitize()
    crossp = vector(np.cross(axial_vector, radial_vector))
    # print(type(coordinates), type(axial_vector), type(radial_vector))
    # print(coordinates, axial_vector, radial_vector, crossp)
    cylindrical_coordinates = change_basis(coordinates, center_vector,
                                           radial_vector, crossp,
                                           axial_vector)

    for n in range(len(cylindrical_coordinates)):

        q = cylindrical_coordinates[n][2]
        d = cylindrical_coordinates[n].project_onto_normal_plane(
            vector([0, 0, 1]))

        radius = d.get_length()
        if radius == 0:
            cylindrical_coordinates[n][1] = 0
        else:
            cylindrical_coordinates[n][1] = convert_angle(
                np.arctan2(d[1], d[0]))
        cylindrical_coordinates[n][2] = radius
        cylindrical_coordinates[n][0] = q
        # print(cylindrical_coordinates[n])

    # i = cylindrical_coordinates[:, 0]
    # j = cylindrical_coordinates[:, 1]
    # k = cylindrical_coordinates[:, 2]

    # plt.figure(figsize=(6, 6))  # Adjust figure size if needed

    # plt.polar(j, k)  # Plot polar coordinates

    return cylindrical_coordinates


def convert_angle(angle):
    normalized_angle = angle % (2 * np.pi)  # Wrap angle to range [0, 2*pi)
    if normalized_angle > np.pi:
        normalized_angle -= 2 * np.pi
    elif normalized_angle <= -np.pi:
        normalized_angle += 2 * np.pi
    return normalized_angle

# %% I/O functions


def import_pdb(fname=None):
    # we are making an numpy array filled with the x y and z for each atom in each
    # row doing it as a fixed array because this is much more memory efficient

    if fname is None:
        fname = '7k3g.pdb'
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
                        'phi': float('nan'), 'psi': float('nan'),
                        'omega': float('nan')
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


def plot_model(coordinates, fig=None, title='structure', plots_per_row=3):

    i = coordinates[:, 0]
    j = coordinates[:, 1]
    k = coordinates[:, 2]

    if fig is None:
        fig = plt.figure()

    num_subplots = len(fig.axes)

    # Calculate subplot dimensions
    rows = (num_subplots // plots_per_row) + 1
    cols = min(num_subplots % plots_per_row + 1, plots_per_row)

    # Add subplot
    ax = fig.add_subplot(rows, cols, num_subplots + 1, projection='3d')
    ax.set_box_aspect([1, 1, 1])

    v = ax.scatter3D(i, j, k, c=np.sqrt(i**2 + k**2))
    # ax.set_xlim(-40, 40)
    # ax.set_ylim(-40, 40)
    # ax.set_zlim(-40, 40)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.set_title(title)
    fig.tight_layout()

    return fig


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
                                residue.get_resname(),
                                'atom_name': a_tuple[4][0],
                                'phi': float('nan'), 'psi': float('nan')
                            }
                            # the last 2 are empty fields where phi &
                            # psi are assigned

                            info_table.append(a_dict)
                            n = n + 1
            coordinates = np.zeros(shape=(n, 3))
            coordinates = coordinates.view(vector)

        except:
            pass


def write_pdb(coordinates, info_table, outfile=None):
    """
    coordinates : Array of float 64
        the atomic coordinates of all atoms
    info_table : list
        list of dictionaries containing data for each coordinate
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


def ramachandran(coordinates, info_table):

    # !!! currently not working will fix later
    get_phi_psi(coordinates, info_table)
    alpha_carbons = select(info_table, atom_name='CA')
    phis = []
    psis = []
    res_nums = []
    for n, res in enumerate(alpha_carbons):
        phis.append(info_table[res]['phi'])
        psis.append(info_table[res]['psi'])
        res_nums.append(info_table[res]['res_num'])
    print(phis)
    return mpl.pyplot.scatter(phis, psis)


def select(info_table, **kwargs):

    # this function will return the indices in a list of all rows in the info
    # list who match the input arguments
    serial_numbers = []
    info = info_table
    db = kwargs.pop('debug_mode', False)

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

# %% getter/setter functions


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
            # if debug_mode:
            # print('get_phi_psi activated in mode 1. Searching based on'
            #       ' serial numbers.')
            ser_nums = kwargs.pop('ser_num')
            if isinstance(ser_nums, int):
                ser_nums = [ser_nums]
            alpha_sns = []
            for sn in ser_nums:

                # if it is not an alpha carbon, find the CA of the residue it
                # is on
                # print(sn)
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
            # print(alpha_sns)
            bb_ser_num_sets = get_dihedral_backbones(info_table, alpha_sns)

            # print(bb_ser_num_sets)
            # check if phi and psi are definable for this residue, then
            # pass backbone atom ser_nums to measure_dihedral to get them
            for group in bb_ser_num_sets:
                # print(f'group = {group}')
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
            # print(f'ca serial number are {ca_ser_nums}')
            bb_ser_num_sets = get_dihedral_backbones(info_table, ca_ser_nums)

            for group in bb_ser_num_sets:
                # print(f'group = {group}')
                if group[0] is not None:  # if there was a prev residue
                    vects = [coordinates[sn]
                             for sn in group[0:4]]  # pass first 4
                    phi = measure_dihedral(vects, debug_mode=True)
                    for sn in group[1:4]:
                        # print(sn)
                        # print('can assign phi')
                        info_table[sn]['phi'] = phi
                else:  # else phi is NaN
                    for sn in group[1:4]:
                        # print(sn)
                        info_table[sn]['phi'] = float('nan')
                if group[-1] is not None:  # if there is a next residue
                    # pass last 4, reversed
                    vects = [coordinates[sn] for sn in group[1:]]
                    psi = measure_dihedral(vects)
                    for sn in group[1:-1]:
                        # print(sn)
                        # print('can assign psi')
                        info_table[sn]['psi'] = psi
                else:  # else psi is NaN
                    for sn in group[1:-1]:
                        # print(sn)
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
            # print(f"find_rotation_angle: {target_angle - current_angle}"
            #       f" ({(target_angle - current_angle)*180/np.pi} degrees)")
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

    # print(f'the kwargs after processing for set_phi_psi are {kwargs}')
    kwargs['atom_name'] = 'CA'
    # if an atom name was specified ignore it
    # and search for 'CA'
    ca_ser_nums = select(info_table, **kwargs)
    # print(f'c alpha ser_nums are {ca_ser_nums}')

    for ndx, ca in enumerate(ca_ser_nums):

        backbone_sn_groups = get_dihedral_backbones(info_table, [ca])
        # print(f'the backbone atom indices for this are: {backbone_sn_groups}')
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
                    # print(rotation_segment)
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
                # print(
                #     f"omega for residue {next_res -1} = "
                # f"{info_table[s]['omega']}")


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
    # print(locals())
    def onion_method():
        print('--------------------------------------')
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
            # print(abs(v - s))
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
            # print(f'number of height index bins: {len(height_ndx_bins)}')
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
    #cylind_coords -= center_v

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
    # print(f'axis vector:     {axis_vector}')
    cutoff_distance = kwargs.get('cutoff_distance', 0.36)
    chainlist = kwargs.get('chainlist', None)
    max_tries = kwargs.get('max_tries', 20)
    center_sn = kwargs.get('center_sn', None)
    radial_sn = kwargs.get('radial_sn', None)

    symmetry_groups = kwargs.get('symmetry_groups', 'all')
    radius = kwargs.get('radius', None)
    multiplicity = kwargs.get('multiplicity', None)
    print(kwargs)
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
        # print( '__________________________________________')

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
    # print(f'anchored_residues are {anchored_residues}')
    # print(f'fixed_residues are {fixed_residues}')

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

    # print(f'after set comprehension, the free residues are\n{free_residues}')

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
        # print(f'the n anchor is {n_anchor}')
        c_anchor = max(anchored_residues)
        # print(f'the c anchor is {c_anchor}')
    # if an anchored segment exists, the free residues are split into those on
    # the n_terminal side of it (c_anchored) or the c terminal side
    # (n_anchored)

    c_anchored_segments = []
    n_anchored_segments = []
    # print(f'there are {len(sn_by_chain)} chains being randomized')

    # for each chain that needs to be checked internally, add a list of
    # its n terminal and c terminal segments serial numbers as a
    # sublist to keep track of associated chains

    # print(anchored_residues)
    for i in range(len(sn_by_chain)):

        # if there are no anchors treat that whole chain of free residues as
        # as being N anchored at the first residue
        if n_anchor is None and c_anchor is None:
            n_anchored_segments.append(free_residues[i])
        else:
            if n_anchor is not None:

                seg = [res for res in free_residues[i] if res < n_anchor]
                # print(f'seg is {seg}')
                # then put it in reverse order so that rotation iterates from C
                # towards n terminus

                # print(f' the {i}th set of c_anchored segments is {seg}')
                seg.sort(reverse=True)
                c_anchored_segments.append(seg)
            if c_anchor is not None:
                # repeat this for n anchored segments
                seg = [res for res in free_residues[i] if res > c_anchor]
                # print(f'this is an n-anchored segment {seg}')
                seg.sort(reverse=False)
                n_anchored_segments.append(seg)

    if len(c_anchored_segments) != len(n_anchored_segments) and \
            len(free_residues) != len(chainlist):
        raise ValueError('the number of lists in c_anchored_segments, '
                         'n_anchored_segments and the number of items in '
                         'chainlist must all be equal')

    # print(f'the chain in this selection are\n{chainlist}\n')
    # print(f'the C-anchored segments for each chain are {c_anchored_segments}')
    # print(f'the N-anchored segments for each chain are {n_anchored_segments}')
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
                        while (symmetry_clash(coordinates)
                               or check_internal_clash(coordinates,
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
                        while (symmetry_clash(coordinates)
                               or check_internal_clash(coordinates,
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
                        while (symmetry_clash(coordinates)
                               or check_internal_clash(coordinates,
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
                        while (symmetry_clash(coordinates)
                               or check_internal_clash(coordinates,
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


def orient(coordinates, info_table, center_sn, axis_sn, radial_sn=None,
           radial_angle=None, reference_vect=None, **kwargs):
    """
    Given a serial number for a central atom and a serial number for an atom
    to define a long axis, orient everything on that chain so that the vector
    from the center atom to the axis atom alligns with the reference vector.

    Parameters
    ----------
    coordinates : numpy ndarray of vectors
        the atomic coordinates of all atoms.
    info_table : list
        list of dictionaries containing data for each coordinate.
    center_sn : int
        serial number of the atom which is to be centered at [0 0 0].
    axis_sn : int
        serial number of the atom which will define the new + z axis relative
        to the center atom.
    radial_sn : int, optional
        serial number of the atom whose projection onto the normal plane of
        the desired z-axis will form a reference vector for rotation about the
        internal longitudinal axis. The default is None.
    radial_angle : float, optional
        if given, the angle about the z axis by which the radial reference
        vector should be rotated. The default is None.
    reference_vect : vector, optional
        If, given, an r3 vector to which the longitudinal axis should be
        alliged instead of the +z axis. The default is None.
....**kwargs : the only option for now is 'recenter' : (True | False) that
        determines whether center_sn should be at [0, 0, 0] in the returned
        structure (True) or if the centering done for the rotations should be
        undone to return it to its original position (after rotation) (False).
        Default is False.

    Returns
    -------
    coordinates : numpy ndarray of vectors.
        coordinates with orientation transformation applied

    points_to_center : vector
        if radial_sn was given this is the vector that points along the radial
        axis. else None is returned.

    """

    # if not specified do not recenter the final coordinates on center_sn
    recenter = kwargs.get('recenter', False)
    if not all(type(x) == int for x in (axis_sn, center_sn, radial_sn)):
        raise ValueError('axis_sn, center_sn, and radial_sn must each be a'
                         'single integer and a valid coordinate serial number')

    x_ax = vector([1, 0, 0])
    y_ax = vector([0, 1, 0])
    z_ax = vector([0, 0, 1])

    # check for correct input
    if radial_sn is not None:
        ref_points = [info_table[center_sn], info_table[axis_sn],
                      info_table[radial_sn]]
    else:
        ref_points = [info_table[center_sn], info_table[axis_sn]]

    ch = ref_points[0]['chain']
    a = all(x['chain'] == ch for x in ref_points)
    mod = ref_points[0]['model_name']
    b = all(x['model_name'] == mod for x in ref_points)
    sub = ref_points[0]['submodel']
    c = all(x['submodel'] == sub for x in ref_points)
    if not a or not b or not c:
        raise ValueError('All coordinate points for orientation must be from'
                         'the same chain, model, and submodel')

    else:

        ser_nums = select(info_table, chain=ch, model_name=mod, submodel=sub)
        # the atom to be made [0,0,0]
        center = coordinates[center_sn].copy()

        coordinates -= center
        # atom to be set to [0,0,z]
        axis_vector = coordinates[axis_sn].copy()
        if np.array_equal(center, axis_vector):
            raise ValueError('axis vector and center vector can not be the'
                             ' same')
        # plot_model(coordinates, None, 'after_centering')
        # find angle between desired axis and yz plane

        xy = axis_vector.project_onto_normal_plane(z_ax)
        crossp = vector(np.cross(xy, y_ax))
        if not xy.is_nonzero():
            yz_ang = 0

        else:
            yz_ang = y_ax.angle_between(xy)

        # if the Z component of the cross product is negative reverse direction
        if crossp[2] < 0:
            yz_ang = -yz_ang

        # rotate everything around the z axis by that angle
        for n in ser_nums:
            if n != center_sn:
                coordinates[n] = coordinates[n].rotate_arround(
                    yz_ang, z_ax)

        # our reference is now in the yz plane, find the angle to +z
        yz = coordinates[axis_sn]
        z_ang = z_ax.angle_between(yz)

        # if the X component of the cross product is negative reverse direction
        crossp = np.cross(yz, z_ax)
        if crossp[0] < 0:
            z_ang = -z_ang

        for n in ser_nums:

            if n != center_sn:
                coordinates[n] = coordinates[n].rotate_arround(z_ang, x_ax)

        # if you specified somnething other than the z axis to align
        # the longitudinal axis along
        if reference_vect is not None:
            crossp = vector(np.cross(z_ax, reference_vect))
            offset_ang = z_ax.angle_between(reference_vect)

            if not np.allclose(crossp, vector([0, 0, 0])):
                test = z_ax.rotate_arround(offset_ang, crossp)
                if np.dot(test, reference_vect) < 0:
                    offset_ang = -offset_ang
                for n in ser_nums:
                    coordinates[n] = coordinates[n].rotate_arround(offset_ang,
                                                                   crossp)
            # plot_model(coordinates, None, 'aligned to rot_ax')
        else:
            reference_vect = z_ax
            crossp = vector([0, 0, 0])
        # the radial projection vector will point to the center of
        # the subunit assembly, when the subunit is rotated by the given
        # angle arround an internal axis. If we want the vector from
        # center atom and the vector to the center of subunit
        # rotation to be e.g. +30 degrees, we specify pi/6 as the argument for
        # radial angle

        if radial_sn is not None and radial_sn != center_sn \
                and radial_sn != axis_sn:

            radial_proj = coordinates[radial_sn].project_onto_normal_plane(
                reference_vect)

            if radial_proj.is_nonzero():
                if radial_angle is not None:
                    if radial_angle != 0:
                        points_to_center = radial_proj.rotate_arround(
                            -radial_angle, reference_vect).copy()
                        points_to_center = points_to_center.unitize()
                    else:
                        points_to_center = radial_proj.unitize()
                else:
                    points_to_center = None
            else:
                raise ValueError('the radial vector has no component vector'
                                 'orthogonal to the axial vector')
        else:
            points_to_center = None
        # plot_model(coordinates, None,
        #            'after rotation arround internal longitudinal axis')
        if not recenter:

            coordinates += center

    return coordinates, points_to_center


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

    # print(f'central_atom is at {coordinates[center_sn]} before translation')
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
                # print(f'n : {n}')
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


def get_struct_orientation(fname):

    # !!! NOTE !!!!! this is not intended as a general use funtion yet
    # only useful for the pdb 7k3g substructure 1
    # The following will perfectly allign a single helix to mei hong's model
    # so that multimerized versions will overlay exactly in the TMD region
    coordinates, info_table = import_pdb(fname)
    center_sn = select(info_table, res_num=23, chain='A', atom_name='CA')[0]
    axis_sn = select(info_table, res_num=13, chain='A', atom_name='CA')[0]
    radial_sn = select(info_table, res_num=21, chain='A', atom_name='O')[0]

    x_ax = vector([1, 0, 0])
    y_ax = vector([0, 1, 0])
    z_ax = vector([0, 0, 1])

    # for some reason it mei hongs structure is tilted over 90 degrees
    # for n in range(len(coordinates)):
    #     coordinates[n] = coordinates[n].rotate_arround(np.pi/2, y_ax)
    points = [coordinates[sn]
              for sn in select(info_table, atom_name='CA', res_num=23)]
    points = np.asarray(points).copy()
    center = sum(points)/5
    ref_point = coordinates[center_sn].copy()
    coordinates = vector(coordinates)

    # plot_model(coordinates[select(info_table, chain='A')])
    points_to_center = -(ref_point-center)
    radius = points_to_center.get_length()
    coordinates -= ref_point
    rad_vect = coordinates[radial_sn].project_onto_normal_plane(z_ax).copy()

    radial_angle = points_to_center.angle_between(rad_vect)
    axis = coordinates[axis_sn].copy()

    cross_prod = vector(np.cross(axis, z_ax))
    axial_offset_ang = axis.angle_between(z_ax)
    # plot_model(coordinates, title='5X29')
    # print(locals())
    return radius, radial_angle, cross_prod, axial_offset_ang

# <codecell>
###### %% SANDBOX AREA FOR TESTING ##########################
# this function loads mei hongs structure and measures the exact radius
# the exact offset of the axial angle from the z axis, and the exact radial
# angle as defined by our set of reference atoms, so that we can apply these
# adjustments to our structures
coordinates, info_table = import_pdb('5x29_S2E_completed.pdb')
center_sn = select(info_table, res_num=23,
                   chain='A', atom_name='CA')[0]
axis_sn = select(info_table, res_num=13, chain='A', atom_name='CA')[0]
radial_sn = select(info_table, res_num=21, chain='A', atom_name='O')[0]
radius, radial_angle, cross_prod, offset_ang = get_struct_orientation(
    '5X29.pdb')
# print(f'ax:      {ax}')
axis_vector = vector([0, 0, 1]).rotate_arround(-offset_ang, cross_prod)
print(
    f'\n\nfor the reference model;\nradius: {radius}\nradial_angle: {radial_angle}\ncross_prod: {cross_prod}\noffset_angle: {offset_ang}\naxis vector: {axis_vector}\n\n')


# plot_model(coordinates, title='reoriented')
multiplicity = 5


for w in range(1):
    coordinates, info_table = import_pdb(
        '5x29_S2E_completed.pdb')
    center_sn = select(info_table, res_num=23,
                        chain='A', atom_name='CA')[0]
    axis_sn = select(info_table, res_num=13, chain='A', atom_name='CA')[0]
    radial_sn = select(info_table, res_num=21, chain='A', atom_name='O')[0]
    # plot_model(coordinates)
    # for i in range(0, 3):

    # fig = plot_model(coordinates, None, 'untransformed hong monomer')
    # # orient the structure to align reference with z
    coordinates, points_to_center = orient(coordinates, info_table, center_sn,
                                            axis_sn, radial_sn, np.pi, recenter=False,
                                            )
    print(f'points_to_center:                    {points_to_center}')

    #plot_model(coordinates)
    if points_to_center is not None:
        # be sure to rotate the radial vector for each subunit in the same way the principle axis was
        points_to_center = points_to_center.rotate_arround(-offset_ang, cross_prod)
        radial_v = points_to_center.copy()*radius
        multimeric_center = coordinates[center_sn] + radial_v
        print(f'multimeric_center:                    {multimeric_center}')
        # coordinates -= multimeric_center
    else:
        print('a good radial coordinate is needed to adjust the internal radial'
              'angles of each subunit when multimerizing')

    print(f'center coordinate: {coordinates[center_sn]}')

    print(f'radius: {radius}')
    # bundle all of these together to be passed along
    rot_sym_vars = {'multiplicity' : multiplicity, 'radial_v' : radial_v,
                    'axis_vector' : axis_vector,
                  'radius' : radius, 'center_sn' : center_sn,
                  'radial_sn' : radial_sn,
                  'multimeric_center' : multimeric_center}
    print(f'\naxis vector {axis_vector}')
    print(f'\nradial_vector {radial_v}')
    print(points_to_center.dot(multimeric_center))
    print(axis_vector.is_orthogonal(radial_v))



#     # the inclusion of rot_sym_vars allows the machine to check symmetry as it is randomizing
    coordinates, info_table = random_backbone(coordinates, info_table, 1,
                                              "residues.txt", max_tries=100,
                                              **rot_sym_vars)

#     # check_method='definitive', sym_clash_check=True, rot_sym_vars=rot_sym_vars)

#     output1 = axial_symmetry(coordinates, info_table,
#                               multiplicity, radius,
#                               center_sn,
#                               axis_sn, radial_sn,
#                               radial_angle=radial_angle, threshold=7)
#     if output1 is not None:
#         coordinates, info_table = output1
#         plot_model(coordinates, title=f'2 pi / {w} radial angle')
#     #     write_pdb(coordinates, info_table,
#     #               f'randomized_pentamers/randomized_pentamer_{w}')

#     # else:
#     #     print(
#     #         f'you cannot have {multiplicity}-fold radial symmetry for model {w}')

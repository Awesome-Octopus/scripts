#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 17:51:39 2024

@author: andrew
"""
import numpy as np
from vector_class import vector


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

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 22:15:30 2024

@author: andrew
"""
import numpy as np
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D
import MDAnalysis as mda
import argparse
from matplotlib import pyplot as plt


#parser = argparse.ArgumentParser()

# the input .gro file
# parser.add_argument('--gro', type=str)

# input trajectory file
# parser.add_argument('--xtc', type=str)

# beginning residue of the helix
# parser.add_argument('-b', type=int)
begin_ctd = 55
# ending residue of the helix
# parser.add_argument('-e', type=int)
end_ctd = 63

# for parsing the flag of whether or not to align to +z or to a second helix
#parser.add_argument('-z', type=int, Optional=True)

# args = parser.parse_args()


begin_tmd = 12
end_tmd = 38
# ################# make it so that ref_axis is assigned based on flag
#  if <you got the argument -z followed by a single number>:
    #    process in for loop down below
    #else if you got a vector:
    #    import directly from arg string
    # else:
    #   set ref_axis to z

# temporarily hard coded to be the first helix given
ref_axis=1


# trj = mda.Universe(args.gro, args.xtc)


###############################

# desired syntax:
#   python helical_angles.py -gro input_gro_file -xtc input_xtc_file
#   -1 start:end -2 start:end -3 start:end -z



##############################################################
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

trj = mda.Universe('step7_production.gro', 'step7_production.xtc')


# a list of slices each containing a protein chain
all_atoms = trj.select_atoms('name *')
protein = trj.select_atoms('name BB or (name SC?)')
chains = [[protein.indices[0]]] if len(protein.indices) > 0 else []

# check if every protein atom is the same resnum or is continguous with the resnum
# of the atom before it, otherwise start a new chain
for curr_ndx in protein.indices[1:]:
    prev_ndx = chains[-1][-1]
    if all_atoms.resnums[curr_ndx] in {all_atoms.resnums[prev_ndx],  all_atoms.resnums[prev_ndx] +1}:
        chains[-1].append(curr_ndx)
    else:
        chains.append([curr_ndx])

theta = np.ndarray((len(trj.trajectory), len(chains)))
phi = np.ndarray((len(trj.trajectory), len(chains)))

# iterate over each frame
for i, ts in enumerate(trj.trajectory):

    #values calculated for each helix this frame


    com = all_atoms.select_atoms(f'resnum {begin_tmd}:{end_tmd}')
    com = vector(com.center_of_geometry(compound='group'))
    com[2] = 0 #project the center on to the xy plane

    z = vector([0, 0, 1])

    for j, chain in enumerate(chains):

        f = all_atoms[chain]
        helix1_bb = f.select_atoms(f'resnum {begin_ctd}:{end_ctd} and name BB')
        helix2_bb = f.select_atoms(f'resnum {begin_tmd}:{end_tmd} and name BB')

        helix1_bb_coords = helix1_bb.positions

        helix1_start = helix1_bb.select_atoms(f'resnum {begin_ctd}')

        radial_vect = vector(helix1_start.positions - com)[0]

        helix1_pca = PCA(n_components=3)
        helix1_pca.fit(helix1_bb_coords)


        helix1_vect = vector(helix1_pca.components_[0])


        if ref_axis is None:

            xy_helix1 = vector([helix1_vect[0], helix1_vect[1], 0]).unitize()
            xy_radial = vector([radial_vect[0], radial_vect[1], 0]).unitize()
            theta[i][j] = xy_radial.angle_between(xy_helix1)
            phi[i][j] = helix1_vect.angle_between(z)


        elif type(ref_axis) == int: # if a single number is given the ref axis is
        # the helical axis with that index, and we want helical angles relative to each other

            # take the exact same code up here, but first rotate the helix 1 vector
            # such that the same rotation applied to helical axis 2 (the tmd) makes
            # it line up onto the +z axis

            helix2_bb_coords = helix2_bb.positions
            ref_helix_pca = PCA(n_components=3)
            ref_helix_pca.fit(helix2_bb_coords)

            # get vector for helix 2 by pca
            ref_helix_vect = vector(ref_helix_pca.components_[0])
            crossp = vector(np.cross(ref_helix_vect, z))

            offset_angle = ref_helix_vect.angle_between(z)

            offset_helix1 = helix1_vect.rotate_arround(offset_angle, crossp)

            offset_rad_vect = radial_vect.rotate_arround(offset_angle, crossp)


            phi[i][j] = offset_helix1.angle_between(z)
            theta[i][j] = offset_rad_vect.angle_between(offset_helix1)

        else: #if we give it an arbitrary axis vector
            ref_axis = vector(np.array(ref_axis)) # convert input to vector


            pass
# *   <-----------  center of mass
#  \
#   \    THETA = helix angle relative to the radial axis
#    \
# *---* <---------xy_projection of the vector from the center of mass to
# ^                the starting atom on the helix n-term
# |
# xy projection of the direction of the helix from n to c term







# Plot the original data and the transformed data
fig = plt.figure()

for i in range(len(chains)):
    ax = fig.add_subplot(3,2,i+1)
    indices = np.arange(len(theta[:,i]))
    sca = ax.scatter(theta[:, i], phi[:, i], c=indices, cmap='viridis', label=f'chain {i+1}')

    plt.xlabel(r'$\theta$ (radians)')
    plt.ylabel(r'$\phi$ (radians)')
    ax.legend()

plt.colorbar(sca, label='timestep')
plt.tight_layout()
plt.show()

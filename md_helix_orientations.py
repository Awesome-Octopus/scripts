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

###############################

# syntax:
#   python helical_angles.py --gro input_gro_file.gro --xtc input_xtc_file.xtc
#   --helix start:end --ref start:end --roll #

# ##############################################################
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
        return shadow

def parse_helix_arg(arg):
    # Check if the argument has the correct format
    if ':' not in arg:
        raise argparse.ArgumentTypeError("format for --helix argument is \"start:end\"")
        # Split the second part by ':'
    start, end = map(int, arg.split(':'))
    return (start, end)

def parse_ref_arg(arg):
    # remove whitespace if present
    arg = arg.replace(" ", "")
    # if a vector is given as imput e.g. [1,2,3]
    if arg[0] == '[' and arg[-1] == ']':
        try:
            return vector(arg)
        except:
            raise ValueError('vector inputs for --ref flag should be in the form [x, y, z]')
    else:
        try:
            start, end = parse_helix_arg(arg)
            return start, end
        except:
            return arg


parser = argparse.ArgumentParser()

parser.add_argument('-g', '--gro', type=str, required=True, help='the input .gro file')

# input trajectory file
parser.add_argument('-x', '--xtc', type=str, required=True, help='input xtc trajectory file')

# collect all arguments for the beginning and end of helices
parser.add_argument('--helix', type=parse_helix_arg, required=True, help='beginning and end residues of the helix in the form begin:end')

# for parsing the flag of whether or not to align to +z or to a second helix
parser.add_argument('-r', '--ref', type=parse_ref_arg, required=True, help='beginning and end residues of the symmetry reference helix in the form begin:end, or \'z\' to use the z axis as reference')


parser.add_argument('--roll', type=int, required=False, help='the residue on the helix whose c alpha c beta bond or BB to SC1 bond will serve as the metric for roll angle')


parser.add_argument('--cg', action='store_true', required=False, help='if cg flag is passed treat the structure as martini cg atoms, else, all atom')

args = parser.parse_args()

############################################

#trj = mda.Universe("test_pentamer.pdb")

trj = mda.Universe(args.gro, args.xtc)

z = vector([0, 0, 1])

# a list of slices each containing a protein chain
all_atoms = trj.select_atoms('name *')
if args.cg:
    protein = trj.select_atoms('name BB or (name SC?)')
else:
    protein = trj.select_atoms('protein')
chains = [[protein.indices[0]]] if len(protein.indices) > 0 else []

# check if every protein atom is the same resnum or is continguous with the resnum
# of the atom before it, otherwise start a new chain
for curr_ndx in protein.indices[1:]:
    prev_ndx = chains[-1][-1]
    if all_atoms.resnums[curr_ndx] in {all_atoms.resnums[prev_ndx],  all_atoms.resnums[prev_ndx] +1}:
        chains[-1].append(curr_ndx)
    else:
        chains.append([curr_ndx])

pitch = np.ndarray((len(trj.trajectory), len(chains)))
yaw = np.ndarray((len(trj.trajectory), len(chains)))
roll = np.ndarray((len(trj.trajectory), len(chains)))

# iterate over each frame
for i, ts in enumerate(trj.trajectory):

    # values calculated for the reference helices this frame if one was given                            com_xy
    if type(args.ref) == tuple:

        # the atoms of all reference helices in all monomers
        if args.cg:
            ref_helix_set = all_atoms.select_atoms(f'resnum {args.ref[0]}:{args.ref[1]} and name BB')
        else:
            ref_helix_set = all_atoms.select_atoms(f'resnum {args.ref[0]}:{args.ref[1]} and backbone')
        # the center of the mass of the set of all reference helices
        try:
            com = vector(ref_helix_set.center_of_geometry(compound='group'))
        except:
            raise ValueError('cannot find center of geometry. (Did you forget --cg flag for a martini topology?)')

        # this is the center projected on to the xy axis
        xy_com = com.copy()
        xy_com[2] =0


        # the axis of symmetry arround which the reference helices are best arranged
        # ^
        # |
        # for each residue across all monomers along the helix we want to define group symmetry
        symmetry_avg = []
        for res in range(int(args.ref[0]), int(args.ref[1])+1):
            # select that residue across all chains
            curr_res = ref_helix_set.select_atoms(f'resnum {res}')
            # add its average to a list
            symmetry_avg.append(vector(curr_res.center_of_geometry(compound='group')))

        symmetry_axis_pca = PCA(n_components=3)
        symmetry_axis_pca.fit(symmetry_avg)

        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # NOTE: this only works if the helices form an overall cylindrical
        # shape that is longer than wide it should be generalized to find
        # the overall axis of symetry between subunits
        symmetry_axis = vector(symmetry_axis_pca.components_[0])



    for j, chain in enumerate(chains):

        f = all_atoms[chain]

        # SELECT THAT ATOMS IN THE tesdt helix on that chain
        if args.cg:
            helix_bb = f.select_atoms(f'resnum {args.helix[0]}:{args.helix[1]} and name BB')
        else:
            helix_bb = f.select_atoms(f'resnum {args.helix[0]}:{args.helix[1]} and backbone')

        helix_bb_coords = helix_bb.positions
        e = helix_bb.select_atoms(f'resnum {args.helix[0]}')
        helix_com = vector(e.center_of_geometry(compound='group'))

        radial_vect = helix_com - com

        helix_pca = PCA(n_components=3)
        helix_pca.fit(helix_bb_coords)


        helix_vect = vector(helix_pca.components_[0])

        if args.roll: # calculate roll if a resnum was given as reference

            hel_vect_cross_z = vector(np.cross(helix_vect, z))
            # print(hel_vect_cross_z)
            # create a vector othogonal to the helical axis, that is the shortest
            # distance to the xy plane
            ortho_hel_vect = helix_vect.rotate_arround(np.pi/2, hel_vect_cross_z)
            # print(f'ortho_hel_vect {ortho_hel_vect}')
            # the positions of the backbone and the side chain of the residue
            # specified for roll form a reference vector to calculate roll

            if args.cg:
                bb_point = f.select_atoms(f'resnum {args.roll} and name BB')
                sc_point = f.select_atoms(f'resnum {args.roll} and name SC1')
            else:
                bb_point = f.select_atoms(f'resnum {args.roll} and name CA')
                sc_point = f.select_atoms(f'resnum {args.roll} and name CB')

            roll_vect = vector(sc_point.positions) - vector(bb_point.positions)

            roll_vect = roll_vect[0].project_onto_normal_plane(helix_vect)
            roll[i][j] = roll_vect.angle_between(ortho_hel_vect)
            print(roll[i][j])

        # select the atoms for the reference helix of that chain if given
        if type(args.ref) == tuple:

            if args.cg:
                ref_helix_bb = f.select_atoms(f'resnum {args.ref[0]}:{args.ref[1]} and name BB')
            else:
                ref_helix_bb = f.select_atoms(f'resnum {args.ref[0]}:{args.ref[1]} and backbone')

            # the helical axis with that index, and we want helical angles relative to each other

            # take the exact same code up here, but first rotate the helix 1 vector
            # such that the same rotation applied to helical axis 2 (the tmd) makes
            # it line up onto the +z axis

            ref_helix_bb_coords = ref_helix_bb.positions
            ref_helix_pca = PCA(n_components=3)
            ref_helix_pca.fit(ref_helix_bb_coords)

            # get vector for helix 2 by pca
            ref_helix_vect = vector(ref_helix_pca.components_[0])
            crossp = vector(np.cross(symmetry_axis, z))

            offset_angle = symmetry_axis.angle_between(z)

            offset_helix = helix_vect.rotate_arround(offset_angle, crossp)
            offset_ref_helix = ref_helix_vect.rotate_arround(offset_angle, crossp)
            offset_rad_vect = radial_vect.rotate_arround(offset_angle, crossp)


            pitch[i][j] = offset_helix.angle_between(z)
            yaw[i][j] = offset_rad_vect.angle_between(offset_helix)
            # print('pitch:\n',pitch[i][j],'\nyaw\n', yaw[i][j])


        # if you wanted to take the angle relative to the standard +z axis
        elif type(args.ref) == str and args.ref in 'xyz':

            if args.ref == 'z':
                xy_helix = vector([helix_vect[0], helix_vect[1], 0])
                xy_radial = vector([radial_vect[0], radial_vect[1], 0])
                yaw[i][j] = xy_radial.angle_between(xy_helix)
                pitch[i][j] = helix_vect.angle_between(z)

            elif args.ref == 'x':
                pass

            elif args.ref == 'y':
                pass

        # !!! ADD LATER
        # support for an arbitrary vector as the reference vector rather than
        # derive it from the average axis of symmetry for the reference helices
        elif type(args.ref) == vector: #if we give it an arbitrary axis vector
            pass

# *   <-----------  center of mass
#  \
#   \    THETA = helix angle relative to the radial axis
#    \
# *---* <---------xy_projection of the vector from the center of mass to
# ^                the starting atom on the helix n-term
# |
# xy projection of the direction of the helix from n to c term



fig = plt.figure()

for j in range(len(chains)):
    ax = fig.add_subplot(3,2,j+1)
    indices = np.arange(len(yaw[:,j]))
    sca = ax.scatter(yaw[:, j], pitch[:, j], c=indices, cmap='viridis',
                      label=f'chain {j+1}')

    plt.xlabel(r'yaw (radians)')
    plt.ylabel(r'pitch (radians)')
    ax.legend()

plt.colorbar(sca, label='timestep')
plt.tight_layout()
plt.show()

# make a separate plot for roll if the flag was given
if args.roll:
    roll_fig = plt.figure()
    for j in range(len(chains)):

        ax = roll_fig.add_subplot(3,2,j+1)
        indices = np.arange(len(roll[:,j]))
        sca = ax.scatter(indices, roll[:, j], label=f'chain {j+1}')

    plt.tight_layout()
    plt.show()
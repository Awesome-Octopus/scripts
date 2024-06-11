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
#   python helical_angles.py --struct input_gro_file.struct --traj input_traj_file.traj
#   --helix start:end --ref start:end --roll # --plot filename.png

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
        raise argparse.ArgumentTypeError("format for --helix argument is \"Nterm:Cterm\"")
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

# <codecell> argument parsing
parser = argparse.ArgumentParser()

# name for output plot
parser.add_argument('-p', '--plot', type=str, required=False, help='filename for output plot')



# input .gro or .psf file
parser.add_argument('-s', '--struct', type=str, required=True, help='the input .gro or .psf file')

# input trajectory file
parser.add_argument('-t', '--traj', type=str, required=True, help='input .xtc or .dcd trajectory file')

# collect all arguments for the beginning and end of helices
parser.add_argument('--helix', type=parse_helix_arg, required=True, help='beginning and end residues of the helix in the form begin:end')

# for parsing the flag of whether or not to align to +z or to a second helix
parser.add_argument('-r', '--ref', type=parse_ref_arg, required=True, help='beginning and end residues of the symmetry reference helix in the form begin:end, or \'z\' to use the z axis as reference')

parser.add_argument('--roll', type=int, required=False, help='the residue on the helix whose c alpha c beta bond or BB to SC1 bond will serve as the metric for roll angle')

parser.add_argument('--cg', action='store_true', required=False, help='if cg flag is passed treat the structure as martini cg atoms, else, all atom')

args = parser.parse_args()

############################################
# <codecell> structure parsing
#trj = mda.Universe("test_pentamer.pdb")

trj = mda.Universe(args.struct, args.traj)

z = vector([0, 0, 1])

# a list of slices each containing a protein chain
all_atoms = trj.select_atoms('name *')

# select protein atoms for cg models if flag set
if args.cg:
    if '.gro' in args.struct:
        protein = trj.select_atoms('name BB or (name SC?)')
        if not protein:
            raise ValueError('there appears to be no coarse grained protein beads present (did you mean for an all-atom structure?)')
    elif '.psf' in args.struct:
        protein = trj.select_atoms('name BAS or (name SI?)')
        if not protein:
            raise ValueError('there appears to be no coarse grained protein beads present (did you mean for an all-atom structure?)')
# otherwise, for all atom models
else:
    protein = trj.select_atoms('protein')
    if not protein:
        raise ValueError('there appears to be no protein atoms present (use --cg for a coarse grained stucture)')
chains = [[protein.indices[0]]] if len(protein.indices) > 0 else []

# check if every protein atom is the same resnum or is continguous with the resnum
# of the atom before it, otherwise start a new chain
for curr_ndx in protein.indices[1:]:
    prev_ndx = chains[-1][-1]
    if all_atoms.resnums[curr_ndx] in {all_atoms.resnums[prev_ndx],  all_atoms.resnums[prev_ndx] +1}:
        chains[-1].append(curr_ndx)
    else:
        chains.append([curr_ndx])

# initialize storage arrays
pitch = np.ndarray((len(trj.trajectory), len(chains)))
yaw = np.ndarray((len(trj.trajectory), len(chains)))
roll = np.ndarray((len(trj.trajectory), len(chains)))


        # the atoms of all reference helices in all monomers
if type(args.ref) == tuple:
    if args.cg:
        if '.gro' in args.struct:
            ref_helix_set = all_atoms.select_atoms(f'resnum {args.ref[0]}:{args.ref[1]} and name BB')
        else:
            ref_helix_set = all_atoms.select_atoms(f'resnum {args.ref[0]}:{args.ref[1]} and name BAS')
    else:
        ref_helix_set = all_atoms.select_atoms(f'resnum {args.ref[0]}:{args.ref[1]} and backbone')

# if you wanted to take the angle relative to a standard axis
elif type(args.ref) == str and args.ref in 'xyz':

    # we take the negative of the vector because by convention, N to C going
    # straight down the membrane normal is considered 0 degrees
    if args.ref == 'z':
        symmetry_axis = -z
    elif args.ref == 'x':
        symmetry_axis = vector([-1,0,0])

    elif args.ref == 'y':
        symmetry_axis = vector([0,-1,0])

# !!! ADD LATER
# support for an arbitrary vector as the reference vector rather than
# derive it from the average axis of symmetry for the reference helices
elif type(args.ref) == vector: #if we give it an arbitrary axis vector
    symmetry_axis = args.ref

##################################################
# <codecell> frame iteration
# iterate over each frame
for i, ts in enumerate(trj.trajectory):

    # if not given as a static vector calculate the axis of symmetry at each frame
    if type(args.ref) == tuple:

        # the center of the mass of the set of all reference helices
        try:
            com = vector(ref_helix_set.center_of_geometry(compound='group'))
        except:
            raise ValueError(f'cannot find reference center of geometry for frame {i}')

        # for each residue across all monomers along the helix we want to
        # define group symmetry
        symmetry_avg = []
        for res in range(int(args.ref[0]), int(args.ref[1])+1):

            # select that residue across all chains
            curr_res = ref_helix_set.select_atoms(f'resnum {res}')
            if not curr_res:
                raise ValueError(f'cannot find atoms matching residue number {res} for the reference helix. It may be because the numbering is wrong, or the --cg flag is set incorrectly.')
            # and add its average to a list
            symmetry_avg.append(vector(curr_res.center_of_geometry(compound='group')))


        symmetry_axis_pca = PCA(n_components=3)
        symmetry_axis_pca.fit(symmetry_avg)

        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # NOTE: this only works if the helices form an overall cylindrical
        # shape that is longer than wide it should be generalized to find
        # the overall axis of symetry between subunits
        symmetry_axis = vector(symmetry_axis_pca.components_[0])

        # since PCA doesnt know have a directionaly to the fit,
        # the vector might be pointing in the wrong direction.
        # we can check this with dot product
        if symmetry_axis.dot(z) < 0:
            symmetry_axis = -symmetry_axis

    else:

        com = vector(protein.center_of_geometry(compound='group'))

    for j, chain in enumerate(chains):

        f = all_atoms[chain]

        # SELECT the ATOMS IN THE test helix on that chain
        if args.cg:
            if '.gro' in args.struct:
                helix_bb = f.select_atoms(f'resnum {args.helix[0]}:{args.helix[1]} and name BB')

            else:
                helix_bb = f.select_atoms(f'resnum {args.helix[0]}:{args.helix[1]} and name BAS')
        else:
            helix_bb = f.select_atoms(f'resnum {args.helix[0]}:{args.helix[1]} and backbone')
        if not helix_bb:
            raise ValueError('No atoms match your input range for helix selection. Inspect your .gro or .psf file.')
        helix_bb_coords = helix_bb.positions
        helix_com = vector(helix_bb.center_of_geometry(compound='group'))


        # print(f'helix_com: {helix_com}')
        radial_vect = helix_com - com
        # consider only the vector component orthogonal to the symmetry axis
        radial_vect = radial_vect.project_onto_normal_plane(symmetry_axis)

        # helix from N to C
        helix_pca = PCA(n_components=3)
        helix_pca.fit(helix_bb_coords)

        helix_vect = vector(helix_pca.components_[0])

        # by getting the first bb position - last backbone position,
        # we get a test vector to tell if PCA is pointing in the
        # right direction or the opposite
        test_vect = vector(helix_bb_coords[1]) - \
            vector(helix_bb_coords[0])

        if helix_vect.dot(test_vect) < 0:
            helix_vect = -helix_vect


        # calculate roll if a resnum was given as reference
        if args.roll:

            hel_vect_cross = vector(np.cross(helix_vect, symmetry_axis))
            # print(hel_vect_cross)

            # create a vector othogonal to the helical axis,
            # that is the shortest
            # distance to the normal plane of the symmetry axis
            ortho_hel_vect = helix_vect.rotate_arround(np.pi/2, hel_vect_cross)

            # the positions of the backbone and the side chain of the residue
            # specified for the --roll argument form a reference vector
            # from which to calculate roll angle relative to the orthogonal
            # helix vector

            if args.cg: # for coarse grained
                if '.gro' in args.struct:  # for gromacs
                    bb_point = f.select_atoms(f'resnum {args.roll} and name BB')
                    sc_point = f.select_atoms(f'resnum {args.roll} and name SC1')
                if '.psf' in args.struct:  # for NAMD
                    bb_point = f.select_atoms(f'resnum {args.roll} and name BAS')
                    sc_point = f.select_atoms(f'resnum {args.roll} and name SI?')
                if not sc_point:
                    raise ValueError('the residue chosen to track helical roll does not have a side chain bead. Choose something other than Ala or Gly')
            else: # for all atom
                bb_point = f.select_atoms(f'resnum {args.roll} and name CA')
                sc_point = f.select_atoms(f'resnum {args.roll} and name CB')
                if not sc_point:
                    raise ValueError('The residue chosen to track helical roll does not have a C beta. Choose a residue other than Gly.')

            # points from backbone to sidechain (CA to CB in all atom)
            roll_vect = vector(sc_point.positions) - vector(bb_point.positions)


            roll_vect = roll_vect[0].project_onto_normal_plane(helix_vect)
            ang = roll_vect.angle_between(ortho_hel_vect)
            cp = vector(np.cross(ortho_hel_vect, roll_vect))
            if cp.dot(helix_vect) <= 0:
                ang = -ang
            roll[i][j] = ang
            #print(roll[i][j])

    ################### unneccesary for now ###############################

        # parse which atoms for the reference helix to take the symmetric average of
        # if type(args.ref) == tuple:

            # will be different if cg vs all atom
            # if args.cg:
                # ref_helix_bb = f.select_atoms(f'resnum {args.ref[0]}:{args.ref[1]} and name BB')
            # else:
                # ref_helix_bb = f.select_atoms(f'resnum {args.ref[0]}:{args.ref[1]} and backbone')
        # select the atoms for the reference helix of that chain if given


        # ref_helix_bb_coords = ref_helix_bb.positions
        # ref_helix_pca = PCA(n_components=3)
        # ref_helix_pca.fit(ref_helix_bb_coords)

        # get vector for helix 2 by pca
        # ref_helix_vect = vector(ref_helix_pca.components_[0])
        # crossp = vector(np.cross(symmetry_axis, z))

        # offset_angle = symmetry_axis.angle_between(z)

        # offset_helix = helix_vect.rotate_arround(offset_angle, crossp)
        # offset_ref_helix = ref_helix_vect.rotate_arround(offset_angle, crossp)
        # offset_rad_vect = radial_vect.rotate_arround(offset_angle, crossp)
    ###########################################################################

        helix_projection = helix_vect.project_onto_normal_plane(symmetry_axis)
        ang = helix_projection.angle_between(radial_vect)
        cp = vector(np.cross(radial_vect, helix_projection))

        # if the cross product of the radial vector and the projection of the helix vector
        # onto the normal plane of the symmetry axis is pointing away from the symmetry
        # axis, the angle is negative
        if cp.dot(symmetry_axis) <= 0:
            ang = -ang

        pitch[i][j] = helix_vect.angle_between(symmetry_axis)
        yaw[i][j] = ang
        # print('pitch:\n',pitch[i][j],'\nyaw\n', yaw[i][j])


# *   <-----------  center of mass
#  \
#   \    Yaw = helix angle relative to the radial axis
#    \
# *---* <---------xy_projection of the vector from the center of mass to
# ^                the center of helix mass
# |
# xy projection of the direction of the helix from n to c term


# <codecell> plotting
fig = plt.figure()

for j in range(len(chains)):
    ax = fig.add_subplot(3,2,j+1)
    indices = np.arange(len(yaw[:,j]))
    sca = ax.scatter(yaw[:, j]*180/np.pi, pitch[:, j]*180/np.pi, c=indices, cmap='viridis', label=f'chain {j+1}')

    plt.xlabel('yaw (degrees)')
    plt.ylabel('pitch (degrees)')
    plt.title(f'chain {j+1}')

plt.colorbar(sca, label='timestep')
plt.tight_layout()

if args.plot:
    if '.' in args.plot:
        name, ext = args.plot.split('.')
        file= name + '_.' + ext
        plt.savefig(file)
    else:
        plt.savefig(args.plot + '.png')
else:
    plt.show()

# make a separate plot for roll if the flag was given
if args.roll:
    roll_fig = plt.figure()
    for j in range(len(chains)):

        ax = roll_fig.add_subplot(3,2,j+1)
        indices = np.arange(len(roll[:,j]))
        sca = ax.scatter(indices, roll[:, j]*180/np.pi, label=f'chain {j+1}')
        plt.xlabel('timestep')
        plt.ylabel('roll (degrees)')
        plt.title(f'chain {j+1}')

    plt.tight_layout()
    if args.plot:
        if '.' in args.plot:
            name, ext = args.plot.split('.')
            roll_plot_file= name + '_roll.' + ext
            plt.savefig(roll_plot_file)
        else:
            plt.savefig(args.plot + '_roll.' + 'png')
    else:
        plt.show()

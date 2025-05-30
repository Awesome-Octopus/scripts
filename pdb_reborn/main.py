#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 01:36:49 2023

@author: andrew
"""


import numpy as np
import sys
from matplotlib import pyplot as plt
from vector_class import vector
from io_funcs import import_pdb, plot_model, write_pdb
from coord_funcs import cart2cylind, cylind2cart, change_basis, convert_angle
from getter_setter import orient, select, random_backbone, get_basis, clone_chain


# these define the residues whose coordinates are used to build interal
# axes for each chain
center_res = 1
axis_res = 1
radial_res = 1


chain = 'A'

# NOT YET IMPLEMENTED
# radial_offset_angle = 0
# tilt_offset_angle = 0
# lean_offset_angle = 0
# offset_matrix = None


coordinates, info_table = import_pdb('test_structs/test_case.pdb')
center_sn = select(info_table, res_num=center_res,
                   chain=chain, atom_name='CA')[0]
axis_sn = select(info_table, res_num=axis_res, chain=chain, atom_name='C')[0]
radial_sn = select(info_table, res_num=radial_res,
                   chain=chain, atom_name='N')[0]


query_struct, query_table = import_pdb('test_structs/coordinate_system.pdb')
query_chains = list(set(query_table[i]['chain']
                    for i in range(len(query_table))))

# CHANGE THIS TO CHANGE FROM A PENTAMER TO ANOTHER NUMBER IF DESIRED
multiplicity = len(query_chains)


query_basis = []

# get the basis for each chain in the symetric multimer you wish to average over
center_pts = []
query_center_sns = []
query_axis_sns = []
query_radial_sns = []
for ch in query_chains:

    c = select(query_table, res_num=center_res, chain=ch, atom_name='CA')[0]
    query_center_sns.append(c)
    center_pts.append(query_struct[c])
    a = select(query_table, res_num=axis_res, chain=ch, atom_name='C')[0]
    query_axis_sns.append(a)
    r = select(query_table, res_num=radial_res, chain=ch, atom_name='N')[0]
    query_radial_sns.append(r)
    query_basis.append(get_basis(query_struct, c, a, r))

# the average of all the atoms defining the center points of each chain
symmetry_center = vector(np.mean(np.asarray(center_pts), axis=0))
print(f'symmetry_center: {symmetry_center}')

group_basis = np.mean(np.asarray(query_basis), axis=0)
symmetry_axis = group_basis[0, :3].view(vector).unitize()
print(f'symmetry_axis: {symmetry_axis}')
# The radial vector for the group basis is somewhat arbitrary since with radial
# perfect symmetry it would average to the 0 vector. So it is assigned the
# coordinates that point from the group center towards the CENTER ATOM first subunit of
# the group.
radial_axis = (query_struct[query_center_sns[0]] - symmetry_center).project_onto_normal_plane(
    symmetry_axis).unitize()
print(f'radial_axis: {radial_axis}')
norm_axis = np.cross(symmetry_axis, radial_axis).view(vector)


# find the average distance of the center coordinate of each chain from the
# center of the group
query_radius = 0
for sn in select(query_table, res_num=center_res, atom_name='CA'):
    query_radius += symmetry_center.distance_between(query_struct[sn])
query_radius = query_radius/len(query_chains)

# how wide you want the radius of the circular arrangement of monomers to be
radius = query_radius
print(f'radius: {radius}')

# --------------DEBUG BLOCK--------------
# b = symmetry_center + symmetry_axis
# c = symmetry_center + radial_axis*radius
# out = np.vstack((symmetry_center, b, c))
# write_pdb(out, info_table[:3], 'test_structs/group_basis')

local_basis_for_query = np.vstack(
    (query_basis[0][:3, 3].T, query_basis[0][0, :3]+query_basis[0][:3, 3].T, query_basis[0][1, :3] + query_basis[0][:3, 3].T, query_basis[0][2, :3] + query_basis[0][:3, 3].T))
write_pdb(local_basis_for_query, info_table[:4], 'test_structs/query_basis')

oriented = orient(coordinates, center_sn, axis_sn, radial_sn, query_struct,
                  target_center_sn=query_center_sns[0],
                  target_radial_sn=query_radial_sns[0],
                  target_axis_sn=query_axis_sns[0])

write_pdb(oriented, info_table, 'test_structs/oriented')
print(f'chain center in main: {oriented[center_sn]}')
print(f'the same atom in query struct is {query_basis[0][:3,3].reshape(-1,1)}')
randomized, info_table = random_backbone(oriented, info_table, 1,
                                         'test_residues.txt', chain_list=query_chains[0],
                                         center_sn=center_sn, radius=radius,
                                         radial_sn=radial_sn, radial_v=radial_axis,
                                         axis_vector=symmetry_axis,
                                         symmetry_center=symmetry_center,
                                         multiplicity=multiplicity,
                                         cutoff_distance=0.36)

# make n-1 new identical chains of your oriented chain
# for i in range(multiplicity-1):
#     oriented, info_table, _ = clone_chain(oriented, info_table, [chain])

# target_chains = set(info_table[i]['chain'] for i in range(len(info_table)))

# for i, ch in enumerate(target_chains):

#     same_chain_sns = select(info_table, chain=ch)

#     # translate to put group center point at 0,0,0
#     oriented[same_chain_sns] = oriented[same_chain_sns] - symmetry_center

#     # find the length of the coordinate of the center atom for the chain to be
#     # manipulated
#     curr_rad = oriented[select(info_table, res_num=center_res, atom_name='CA',
#                                chain=ch)].get_length()

#     # find the diference between how much you want the radius to be and its
#     # current value
#     delta_rad = radius - curr_rad

#     # make a vector that moves the current position so that the radius
#     # from the center point will be as desired
#     delta_vect = delta_rad*oriented[select(info_table, res_num=center_res,
#                                            atom_name='CA', chain=ch)].copy().unitize()

#     for sn in same_chain_sns:
#         oriented[sn] = oriented[sn] + delta_vect
#         oriented[sn] = oriented[sn].rotate_arround(
#             2*np.pi/multiplicity*(i+1), vector(group_basis[0, :3]))

#     oriented[same_chain_sns] = oriented[same_chain_sns] + symmetry_center
# write_pdb(oriented, info_table, 'test_structs/test')
# # write_pdb(oriented, info_table, 'test_structs/randomized_pentamer_')

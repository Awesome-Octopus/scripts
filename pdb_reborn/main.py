#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 01:36:49 2023

@author: andrew
"""

# def main():
import numpy as np
import sys
from matplotlib import pyplot as plt
from vector_class import vector
from io_funcs import import_pdb, plot_model, write_pdb
from coord_funcs import cart2cylind, cylind2cart, change_basis, convert_angle
from getter_setter import orient, select, random_backbone, axial_symmetry, get_basis

center_res = 22
axis_res = 10
radial_res = 22
radial_offset_angle = 0
tilt_offset_angle = 0
lean_offset_angle = 0
offset_matrix = None

coordinates, info_table = import_pdb('5x29_S2E_completed.pdb')
center_sn = select(info_table, res_num=center_res, chain='A', atom_name='CA')[0]
# center_sn = 1
axis_sn = select(info_table, res_num=axis_res, chain='A', atom_name='CA')[0]
# axis_sn = 0
radial_sn = select(info_table, res_num=radial_res, chain='A', atom_name='CB')[0]
# radial_sn = 2 
query_struct, query_table = import_pdb('5X29.pdb')


multiplicity = len(set(d['chain'] for d in query_table))
# randomized_coords, info_table = random_backbone(coordinates, info_table, 1, '/home/andrew/scripts/pdb_reborn/residues.txt', 
#                 center_sn=center_sn, radial_sn=radial_sn, axis_vector=axis_vect,
#                 radial_v=radial_vect)axis_sn = select(info_table, res_num=10, chain='A', atom_name='CA')[0]

basis = get_basis(query_struct, center_sn, axis_sn, radial_sn)
# get the orientation of the overall multimer, arround which the monomers have 
# radial symmetry


center_pts = np.zeros((multiplicity,3))
axis_pts = np.zeros((multiplicity,3))
center_sns = select(query_table, res_num=center_res, atom_name='CA')
axis_sns = select(query_table, res_num=axis_res, atom_name='CA')
center_avg = vector([0,0,0])
axis_avg = vector([0,0,0])

# get the mean position of the same atom across the chains of the referece structure
# use the difference in mean position between a center and axis points to define
# an axis of symmetry
for n in range(multiplicity):
    center_avg = center_avg + query_struct[center_sns[n]]
    axis_avg = axis_avg + query_struct[axis_sns[n]]
group_axis = axis_avg/multiplicity - center_avg/multiplicity
group_axis = group_axis.unitize()


target_center_sn = select(query_table, res_num=center_res, chain='A', atom_name='CA')
translation_vector = query_struct[target_center_sn] - coordinates[center_sn]

# align your randomided struture to the quesry structure first
coordinates, rot_mat = orient(coordinates, center_sn, axis_sn, radial_sn, 
                                 target_basis=basis, translate=translation_vector)


write_pdb(coordinates, info_table, 'test_structs/test')
#if we wish to rotate each monomer arround its own long axis paralell to 
# # the group symmetry axis
# if radial_offset_angle != 0:
#     offset_matrix = np.array([[1, 0, 0], [0, np.cos(radial_offset_angle), -np.sin(radial_offset_angle)],
#                                [0, np.sin(radial_offset_angle), np.cos(radial_offset_angle)]])
#     offset_matrix = offset_matrix.T

# # the same logic applies to the other 2 axes
# if tilt_offset_angle != 0:
#     tilt_rot_mat = np.array([[np.cos(tilt_offset_angle), 0, np.sin(tilt_offset_angle)],
#                              [0, 1, 0],
#                              [-np.sin(tilt_offset_angle), 0, np.cos(tilt_offset_angle)]])
#     if offset_matrix is not None:
#         offset_matrix = offset_matrix @ tilt_rot_mat.T
#     else:
#         offset_matrix = tilt_rot_mat.T

# if lean_offset_angle != 0:
#     lean_rot_mat = np.array([[np.cos(lean_offset_angle), -np.sin(lean_offset_angle), 0],
#                              [np.sin(tilt_offset_angle), np.cos(tilt_offset_angle), 0],
#                              [0, 0, 1]])
#     if offset_matrix is not None:
#         offset_matrix = offset_matrix @ lean_rot_mat.T
#     else:
#         offset_matrix = lean_rot_mat.T

# if offset_matrix is not None:
#     basis = basis @ offset_matrix
    
# pentamer_coords, info_table = axial_symmetry(coordinates, info_table, 5, 9, 
#                                              center_sn, axis_sn, radial_sn, 
#                                              translate=False, 
#                                              rot_ax=group_axis,
#                                              target_basis=basis, threshold=0.1)




# if pentamer_coords is not None and info_table is not None:
#     pass
#     # write_pdb(pentamer_coords, info_table, 'test_structs/aligned')
# else:
#     print("the requested radius for multimerization is too small")
    
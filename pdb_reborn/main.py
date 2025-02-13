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
from getter_setter import orient, select, random_backbone, axial_symmetry, get_basis



center_res = 22
axis_res = 10
radial_res = 22
radial_offset_angle = 0
tilt_offset_angle = 0
lean_offset_angle = 0
radius = 1.3
offset_matrix = None
chain = 'A'

coordinates, info_table = import_pdb('5x29_S2E_completed.pdb')
center_sn = select(info_table, res_num=center_res, chain=chain, atom_name='CA')[0]
axis_sn = select(info_table, res_num=axis_res, chain=chain, atom_name='CA')[0]
radial_sn = select(info_table, res_num=radial_res, chain=chain, atom_name='C')[0]

query_struct, query_table = import_pdb('5X29.pdb')
query_center_sn = select(query_table, res_num=center_res, chain=chain, atom_name='CA')[0]
query_axis_sn = select(query_table, res_num=axis_res, chain=chain, atom_name='CA')[0]
query_radial_sn = select(query_table, res_num=radial_res, chain=chain, atom_name='C')[0]
query_basis = get_basis(query_struct, query_center_sn, query_axis_sn, query_radial_sn)


multiplicity = len(set(d['chain'] for d in query_table))


# get the orientation of the overall multimer, arround which the monomers have 
# radial symmetry
# print('the monomer basis inside the multimer was:\n', query_basis)
initial_basis = get_basis(coordinates, center_sn, axis_sn, radial_sn)
# print('while our initial monomer was:\n', initial_basis)


############### THE FOLLOWING CODE IS VERIFIED TO PRODUCE THE CORRECT GROUP AXIS
center_pts = np.zeros((multiplicity,3))
axis_pts = np.zeros((multiplicity,3))
center_sns = select(query_table, res_num=center_res, atom_name='CA')
# print(center_sns)
axis_sns = select(query_table, res_num=axis_res, atom_name='CA')
center_avg = vector([0,0,0])
axis_avg = vector([0,0,0])

# get the mean position of the same atom across the chains of the referece structure
# use the difference in mean position between a center and axis points to define
# an axis of symmetry
for n in range(multiplicity):
    center_avg = center_avg + query_struct[center_sns[n]]
    axis_avg = axis_avg + query_struct[axis_sns[n]]
group_axis = (axis_avg/multiplicity - center_avg/multiplicity).unitize()

#### DEFINE GROUP BASIS########################################################
# we then use this to define a basis for the symmetry group
# first a radial axis, which is the group center to the center atom of chain A
group_rad_ax = (query_struct[center_sns[0]] - 
                center_avg).project_onto_normal_plane(group_axis).unitize()
# then the norm, which is the cross product
group_norm_ax = vector(np.cross(group_axis, group_rad_ax))
group_basis = np.vstack((group_axis, group_rad_ax, group_norm_ax)).T
# print('the principle axis for the multimer group:\n', group_axis)
###############################################################################

translation_vector = query_struct[query_center_sn] - coordinates[center_sn]

#if we wish to rotate each monomer arround its own long axis paralell to 
# the group symmetry axis
#!!! change this to a matrix multiplication operation
if radial_offset_angle !=0:
    for ndx, row in enumerate(coordinates):
        coordinates[ndx]  = vector(coordinates[ndx]).rotate_arround(
            radial_offset_angle, vector(query_basis[0]))

# print('Offset Matrix:\n', offset_matrix)
# if offset_matrix is not None:
#     coordinates = coordinates @ offset_matrix.T

pentamer_coords, info_table = axial_symmetry(coordinates, info_table, 
                                             multiplicity, radius, 
                                             center_sn, axis_sn, radial_sn,
                                             target_basis=query_basis,
                                             cofr=center_avg, 
                                             rot_axis=group_axis,
                                             threshold=0.36)

if pentamer_coords is not None and info_table is not None:

    write_pdb(pentamer_coords, info_table, 'test_structs/aligned')
else:
    print("the requested radius for multimerization is too small")
    
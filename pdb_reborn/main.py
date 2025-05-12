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



center_res = 1
axis_res = 1
radial_res = 1
radial_offset_angle = 0
tilt_offset_angle = 0
lean_offset_angle = 0
radius = 2
offset_matrix = None
chain = 'A'

#
coordinates, info_table = import_pdb('test_structs/test_system2.pdb')
center_sn = select(info_table, res_num=center_res, chain=chain, atom_name='CA')[0]
axis_sn = select(info_table, res_num=axis_res, chain=chain, atom_name='C')[0]
radial_sn = select(info_table, res_num=radial_res, chain=chain, atom_name='N')[0]
test_basis = get_basis(coordinates, center_sn, axis_sn, radial_sn)
# print(test_basis)

query_struct, query_table = import_pdb('test_structs/coordinate_system2.pdb')
query_center_sn = select(query_table, res_num=center_res, chain=chain, atom_name='CA')[0]
query_axis_sn = select(query_table, res_num=axis_res, chain=chain, atom_name='C')[0]
query_radial_sn = select(query_table, res_num=radial_res, chain=chain, atom_name='N')[0]
query_basis = get_basis(query_struct, query_center_sn, query_axis_sn, query_radial_sn)
print(f'query_basis: {query_basis}')

multiplicity = len(set(d['chain'] for d in query_table))
#!!! everything on this line verified to be correct


# ###### THE FOLLOWING CODE IS VERIFIED TO PRODUCE THE CORRECT GROUP AXIS #####
center_pts = np.zeros((multiplicity,3))
# print(f'center_pts: {center_pts}')
axis_pts = np.zeros((multiplicity,3))
center_sns = select(query_table, res_num=center_res, atom_name='CA')
# print(f'center_sns: {center_sns}')
axis_sns = select(query_table, res_num=axis_res, atom_name='C')
# print(f'axis sns {axis_sns}')
center_avg = vector([0,0,0])
axis_avg = vector([0,0,0])

# get the mean position of the same atom across the chains of the reference structure
# use the difference in mean position between a center and axis points to define
# an axis of symmetry
for n in range(multiplicity):
    center_avg = center_avg + query_struct[center_sns[n]]
    axis_avg = axis_avg + query_struct[axis_sns[n]]
    
    
# If you picked center atoms and axis atoms such that the average of the long 
# axis vector for each monomer is 0, it causes a problem. In such case at least
# one atom was defined wrong. warn the user to correct.     
    try:   
        group_axis = (axis_avg/multiplicity - center_avg/multiplicity).unitize()
        # print(f'group axis:\n{group_axis}')
    except ZeroDivisionError:
        print('the atoms you have selected define the long axes of each monomer to',
          'cancel each other out.\n Are you sure you defined your axes correctly?')
# #/### DEFINE GROUP BASIS########################################################
# we then use this to define a basis for the symmetry group
# first a radial axis, which is the group center to the center atom of chain A
if multiplicity > 1:
    group_rad_ax = (query_struct[center_sns[0]] - 
                    center_avg).project_onto_normal_plane(group_axis).unitize()

# if the reference structure is only one chain just use the projection of the 
# vector formed by the radial reference atom coordinate - center reference atom
# coordinate onto the normal plane of the long symmetry axis
else:
    group_rad_ax = (query_struct[query_radial_sn] - \
    center_avg).project_onto_normal_plane(group_axis).unitize()
        
# print(group_rad_ax, group_axis)
# then the norm, which is the cross product
group_norm_ax = vector(np.cross(group_axis, group_rad_ax))
group_basis = np.vstack((group_axis, group_rad_ax, group_norm_ax))

# #############################################################################

translation_vector = center_avg - coordinates[center_sn] + group_rad_ax*radius
coordinates += translation_vector

#if we wish to rotate each monomer arround its own long axis paralell to 
# the group symmetry axis
# !!! change this to a matrix multiplication operation
if radial_offset_angle !=0:
    for ndx, row in enumerate(coordinates):
        coordinates[ndx]  = vector(coordinates[ndx]).rotate_arround(
            radial_offset_angle, vector(query_basis[0]))


if offset_matrix is not None:
    coordinates = coordinates @ offset_matrix.T


multimer_coords, info_table = axial_symmetry(coordinates, info_table, 
                                             multiplicity, radius, 
                                             center_sn, axis_sn, radial_sn,
                                             target_basis=query_basis,
                                             cofr=center_avg, 
                                             rot_axis=group_axis,
                                             threshold=0.36)

if multimer_coords is not None and info_table is not None:
    write_pdb(multimer_coords, info_table, 'test_structs/test')
else:
    print("the requested radius for multimerization is too small")
    
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



center_res =22
axis_res = 22
radial_res = 22
radial_offset_angle = 0
tilt_offset_angle = 0
lean_offset_angle = 0
radius = 2
offset_matrix = None
chain = 'A'
multiplicity = 5

#
coordinates, info_table = import_pdb('5x29_S2E_translated.pdb')
center_sn = select(info_table, res_num=center_res, chain=chain, atom_name='CA')[0]
axis_sn = select(info_table, res_num=axis_res, chain=chain, atom_name='C')[0]
radial_sn = select(info_table, res_num=radial_res, chain=chain, atom_name='N')[0]
test_basis = get_basis(coordinates, center_sn, axis_sn, radial_sn)
# print(f'test_basis:\n{test_basis}')

query_struct, query_table = import_pdb('5X29.pdb')

# store the index atom for each chain in the reference structure in a dictionary
query_center_sn = {}
query_axis_sn = {}
query_radial_sn = {}
query_basis = {}



test = orient(coordinates, center_sn, axis_sn, radial_sn, query_struct)


# get the basis for each chain in the symetric multimer you wish to average over
for ch in set(query_table[i]['chain'] for i in range(len(query_table))):
    query_center_sn = select(query_table, res_num=center_res, chain=ch, atom_name='CA')[0]
    query_axis_sn = select(query_table, res_num=axis_res, chain=ch, atom_name='C')[0]
    query_radial_sn = select(query_table, res_num=radial_res, chain=ch, atom_name='N')[0]
    query_basis[ch] = get_basis(query_struct, query_center_sn, query_axis_sn, query_radial_sn)

# the basis for the whole symmetry group is the average of their bases
# this provides the long axis and the center of geometry, but
# the radial and norm axis are taken from the first chain of the symetry system
# because they usually otherwise end up cancelling each other out.
group_basis = np.zeros((4,4))
for ch in query_basis.keys():  
    group_basis += query_basis[ch]
    
group_basis = group_basis/len(query_basis)
group_basis[1:3, :3] = query_basis[chain][1:3, :3]
    
print(group_basis)
chains_added = []
for i in range(multiplicity-1):
    test, info_table, ch = clone_chain(test, info_table, [chain])
    chains_added.append(ch[0])
    
group_center = group_basis[:3, 3].T
for i, ch in enumerate(chains_added):
    same_chain_sns = select(info_table, chain=chains_added[i])
    for sn in same_chain_sns:
        test[sn] -= group_center
        test[sn] = test[sn].rotate_arround(2*np.pi/multiplicity*(i+1), vector(group_basis[0, :3]))
        test[sn] += group_center
    
write_pdb(test, info_table, 'test_structs/test')


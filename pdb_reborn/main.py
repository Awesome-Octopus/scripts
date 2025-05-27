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
center_res = 34
axis_res = 22
radial_res = 40


chain = 'A'
# CHANGE THIS TO CHANGE FROM A PENTAMER TO ANOTHER NUMBER IF DESIRED
multiplicity = 5


# NOT YET IMPLEMENTED
radial_offset_angle = 0
tilt_offset_angle = 0
lean_offset_angle = 0
offset_matrix = None



coordinates, info_table = import_pdb('test_structs/translated.pdb')
center_sn = select(info_table, res_num=center_res, chain=chain, atom_name='CA')[0]
axis_sn = select(info_table, res_num=axis_res, chain=chain, atom_name='CA')[0]
radial_sn = select(info_table, res_num=radial_res, chain=chain, atom_name='CA')[0]


query_struct, query_table = import_pdb('5X29.pdb')
query_chains = set(query_table[i]['chain'] for i in range(len(query_table)))

query_basis = []

# get the basis for each chain in the symetric multimer you wish to average over
center_pts = []
for ch in query_chains:
    
    query_center_sn = select(query_table, res_num=center_res, chain=ch, atom_name='CA')[0]
    center_pts.append(query_struct[query_center_sn])
    query_axis_sn = select(query_table, res_num=axis_res, chain=ch, atom_name='CA')[0]
    query_radial_sn = select(query_table, res_num=radial_res, chain=ch, atom_name='CA')[0]
    query_basis.append(get_basis(query_struct, query_center_sn, query_axis_sn, query_radial_sn))

# the average of all the atoms defining the center points of each chain    
center_pt = vector(np.mean(np.asarray(center_pts), axis=0))

# the group basis is the average of all the bases of the query chains
a = np.asarray(query_basis)
# unitize the inner matrix
group_basis = np.mean(a, axis=0)
for i in range(3):
    group_basis[i, :3] = group_basis[i, :3].view(vector).unitize()
    
oriented = orient(coordinates, center_sn, axis_sn, radial_sn, query_struct,
              target_center_sn=query_center_sn,
              target_radial_sn=query_radial_sn, target_axis_sn=query_axis_sn)
# write_pdb(test, info_table, 'test_structs/oriented')

# find the average distance of the center coordinate of each chain from the 
# center of the group
query_radius = 0
for sn in select(query_table, res_num=center_res, atom_name='CA'):
    query_radius += center_pt.distance_between(query_struct[sn])
query_radius = query_radius/len(query_chains)

# how wide you want the radius of the circular arrangement of monomers to be
radius = query_radius

# make n-1 new identical chains of your oriented chain
for i in range(multiplicity-1):
    oriented, info_table, _ = clone_chain(oriented, info_table, [chain])
    
target_chains = set(info_table[i]['chain'] for i in range(len(info_table)))

for i, ch in enumerate(target_chains):
    
    same_chain_sns = select(info_table, chain=ch)
    
    # translate to put group center point at 0,0,0
    oriented[same_chain_sns] = oriented[same_chain_sns] - center_pt
    
    # find the length of the coordinate of the center atom for the chain to be 
    # manipulated
    curr_rad = oriented[select(info_table, res_num=center_res, atom_name='CA', 
                               chain=ch)].get_length()
    
    # find the diference between how much you want the radius to be and its 
    # current value
    delta_rad = radius - curr_rad
    
    # make a vector that moves the current position so that the radius 
    # from the center point will be as desired
    delta_vect = delta_rad*oriented[select(info_table, res_num=center_res,
                                       atom_name='CA', chain=ch)].copy().unitize()
    
    for sn in same_chain_sns:
        oriented[sn] = oriented[sn] + delta_vect
        oriented[sn] = oriented[sn].rotate_arround(
            2*np.pi/multiplicity*(i+1), vector(group_basis[0, :3]))
    
    oriented[same_chain_sns] = oriented[same_chain_sns] + center_pt

write_pdb(oriented, info_table, 'test_structs/test')
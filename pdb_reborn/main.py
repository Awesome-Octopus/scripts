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
from getter_setter import orient, select, get_transform_mat, random_backbone

coordinates, info_table = import_pdb('test_structs/coordinate_system.pdb')
center_sn = select(info_table, res_num=1, chain='A', atom_name='CA')[0]

axis_sn = select(info_table, res_num=1, chain='A', atom_name='N')[0]
radial_sn = select(info_table, res_num=1, chain='A', atom_name='C')[0]
query_struct, info_table = import_pdb('test_structs/coordinate_system2.pdb')
multiplicity = len(set(d['chain'] for d in info_table))

# trans_mat = get_transform_mat(
#                     coordinates, query_struct, center_sn, axis_sn, radial_sn, 
#                     target_center_sn=center_sn, target_axis_sn=axis_sn, 
#                     target_radial_sn=radial_sn)

new_coords = orient(coordinates, query_struct, center_sn, axis_sn, radial_sn, 
                    target_center_sn=center_sn, target_axis_sn=axis_sn, 
                    target_radial_sn=radial_sn)

write_pdb(new_coords, info_table, 'test_structs/test_alignment')
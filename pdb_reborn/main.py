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

coordinates, info_table = import_pdb('5x29_S2E_completed.pdb')
center_sn = select(info_table, res_num=22, chain='A', atom_name='CA')[0]

axis_sn = select(info_table, res_num=10, chain='A', atom_name='CA')[0]
radial_sn = select(info_table, res_num=22, chain='A', atom_name='CB')[0]
query_struct, _ = import_pdb('5X29.pdb')
multiplicity = len(set(d['chain'] for d in info_table))

# randomized_coords, info_table = random_backbone(coordinates, info_table, 1, '/home/andrew/scripts/pdb_reborn/residues.txt', 
#                 center_sn=center_sn, radial_sn=radial_sn, axis_vector=axis_vect,
#                 radial_v=radial_vect)
basis = get_basis(query_struct, center_sn, axis_sn, radial_sn)
print(basis)
# new_coords, rot_mat = orient(coordinates, center_sn, axis_sn, radial_sn, 
#                     target_basis=np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]]), 
#                     translate=False)

pentamer_coords, info_table = axial_symmetry(coordinates, info_table, 5, 50, 
                                             center_sn, axis_sn, radial_sn,  
                                             translate=vector([1, 1, 1]))

write_pdb(pentamer_coords, info_table, 'test_structs/test_output')
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
from getter_setter import orient, select, random_backbone, axial_symmetry

coordinates, info_table = import_pdb('test_structs/poly_gly_aligned_helix.pdb')
center_sn = select(info_table, res_num=1, chain='A', atom_name='CA')[0]

axis_sn = select(info_table, res_num=1, chain='A', atom_name='N')[0]
radial_sn = select(info_table, res_num=1, chain='A', atom_name='C')[0]
query_struct, info_table = import_pdb('test_structs/poly_gly_tilted_helix.pdb')
multiplicity = len(set(d['chain'] for d in info_table))

axis_vect = coordinates[axis_sn] - coordinates[center_sn]
radial_vect = coordinates[radial_sn] - coordinates[center_sn]
randomized_coords, info_table = random_backbone(coordinates, info_table, 1, '/home/andrew/scripts/pdb_reborn/residues.txt', 
                center_sn=center_sn, radial_sn=radial_sn, axis_vector=axis_vect,
                radial_v=radial_vect)


new_coords, _ = orient(randomized_coords, center_sn, axis_sn, radial_sn, 
                    target_center_sn=center_sn, target_axis_sn=axis_sn, 
                    target_radial_sn=radial_sn, target_coords=query_struct)

# pentamer_coords, info_table = axial_symmetry(new_coords, info_table, 5, 5, 
#                                              center_sn, axis_sn, radial_sn)

write_pdb(pentamer_coords, info_table, 'test_structs/test_output')
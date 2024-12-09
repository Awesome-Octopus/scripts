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
from getter_setter import orient, select, get_struct_orientation, random_backbone

coordinates, info_table = import_pdb('test_structs/poly_gly_aligned_helix.pdb')
center_sn = select(info_table, res_num=1, chain='A', atom_name='CA')[0]

axis_sn = select(info_table, res_num=5, chain='A', atom_name='N')[0]
radial_sn = select(info_table, res_num=1, chain='A', atom_name='C')[0]
multiplicity = len(set(d['chain'] for d in info_table))

radius, radial_angle, cross_prod, axial_offset_ang = get_struct_orientation('test_structs/poly_gly_aligned_helix.pdb', center_sn, axis_sn, radial_sn)
axis = vector([0,0,1]).rotate_arround(-axial_offset_ang, cross_prod)
print(axis)
# coordinates, info_table = import_pdb('test_structs/poly_gly_tilted_helix.pdb')

coordinates, points_to_center = orient(coordinates, info_table, center_sn, axis_sn, radial_sn, radial_angle=radial_angle, reference_vect=axis)

# write_pdb(coordinates, info_table, 'test_structs/test')

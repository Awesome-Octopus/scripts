#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 22:16:59 2024

@author: andrew
"""
res_num = 48
prev_ca = select(info_table, res_num=res_num, atom_name='CA')
prev_c = select(info_table, res_num=res_num, atom_name='C')
n = select(info_table, res_num=res_num+1, atom_name='N')
ca = select(info_table, res_num=res_num+1, atom_name='CA')
current = measure_dihedral(
    [coordinates[prev_ca], coordinates[prev_c], coordinates[n], coordinates[ca]])
rotation = np.pi/2 - current
omega = coordinates[n]-coordinates[prev_c]
for i in range(n, len(coordinates)):
    coordinates[i] -= coordinates[prev_c]
    coordinates[i] = coordinates[i].rotate_arround(rotation, omega)
    coordinates[i] += coordinates[prev_c]
print(measure_dihedral([coordinates[prev_ca],
      coordinates[prev_c], coordinates[n], coordinates[ca]]))

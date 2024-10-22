#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 01:36:49 2023

@author: andrew
"""

# def main():
import numpy as np
import sys
from vector_class import vector
from io_funcs import import_pdb, plot_model, write_pdb
from coord_funcs import cart2cylind, cylind2cart, change_basis, convert_angle
from getter_setter import orient, select, get_struct_orientation, random_backbone


#from getter_setter import select, get_struct_orientation, orient, \
#    random_backbone


coordinates, info_table = import_pdb('test_structs/orthogonal.pdb')
center_sn = select(info_table, res_num=1, chain='A', atom_name='CA')[0]

axis_sn = select(info_table, res_num=1, chain='A', atom_name='N')[0]
radial_sn = select(info_table, res_num=1, chain='A', atom_name='C')[0]
print(f' this should be [2, 0, 1]: {coordinates[axis_sn]}')
print(f' this should be [2, 1, 0]: {coordinates[radial_sn]}')
print(f' this should be [2, 0, 0]: {coordinates[center_sn]}')
multiplicity = len(set(d['chain'] for d in info_table))


radius, radial_angle, cross_prod, offset_ang = get_struct_orientation('test_structs/orthogonal.pdb', center_sn, axis_sn, radial_sn)
print(f' axial offset angle: {offset_ang} \n cross product vector: {cross_prod}')
print(f'for ref struct: radial angle {radial_angle}')
# rotate the principle axis so that it is parallel with +z
if not np.allclose(cross_prod, vector([0, 0, 0])):

    axis_vector = vector([0, 0, 1]).rotate_arround(offset_ang, cross_prod)
else:
    cross_prod = vector([0, 1, 0])
    axis_vector = vector([0, 0, 1])

print(f'axis vector {axis_vector}')


for w in range(1):

    print('--------- query structure --------------\n')
    coordinates, info_table = import_pdb(
        'test_structs/internal_90deg_rot_about_z.pdb')
    center_sn = select(info_table, res_num=1,
                        chain='A', atom_name='CA')[0]
    axis_sn = select(info_table, res_num=1, chain='A', atom_name='N')[0]
    radial_sn = select(info_table, res_num=1, chain='A', atom_name='C')[0]


    # orient the structure to align reference with z
    coordinates, points_to_center = orient(coordinates, info_table, center_sn,
                                             axis_sn, radial_sn,
                                             radial_angle=radial_angle,
                                             reference_vect=axis_vector,
                                             recenter=False)
    print(f'points_to_center:                    {points_to_center}')
    write_pdb(coordinates, info_table, 'test_structs/test')
    if points_to_center is not None:



        # points_to_center = points_to_center.rotate_arround(offset_ang, cross_prod)

        radial_v = points_to_center.copy()*radius

        print(f'radial_v {radial_v}')
        # multimeric_center = coordinates[center_sn] - radial_v
        # print(f'multimeric_center:                    {multimeric_center}')

#         set multimeric center to origin
        # coordinates -= multimeric_center
    # else:
        # print('a good radial coordinate is needed to adjust the internal radial'
              # 'angles of each subunit when multimerizing')

    # print(f'center coordinate: {coordinates[center_sn]}')
#     print(radial_sn)
#     print(f'radius: {radius}')
#     # bundle all of these together to be passed along
#     rot_sym_vars = {'multiplicity' : multiplicity, 'radial_v' : radial_v, \
#                    'axis_vector' : axis_vector,
#                    'radius' : radius,
#                    'center_sn' : center_sn,
#                    'cutoff_distance' : 0.01,
#                    'radial_sn' : radial_sn,
#                    'multimeric_center' : multimeric_center}
#     print(f'\naxis vector {axis_vector}')
#     print(f'\nradial_vector {radial_v}')
    # print(points_to_center.angle_between(multimeric_center))
    # print(axis_vector.is_orthogonal(radial_v))


#     # the inclusion of rot_sym_vars allows the machine to check symmetry as it is randomizing
# #     coordinates, info_table = random_backbone(coordinates, info_table, 1,
# #                                               'residues.txt', max_tries=100,
# #     check_method='definitive', sym_clash_check=True, **rot_sym_vars)

# # #     output1 = axial_symmetry(coordinates, info_table,
# # #                               multiplicity, radius,
# # #                               center_sn,
# # #                               axis_sn, radial_sn,
# # #                               radial_angle=radial_angle, threshold=7)
# # #     if output1 is not None:
# # #         coordinates, info_table = output1
# # #         plot_model(coordinates, title=f'2 pi / {w} radial angle')
# # #     #     write_pdb(coordinates, info_table,
# # #     #               f'randomized_pentamers/randomized_pentamer_{w}')

# # #     # else:
# # #     #     print(
# # #     #         f'you cannot have {multiplicity}-fold radial symmetry for model {w}')
        # write_pdb(coordinates, info_table, 'test_structs/reoriented')


# # return coordinates, info_table, center_sn


# # if __name__ == 'main':
# #     main()

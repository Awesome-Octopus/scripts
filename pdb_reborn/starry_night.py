#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 14 21:48:54 2025

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

model_group = []
star, table = import_pdb('test_structs/star.pdb')
model_group.append([star.copy(), table.copy()])
for i in range(1, 6):
    star, table, new_chain
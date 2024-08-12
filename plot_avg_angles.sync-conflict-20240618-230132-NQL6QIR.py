#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 23:06:13 2024

@author: andrew
"""

import numpy as np
import MDAnalysis as mda
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

parser.add_argument('-f', '--file', required=True, type=str)
args = parser.parse_args()


data = np.loadtxt(args.file)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.errorbar(data[:,2], data[:,0], xerr = data[:,1], yerr= data[:,3], fmt= 'o')
ax.set_xlabel('yaw (degrees)')
ax.set_ylabel('pitch (degrees)')
plt.show()
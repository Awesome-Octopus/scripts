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

# sometimes the roll angle couldn't be calculated because the helix axis is
# essentially coolinear with the symmetry axis, in such a case,
# nan is reported. we will tell it to treat it as 0
data[:,4] = np.where(np.isnan(data[:,4]), 0, data[:,4])
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
sca = ax.scatter(data[:,2], data[:,0], c=data[:,4], cmap='viridis')
ax.errorbar(data[:,2], data[:,0], xerr = data[:,3], yerr= data[:,1],
            fmt='none', ecolor='black', alpha=0.2, capsize=2)
ax.set_xlabel('yaw (degrees)')
ax.set_ylabel('pitch (degrees)')
plt.colorbar(sca, label='mean roll angle (deg.)')
plt.show()

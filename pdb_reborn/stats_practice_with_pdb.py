#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 18:41:59 2025

@author: andrew
"""

import numpy as np
import seaborn as sns
import scipy.stats as stats
import pandas as pd
import sys
import sklearn as skl
import matplotlib as mpl
from matplotlib import pyplot as plt
from vector_class import vector
from io_funcs import import_pdb, plot_model, write_pdb
from coord_funcs import cart2cylind, cylind2cart, change_basis, convert_angle
from getter_setter import orient, select, random_backbone, axial_symmetry, get_basis

coords, info = import_pdb('/home/andrew/scripts/pdb_reborn/5X29.pdb')
df = pd.DataFrame({'x': coords[:,0], 'y': coords[:,1], 'z': coords[:,2]})

# n =stats.norm.pdf(np.linspace(xmin, xmax, 100), loc=df['x'].mean(), scale=df['x'].std())
hing  = skl.datasets.fetch_california_housing()
hing = pd.DataFrame(data=hing['data'], columns=hing['feature_names'])

mock = stats.norm.pdf(x=np.linspace(np.min(hing['MedInc']), np.max(hing['MedInc']), len(hing['MedInc'])), loc=hing['MedInc'].median(), scale=hing['MedInc'].std())
sns.histplot(hing['MedInc'])
plt.plot(mock*1000)

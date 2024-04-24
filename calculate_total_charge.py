#!/usr/bin/env python3

import sys, re

# f = open(sys.argv[1], 'r')
f = open('residues.txt', 'r')

# we want anything that matches a list of positively charged lipids
a = open('anions.txt', 'r')
lines = a.readlines()
anion_dict = {}
for anion in lines:
    anion = line.strip()


lines = f.readlines()
lipid_charges = {''}
for line in lines:
    line = line.strip()


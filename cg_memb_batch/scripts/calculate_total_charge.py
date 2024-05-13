#!/usr/bin/env python3

# make sure the names of ions in your ion dictionary match names in topology file
# exactly!!!


import sys
import re
import os

if not len(sys.argv) == 4:
    print("Usage: python calculate_total_charge.py /path/to/ion_list path/to/system_top protein_net_charge")
    sys.exit(1)

if not all((type(sys.argv[1]) is str, type(sys.argv[2]) is str)):
    print("Usage: python calculate_total_charge.py /path/to/ion_list path/to/system_top protein_net_charge")
    sys.exit(1)


# we want a dictionary of ion names and their charges
a = open(sys.argv[1], 'r')
lines = a.readlines()
ion_dict = {}
for ion in lines:
    ion = ion.strip().split()
    if ion != '' and len(ion) == 2:
        ion_dict.update({ion[0]: int(ion[1])})

    else:
        print("ion list should be a text file with one ion per line followed by the charge. e.g. : Na 1")
        sys.exit(1)

top_file = open(sys.argv[2], 'r')

lines = top_file.readlines()
valid = set("-+0123456789")
if not all(char in valid for char in sys.argv[3]):
    print("protein net charge must be an integer")
    sys.exit(1)
else:
    total_charge = int(sys.argv[3])

past_header = False
for line in lines:
    line = line.strip()
    # ignore everything until we hit a line that matches the header for molecules
    if line == "[ molecules ]":
        past_header = True

    # if the line doesnt start with a semicolon and isnt empty

    if past_header and line != [] and line[0] != ";" and line[0] != "[":
        line = line.split()
        if len(line) != 2:
            print(
                "Each molecule designation in the top file must be followed by a number, e.g. POPE 1")
            exit(1)
        else:
            molecule = line[0]
        if molecule in ion_dict.keys():
            total_charge += (int(ion_dict[molecule]) * int(line[1]))

print(total_charge)

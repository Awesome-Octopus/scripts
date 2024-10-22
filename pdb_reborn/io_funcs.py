#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 00:47:16 2024

@author: andrew
"""

# %% I/O functions
import numpy as np
from vector_class import vector
from Bio.PDB.PDBParser import PDBParser
# how we will import the raw text
# from Bio.PDB import PDBIO

import matplotlib.pyplot as plt

def import_pdb(fname=None):
    # we are making an numpy array filled with the x y and z for each atom in each
    # row doing it as a fixed array because this is much more memory efficient

    if fname is None:
        fname = 'protein.pdb'
    outfile_name = fname[0:-4]
    structure_id = "mode_7"
    parser = PDBParser(PERMISSIVE=1)
    # permissive=1 allows it to accept incomplete structures without error
    whole_pdb = parser.get_structure(outfile_name, fname)
    print(f'Coordinates for {fname} imported.')
    n = 0
    info_table = []
    coordinates = []
    for model in whole_pdb:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    a_tuple = atom.full_id
                    a_dict = {
                        'ser_num': n, 'model_name': a_tuple[0],
                        'submodel': a_tuple[1], 'chain': a_tuple[2],
                        'res_num': a_tuple[3][1], 'res_name':
                        residue.get_resname(), 'atom_name': a_tuple[4][0],
                        'phi': float('nan'), 'psi': float('nan'),
                        'omega': float('nan')
                    }
                    # the last 2 are empty fields where phi & psi are assigned

                    info_table.append(a_dict)
                    n = n + 1
    coordinates = np.zeros((3, 3))
    coordinates = np.zeros(shape=(n, 3))
    coordinates = coordinates.view(vector)
    n = 0
    for model in whole_pdb:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coordinates[n] = atom.get_coord()
                    n = n + 1
    return coordinates, info_table


def plot_model(coordinates, fig=None, title='structure', plots_per_row=3):

    i = coordinates[:, 0]
    j = coordinates[:, 1]
    k = coordinates[:, 2]

    if fig is None:
        fig = plt.figure()

    num_subplots = len(fig.axes)

    # Calculate subplot dimensions
    rows = (num_subplots // plots_per_row) + 1
    cols = min(num_subplots % plots_per_row + 1, plots_per_row)

    # Add subplot
    ax = fig.add_subplot(rows, cols, num_subplots + 1, projection='3d')
    ax.set_box_aspect([1, 1, 1])

    v = ax.scatter3D(i, j, k, c=np.sqrt(i**2 + k**2))
    # ax.set_xlim(-40, 40)
    # ax.set_ylim(-40, 40)
    # ax.set_zlim(-40, 40)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.set_title(title)
    fig.tight_layout()

    return fig


def import_from_vmd(filename=None):

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #  UNDER CONSTRUCTION
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if filename is None:
        pass
    else:
        try:
            w = 1

            filename = f'cg-trajectories/S2E-processed-cg-input-{w}.pdb'
            # outfile_name = f'randomized_S2E/randomized-S2E-oriented-{w}.pdb'
            whole_pdb = PDBParser.get_structure(w, filename)
            n = 0
            global info_table
            info_table = []
            for model in whole_pdb:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            a_tuple = atom.full_id
                            a_dict = {
                                'ser_num': n, 'model_name': a_tuple[0],
                                'submodel': a_tuple[1], 'chain': a_tuple[2],
                                'res_num': a_tuple[3][1], 'res_name':
                                residue.get_resname(),
                                'atom_name': a_tuple[4][0],
                                'phi': float('nan'), 'psi': float('nan')
                            }
                            # the last 2 are empty fields where phi &
                            # psi are assigned

                            info_table.append(a_dict)
                            n = n + 1
            coordinates = np.zeros(shape=(n, 3))
            coordinates = coordinates.view(vector)

        except:
            pass


def write_pdb(coordinates, info_table, outfile=None):
    """
    coordinates : Array of float 64
        the atomic coordinates of all atoms
    info_table : list
        list of dictionaries containing data for each coordinate
    Output a pdb text file of all currently loaded models with a default
    filename if none is given

    Returns
    -------
    None.

    """

    if outfile is None:
        model_name = info_table[0]['model_name']
    else:
        model_name = outfile

    # if the filename is already taken, augment the number at the end of the
    # filename

    try:
        f = open(f'{model_name}.pdb', 'x')
    except FileExistsError:
        count = 0
        while True:
            try:
                f = open(f'{model_name}_{count}.pdb', 'x')
            except FileExistsError:
                count = count + 1
            else:
                break

    f.write(
        'MODEL        1                                                 '
        '               \n')
    num_termini = 0
    new_ser_num = 0
    model_count = 1
    # we need to keep track of how many times we hit the end of a chain because
    # the TER line gets its own serial number and we have to increase all
    # subsequent serial numbers by 1 for each TER we have encountered

    for n, a in enumerate(info_table):
        new_ser_num = new_ser_num + 1
        element = a['atom_name'][0]

        # the element is always the first letter in the name eg HG1 for gamma
        # hydrogen

        if element not in ['C', 'N', 'O', 'S', 'H']:

            # consider anything that is not one of these a heteroatom, note
            # this doesn't include phosphorous
            tag = 'HETATM'

        else:
            tag = 'ATOM  '

        if n != 0:
            if info_table[n-1]['chain'] != a['chain']:

                # if the chain for this is different from the last

                num_termini = num_termini + 1
                f.write(
                    f'TER   {n+num_termini:>5d}      '
                    f'{info_table[n-1]["res_name"]:<} '
                    f'{info_table[n-1]["chain"]}'
                    f'{info_table[n-1]["res_num"]: > 4d}'
                    f'{chr(10)}')
                new_ser_num = new_ser_num + 1
            if info_table[n-1]['submodel'] != a['submodel']:
                model_count = model_count + 1
                f.write(
                    'ENDMDL                                               '
                    '                           \n')
                f.write(
                    f'MODEL     {model_count:>4d}                           '
                    f'                                       {chr(10)}')
                num_termini = 0
                new_ser_num = 1
                # reset the numbering when you encounter a new model
        f.write(
            f'{tag}{new_ser_num:>5d} {a["atom_name"]:<4s} {a["res_name"]:<4s}'
            f'{a["chain"]}{a["res_num"]:>4d}{coordinates[n][0]:>12.3f}'
            f'{coordinates[n][1]:>8.3f}{coordinates[n][2]:>8.3f}  1.00  '
            f'0.00          {element:>2s}{chr(10)}')

    # after the very last atom print the final footer

    f.write(
        f'TER   {new_ser_num + 1:>5d}      {info_table[n-1]["res_name"]:<} '
        f'{info_table[n-1]["chain"]}{info_table[n-1]["res_num"]:>4d}         '
        f'                                               {chr(10)}')

    # this follows pdb guidelines for fixed width columns

    f.write(
        'ENDMDL                                                             '
        '             \n')
    f.write(
        'END                                                                '
        '             ')

    f.close()

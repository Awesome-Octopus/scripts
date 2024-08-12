# -*- coding: utf-8 -*-

import re, numpy

# for each seq in seqs an itp-topology will be generated:
seqs = ['MSSMSSMSSSMSSMSSSSSSSSMSSSSSSMSSSSSSSSSSMSMSSSSSMSSM', 'MSSMSMSSMSSMSSSMSSMSMSMSMSMSMSSMSSSSMSSSSSSSSSSMSSMS', 'SSSMSMSSSSSMSSMSSSMSSSSMSMSSMSSMSSSSSMSSSMSSSMSSSMSM', 'MSMSSMSSSMSSSSSMSSSMSSSSSSMSSMSSSSSSSMSSSSMSSSSMSSSM']

class Atom:
    '''Atom object'''
    def __init__(self, id, type, resid, resname, name, cgnr, charge):
        '''Init function for Atom object'''
        self.id = int(id)
        self.resid = int(resid)
        self.resname = str(resname)
        self.type = str(type)
        self.name = str(name)
        self.charge = float(charge)
        self.cgnr = int(cgnr)
    
    def fmt(self):
        '''formatted line for Atom object'''
        return "%3d %5s %3d %5s %5s %3d %1.3f" \
            % (self.id, self.type, self.resid, self.resname, self.name, self.cgnr, self.charge)


class Bond:
    '''Bond object'''
    def __init__(self, atom1, atom2, l, k, id):
        '''Init function for Bond object'''
        self.id = id # unique internal id of the bond
        self.atom1 = int(atom1) # atom1 unique id as in the itp
        self.atom2 = int(atom2) # atom2 unique id as in the itp
        self.type = 1 # type of bond
        self.k = float(k) # current k
        self.l = float(l) # current l

    def fmt(self):
        '''formatted line for Bond object'''
        return "%3d %3d %1d %5.3f %5.3f" \
            % (self.atom1, self.atom2, self.type, self.l, self.k)

class Constraint:
    '''Constraint object'''
    def __init__(self, atom1, atom2, l, id):
        '''Init function for Bond object'''
        self.id = id # unique internal id of the bond
        self.atom1 = int(atom1) # atom1 unique id as in the itp
        self.atom2 = int(atom2) # atom2 unique id as in the itp
        self.type = 1 # type of bond
        self.l = float(l) # current l

    def fmt(self):
        '''formatted line for Constraint object'''
        return "%3d %3d %1d %5.3f" \
            % (self.atom1, self.atom2, self.type, self.l)


class Angle:
    '''Angle object'''
    def __init__(self, atom1, atom2, atom3, phi, k, id):
        '''Init function for Angle object'''
        self.id = id
        self.atom1 = int(atom1)
        self.atom2 = int(atom2)
        self.atom3 = int(atom3)
        self.type = 2
        self.k = float(k)
        self.phi = float(phi)
  
    def fmt(self):
        '''formatted line for Angle object'''
        return "%3d %3d %3d %1d %5.3f %5.3f" \
            % (self.atom1, self.atom2, self.atom3, self.type, self.phi, self.k)

class Dihedral:
    '''Dihedral object'''
    def __init__(self, atom1, atom2, atom3, atom4, type, phi, k, id, mult):
        '''Init function for Angle object'''
        self.id = id
        self.atom1 = int(atom1)
        self.atom2 = int(atom2)
        self.atom3 = int(atom3)
        self.atom4 = int(atom4)
        self.type = int(type)
        self.k = float(k)
        self.phi = float(phi)
        self.mult = int(mult)
  
    def fmt(self):
        '''formatted line for Dihedral object'''
        return "%3d %3d %3d %3d %1d %5.3f %5.3f %1d" \
            % (self.atom1, self.atom2, self.atom3, self.atom4, self.type, self.phi, self.k, self.mult)

class Molecule:
    def __init__(self):
        '''Init function for Molecule object'''
        self.molname = ''
        self._atoms = {}
        self._bonds = []
        self._constraints = []
        self._angles = []
        self._dihedrals = []
        self._exclusions = []
        self.atom_num = 0
        self.res_num = 0
        
    def write_itp(self, fname):
        '''Write an itp file'''
        out = ['; AUTOMATICALLY GENERATED TOPOLOGY FILE AT ITERATION']
        out.append(' [ moleculetype ] ')
        out.append('; molname         nrexcl')
        out.append(str(self.molname) + '   1')
        out.append(' [ atoms ] ')
        out.append(';id type resnr residu atom cgnr charge')
        out.extend([i.fmt() for i in self._atoms.values()])
        out.append(' [ bonds ] ')
        out.append(';atom_id1 atom_id2 func_type length force_c')
        out.extend([i.fmt() for i in self._bonds])
        out.append(' [ constraints ] ')
        out.append(';atom_id1 atom_id2 func_type length')
        out.extend([i.fmt() for i in self._constraints])
        out.append(' [ angles ] ')
        out.append(';atom_id1 atom_id2 atom_id3 func_type angle force_c')
        out.extend([i.fmt() for i in self._angles])
        out.append(' [ dihedrals ] ')
        out.append(';atom_id1 atom_id2 atom_id3 atom_id4 func_type dih force_c')
        out.extend([i.fmt() for i in self._dihedrals])
        out.append(' [ exclusions ] ')
        out.extend(self._exclusions)
        fout = open(fname, 'w')
        fout.write("\n".join(out))
        fout.close()

for seq_id in range(len(seqs)):
    chain = seqs[seq_id]
    
    m = Molecule()
    m.molname = 'SM' + str(seq_id + 1)
    
    par_b_out = []
    par_an_out = []
    par_dih_out = []
    
    for i,j in enumerate(chain):
        if j == 'M':
            m.res_num += 1
            m.atom_num += 1
            m._atoms[int(m.atom_num)] = Atom(int(m.atom_num), 'SC3', m.res_num, 'MAL', 'BB', int(m.atom_num), 0)
            m.atom_num += 1
            m._atoms[int(m.atom_num)] = Atom(int(m.atom_num), 'P1', m.res_num, 'MAL', 'SC1', int(m.atom_num), 0)
            m.atom_num += 1
            m._atoms[int(m.atom_num)] = Atom(int(m.atom_num), 'Qa', m.res_num, 'MAL', 'SC2', int(m.atom_num), -1)
            
            m._bonds.append(Bond(m.atom_num - 2, m.atom_num - 1, 0.214, 1079, 1))
            par_b_out.append(str(m.atom_num - 2) + ' ' + str(m.atom_num - 1) + ' 1')
            m._bonds.append(Bond(m.atom_num - 2, m.atom_num    , 0.214, 1079, 1))
            par_b_out.append(str(m.atom_num - 2) + ' ' + str(m.atom_num) + ' 1')
            m._bonds.append(Bond(m.atom_num - 1, m.atom_num    , 0.374, 6030, 2))
            par_b_out.append(str(m.atom_num - 1) + ' ' + str(m.atom_num) + ' 2')
            
            # bond to the next unit
            if i < (len(chain) - 1):
                m._bonds.append(Bond(m.atom_num - 2, m.atom_num + 1, 0.210, 2696, 3))
                par_b_out.append(str(m.atom_num - 2) + ' ' + str(m.atom_num + 1) + ' 3')
            
            # angle to the previous unit
            if i > 0:
                 if chain[i - 1] == 'S':
                    m._angles.append(Angle(m.atom_num - 1, m.atom_num - 2, m.atom_num - 6, 94.5, 121, 2))
                    par_an_out.append(str(m.atom_num - 1) + ' ' + str(m.atom_num - 2) + ' ' + str(m.atom_num - 6) + ' 1')
                    m._angles.append(Angle(m.atom_num    , m.atom_num - 2, m.atom_num - 6, 165, 121, 2))
                    par_an_out.append(str(m.atom_num) + ' ' + str(m.atom_num - 2) + ' ' + str(m.atom_num - 6) + ' 2')
            
            # angles to the next unit
            if i < (len(chain) - 1):
                if chain[i + 1] == 'S':
                    m._angles.append(Angle(m.atom_num - 1, m.atom_num - 2, m.atom_num + 1, 133, 118, 4))
                    par_an_out.append(str(m.atom_num - 1) + ' ' + str(m.atom_num - 2) + ' ' + str(m.atom_num + 1) + ' 3')
                    m._angles.append(Angle(m.atom_num    , m.atom_num - 2, m.atom_num + 1, 93, 118, 4))
                    par_an_out.append(str(m.atom_num) + ' ' + str(m.atom_num - 2) + ' ' + str(m.atom_num + 1) + ' 4')
    
        if j == 'S':
            m.res_num += 1
            m.atom_num += 1
            m._atoms[int(m.atom_num)] = Atom(int(m.atom_num), 'SC1', m.res_num, 'STY', 'BB', int(m.atom_num), 0)
            m.atom_num += 1
            m._atoms[int(m.atom_num)] = Atom(int(m.atom_num), 'SC5', m.res_num, 'STY', 'SC1', int(m.atom_num), 0)
            m.atom_num += 1
            m._atoms[int(m.atom_num)] = Atom(int(m.atom_num), 'SC5', m.res_num, 'STY', 'SC2', int(m.atom_num), 0)
            m.atom_num += 1
            m._atoms[int(m.atom_num)] = Atom(int(m.atom_num), 'SC5', m.res_num, 'STY', 'SC3', int(m.atom_num), 0)
            
            m._bonds.append(Bond(m.atom_num - 3, m.atom_num - 2, 0.202, 7512, 4))
            par_b_out.append(str(m.atom_num - 3) + ' ' + str(m.atom_num - 2) + ' 4')
            m._constraints.append(Constraint(m.atom_num - 2, m.atom_num - 1, 0.27, 5))
            m._constraints.append(Constraint(m.atom_num - 2, m.atom_num    , 0.27, 5))
            m._constraints.append(Constraint(m.atom_num - 1, m.atom_num    , 0.27, 5))
            
            # bond to the next unit
            if i < (len(chain) - 1):
                if chain[i + 1] == 'M':
                    m._bonds.append(Bond(m.atom_num - 3, m.atom_num + 1, 0.210, 2696, 3))
                    par_b_out.append(str(m.atom_num - 3) + ' ' + str(m.atom_num + 1) + ' 3')
                if chain[i + 1] == 'S':
                    m._bonds.append(Bond(m.atom_num - 3, m.atom_num + 1, 0.191, 10000, 3))
                    par_b_out.append(str(m.atom_num - 3) + ' ' + str(m.atom_num + 1) + ' 5')
            
            # angles in the current unit
            m._angles.append(Angle(m.atom_num - 3, m.atom_num - 2, m.atom_num - 1, 150, 52, 5))
            m._angles.append(Angle(m.atom_num - 3, m.atom_num - 2, m.atom_num    , 150, 52, 5))
            
            # angle to the previous unit
            if i > 0:
                if chain[i - 1] == 'M':
                    m._angles.append(Angle(m.atom_num - 2, m.atom_num - 3, m.atom_num - 6, 107, 121, 6))
                    par_an_out.append(str(m.atom_num - 2) + ' ' + str(m.atom_num - 3) + ' ' + str(m.atom_num - 6) + ' 5')
                if chain[i - 1] == 'S':
                    m._angles.append(Angle(m.atom_num - 2, m.atom_num - 3, m.atom_num - 7, 125, 53, 7))
                    par_an_out.append(str(m.atom_num - 2) + ' ' + str(m.atom_num - 3) + ' ' + str(m.atom_num - 7) + ' 6')
            
            # angles to the next unit
            if i < (len(chain) - 1):
                if chain[i + 1] == 'M':        
                    m._angles.append(Angle(m.atom_num - 2, m.atom_num - 3, m.atom_num + 1, 92, 104, 8))
                    par_an_out.append(str(m.atom_num - 2) + ' ' + str(m.atom_num - 3) + ' ' + str(m.atom_num + 1) + ' 7')
                if chain[i + 1] == 'S':        
                    m._angles.append(Angle(m.atom_num - 2, m.atom_num - 3, m.atom_num + 1, 120, 50, 9))
                    par_an_out.append(str(m.atom_num - 2) + ' ' + str(m.atom_num - 3) + ' ' + str(m.atom_num + 1) + ' 8')
        
    m.write_itp('sma_top_' + str(seq_id + 1) + '.itp')

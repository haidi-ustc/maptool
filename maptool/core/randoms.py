#!/usr/bin/env python
import os
import time
import random
import numpy as np
from pymatgen.io.cif import CifWriter
from maptool.util.utils import wait,wait_sep,multi_structs
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from pymatgen import Structure,Molecule,Composition
from spglib import get_symmetry_dataset

from pyxtal.symmetry import get_symbol_and_number
from pyxtal.crystal import (random_crystal, 
                           random_crystal_1D, 
                           random_crystal_2D, 
                           random_cluster )

def random_operation(choice):
    assert choice in ["1","2","3","4"]
    if choice=="1":
        return random_structure()
    elif choice=="2":
        pass
    elif choice=="3":
        pass
    else:
        pass

def get_random_spg(dim):
    if dim==3:
        return np.random.randint(1,231)
    elif dim==2:
        return np.random.randint(1,81)
    elif dim==1:
        return np.random.randint(1,76)
    else:
        return np.random.randint(1,57)

def get_random_structure(elem,num_atom,sg,dim,thickness=0):
    max_try=100
    factor=1.0
    i=0
    while i<max_try:
        random.seed(time.time()*1e8)
        if sg==0:
           sg=get_random_spg(dim)
        symbol, sg= get_symbol_and_number(sg,dim)
        if dim == 3:
            rand_crystal = random_crystal(sg, elem, num_atom, factor)
        elif dim == 2:
            rand_crystal = random_crystal_2D(sg, elem, num_atom, thickness, factor)
        elif dim == 1:
            rand_crystal = random_crystal_1D(sg, elem, num_atom, thickness, factor)
        if dim == 0:
            rand_crystal = random_cluster(sg, elem, num_atom, factor)
        if rand_crystal.valid:
           comp = str(rand_crystal.struct.composition)
           comp = comp.replace(" ", "")
           if dim > 0:
               outpath =  comp + '.cif'
               CifWriter(rand_crystal.struct, symprec=0.1).write_file(filename = outpath)
               ans = get_symmetry_dataset(rand_crystal.spg_struct, symprec=1e-1)['international']
               print('Symmetry requested: {:d}({:s}), generated: {:s}'.format(sg, symbol, ans))
               return True
           else:
               outpath =  comp + '.xyz'
               rand_crystal.to_file(filename = outpath, fmt='xyz')
               ans = PointGroupAnalyzer(rand_crystal.molecule).sch_symbol
               print('Symmetry requested: {:d}({:s}), generated: {:s}'.format(sg, symbol, ans))
               return True

        i+=1

def random_structure():
    print("Input the formula, space group and dimension")
    print("Ff the space group sets to 0, it will use random one")
    print("the dimension can be 0,1,2 and 3")
    print("(Fe3O4)4 20 3")
    print("As for 2D or 1D, extral number should be supplied,")
    print("which stands for thickness for 2D or area for 1D")
    print("Si4O8 20 2 2.0")
    wait_sep()
    in_str=wait()
    if in_str=="0":
        return None
    in_str=in_str.split()
    
    comp=Composition(in_str[0])
    elem=[el.symbol for el in comp]
    num_atom=[int(comp[el]) for el in elem]
    sg=int(in_str[1])
    dim=int(in_str[2])
    
    if dim<=2:
       thickness=float(in_str[3])
       return get_random_structure(elem,num_atom,sg,dim,thickness=thickness)
    else:
       return get_random_structure(elem,num_atom,sg,dim)

if __name__=="__main__":
    random_structure()

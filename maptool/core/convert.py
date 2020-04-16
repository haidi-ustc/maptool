#!/usr/bin/env python

import numpy as np
from maptool.io.read_structure import pmg2ase,read_structures
from maptool.util.utils import wait,wait_sep,multi_structs
from dpdata import System
from maptool import NAME
from ase.io import write

def structure2system(structure):
    atom_numbs=[]
    atom_names=[]
    comp=structure.composition
    symbol=sorted(structure.symbol_set)
    el_map=dict(zip(symbol,range(len(symbol)))) 
    atom_types=[el_map[ii.symbol] for ii in structure.species]
    cell=structure.lattice.matrix
    coords=structure.cart_coords
    for el in structure.symbol_set:
        atom_numbs.append(int(comp[el]))
        atom_names.append(el)
    data={
        "atom_numbs":atom_numbs,
        "atom_names":atom_names,
        "atom_types":atom_types,
        "cells":cell.reshape(1,3,3),
        "coords":coords.reshape(1,sum(atom_numbs),3),
        "orig":np.array([0,0,0])
            }
    return System(data=data)

def covert_operation():
    structs,fnames=read_structures()
    if structs is None:
        return None
    multi_structs(structs,fnames)
    print('input the target file format')
    print("supported format: vasp lammps xsf cif nc json yaml xyz ...")
    print("for more information see ASE and Pymatgen manual.")
    wait_sep()
    fmt=wait()

    flag=True
    if fmt.lower()=='vasp':
        fmt='poscar'

    for struct,fname in zip(structs,fnames):
        filename=NAME+'_'+fname.lower()+'.'+fmt

        if  "POSCAR" in fname and flag:
            # it is a bug in pymatgen
            print("The string POSCAR in filename will be replaced by poscar")
            flag=False

        if fmt=='lammps':
            system=structure2system(struct)
            system.to_lammps_lmp(filename)
            continue

        try:
            print("Write {0:20s} {1:s}".format(fname,'Pymatgen:Struture'))
            struct.to(fmt,filename)
        except:
            print("Write {0:20s} {1:s} ".format(fname,'Atoms:Atoms'))
            atoms=pmg2ase(struct)
            write(filename)
    return True

if __name__=="__main__":
    from pymatgen import Structure
    structure=Structure.from_file('POSCAR')
    system=structure2system(structure)
    print(system)
    print('--')
    covert_operation([structure],["POSCAR"],'lammps')


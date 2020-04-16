#!/usr/bin/env python
from maptool import NAME
from maptool.util.utils import sepline,multi_structs
from maptool.io.read_structure import read_structures
from pymatgen import Structure,Molecule
from pymatgen.symmetry.analyzer import PointGroupAnalyzer,SpacegroupAnalyzer
from pymatgen.analysis.molecule_matcher import MoleculeMatcher

def structure_symmetry():
    structs,fnames=read_structures()
    multi_structs(structs,fnames)
    for struct,fname in zip(structs,fnames):
        if isinstance(struct,Structure):
            sa=SpacegroupAnalyzer(struct)
            print("file name: {}".format(fname))
            print("{} : {}".format('Structure Type','periodicity'))
            print("{} : {}".format('Lattice Type',sa.get_lattice_type()))
            print("{} : {}".format('Space Group ID',sa.get_space_group_number()))
            print("{} : {}".format('International Symbol',sa.get_space_group_symbol()))
            print("{} : {}".format('Hall Symbol',sa.get_hall()))
            sepline()
        if isinstance(struct,Molecule):
            print("file name: {}".format(fname))
            sa=PointGroupAnalyzer(struct)
            print("{} : {}".format('Structure Type','non-periodicity'))
            print("{} : {}".format('International Symbol',ast.get_pointgroup()))
    return True

def get_primitive_cell():
    structs,fnames=read_structures()
    multi_structs(structs,fnames)
    for struct,fname in zip(structs,fnames):
        sepline(ch='Primitive Cell',sp='-')
        ast=SpacegroupAnalyzer(struct)
        prim_st=ast.find_primitive()
        print(prim_st)
        sepline()
        print('save to '+NAME+'_primitive_'+fname+'.vasp')
        prim_st.to(filename=NAME+'_primitive_'+fname+'.vasp',fmt='poscar')
        sepline()
    return True

def get_conventional_cell():
    structs,fnames=read_structures()
    multi_structs(structs,fnames)
    for struct,fname in zip(structs,fnames):
        sepline(ch='conventional  cell',sp='-')
        ast=SpacegroupAnalyzer(struct)
        conv_st=ast.get_conventional_standard_structure()
        print(conv_st)
        sepline()
        print('save to '+NAME+'_convention_'+fname+'.vasp')
        conv_st.to(filename=NAME+'_conventional_'+fname+'.vasp',fmt='poscar')
        sepline()
    return True

def structure_finger_print():
    return None

def structures_difference(distance_tolerance=0.1,rcut=30):
    return None
       
def distance(struct1,struct2,rcut,pbc=False,):
    return None

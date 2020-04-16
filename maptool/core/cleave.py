#!/usr/bin/env python

import os
import numpy as np
from pymatgen import Structure,Molecule
from maptool.util.utils import wait,wait_sep,multi_structs
from maptool.io.read_structure import read_structures,ase2pmg
from pymatgen.core.surface import generate_all_slabs, SlabGenerator
from ase.build.tube import nanotube
from ase.io import read,write

def generate_selected_slab(structs,fnames,miller_index,min_slab_size,min_vac_size):
    for struct,fname in zip(structs,fnames):
        slab=SlabGenerator(struct,miller_index,min_slab_size=min_slab_size,min_vacuum_size=min_vac_size,lll_reduce=True)
        slab_struct=slab.get_slab()
        slab_struct.sort()
        miller_str=[str(i) for i in miller_index]
    
        filename='_'.join(miller_str)+"_"+fname+'.vasp'
        slab_struct.to(filename=filename,fmt='POSCAR')
    
def generate_all_slab(structs,fnames,max_index,min_slab_size,min_vac_size):
    for struct,fname in zip(structs,fnames):
        slabs=generate_all_slabs(struct,max_index=max_index,min_slab_size=min_slab_size,min_vacuum_size=min_vac_size,lll_reduce=True)
        for slab_struct in slabs:
            slab_struct.sort()
            miller_str=[str(i) for i in slab_struct.miller_index]
            filename='_'.join(miller_str)+"_"+fname+'.vasp'
            slab_struct.to(filename=filename,fmt='POSCAR')

def generate_shell(structs,fnames,center_atom,radius,shell=None,vacuum=15):
    for struct,fname in zip(structs,fnames):
        center_coord=struct[center_atom].coords
        if shell is None:
            sites=struct.get_neighbors_in_shell(center_coord,0,radius)
            file_name="sphere_"+fname+".vasp"
        else:
            sites=struct.get_neighbors_in_shell(center_coord,radius,shell)
            file_name="shell_"+fname+".vasp"
        coords=[site[0].coords for site in sites]
        species=[site[0].specie for site in sites]
        mol=Molecule(coords=coords,species=species)
        max_dist=np.max(mol.distance_matrix)
        a=b=c=max_dist+vacuum
        box_struct=mol.get_boxed_structure(a,b,c)
        box_struct.to(filename=file_name,fmt='poscar')

def cleave_operation(choice):
    # cleava surface
    structs,fnames=read_structures()
    multi_structs(structs,fnames)
    for struct in structs:
        if isinstance(struct,Molecule):
           print("cleave operation is only supported for periodic structure, skip !!")
    if choice=="1":
        print("Input the miller index, minimum size in angstroms of layers containing atomssupercell")
        print("and Minimize size in angstroms of layers containing vacuum like this:")
        print('1 0 0 | 5 | 5')
        print('it means miller index is [1,0,0]')
        print("min_slab_size is 5 Ang ")
        print("min_vacum_size is 5 Ang ")
        print("or like this : ")
        print('2 | 5 | 5')
        print('it will generate all slab with miller index less than 2')
        wait_sep()
        in_str=wait()
        len_para=len(in_str.split('|')[0].split())
        if len_para==3:
           tmp_list=in_str.split('|')
           miller_index=[int(x) for x in tmp_list[0].strip().split() ]
           min_slab_size=float(tmp_list[1])
           min_vac_size=float(tmp_list[2])
           generate_selected_slab(structs,fnames,miller_index,min_slab_size,min_vac_size)
           return True
        elif len_para==1:
           tmp_list=in_str.split('|')
           max_index=int(tmp_list[0])
           min_slab_size=float(tmp_list[1])
           min_vac_size=float(tmp_list[2])
           generate_all_slab(structs,fnames,max_index,min_slab_size,min_vac_size)
           return True
        else:
           print("unknow format")
           return None
 
    #cleave sphere
    elif choice=="2":
        print("Input the center atom index, sphere radius and vacuum layer thickness")
        print('1 3.5 15')
        print('it means the sphere will be selected according to the 1st atom')
        print("with the radius equals 5Ang, and vacuum layer thickness is 15 Ang")
        wait_sep()
        in_str=wait()
        para=in_str.split()
        center_atom=int(para[0])-1 
        radius=float(para[1])
        vacuum=float(para[2])
        generate_shell(structs,fnames,center_atom,radius,shell=None,vacuum=vacuum)
        return True

    #cleave shell
    elif choice=="3":
        print("Input the center atom index, start radius, shell thickness and")
        print("vacuum layer thickness")
        print('1 5 10  15')
        print('it means the ball shell will be selected according to the 1st atom')
        print("with the 5< r <15Ang, and vacuum layer thickness is 15 Ang")
        wait_sep()
        in_str=""
        while in_str=="":
           in_str=input().strip()
        para=in_str.split()
        center_atom=int(para[0])-1 
        radius=float(para[1])
        shell=float(para[2])
        vacuum=float(para[3])
        generate_shell(structs,fnames,center_atom,radius,shell=shell,vacuum=vacuum)
        return True

    else:
        print("unkown choice")
        return None


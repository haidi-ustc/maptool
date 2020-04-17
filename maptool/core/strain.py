#!/usr/bin/env python
import os
import numpy as np
from pymatgen.core import Structure,Molecule,Lattice
from maptool.util.utils import wait,wait_sep,multi_structs,your_choice
from maptool.io.read_structure import read_structures,ase2pmg

from maptool.core.selection import atom_selection
from maptool.io.data_io import DataIO
from pymatgen.core.surface import generate_all_slabs, SlabGenerator
from ase.io import read,write
from  goto import with_goto

@with_goto
def strain_operation():
    structs,fnames=read_structures()
    multi_structs(structs,fnames)
    for struct in structs:
        if isinstance(struct,Molecule):
           print("cleave operation is only supported for periodic structure, skip !!") 
    label .input
    print('input the strain component like :')
    print("0.01")
    print("it means aplly a strain of 1% along all directions")
    print("or")
    print("0.01 0.0 0.0")
    print("it means apply a strain of 1% along the x direction")
    print("or")
    print("0.01:0.03:5 0.0 0.0")
    print("it means to devide strain range into 5 parts")
    wait_sep()
    strain_str=wait()
    if strain_str=="0":
       return None
    tmp_list=strain_str.split(" ")
    if len(tmp_list)==1:
        strain=[[float(tmp_list[0]),float(tmp_list[0]),float(tmp_list[0])]]
    elif len(tmp_list)==3 and not ":" in strain_str:
        strain=[[float(x) for x in tmp_list]]
    elif len(tmp_list)==3 and strain_str.count(":")%2==0:
        tmp1=[float(x) for x in tmp_list[0].split(":")]
        tmp2=[float(x) for x in tmp_list[1].split(":")]
        tmp3=[float(x) for x in tmp_list[2].split(":")]
        strain=[]
        if len(tmp1)==3:
            x_range=np.linspace(tmp1[0],tmp1[1],int(tmp1[2]))
        else:
            x_range=tmp1

        if len(tmp2)==3:
            y_range=np.linspace(tmp2[0],tmp2[1],int(tmp2[2]))
        else:
            y_range=tmp2
        if len(tmp3)==3:
            z_range=np.linspace(tmp3[0],tmp3[1],int(tmp3[2]))
        else:
            z_range=tmp3
        for ix in x_range:
            for iy in y_range:
                for iz in z_range:
                        strain.append([ix,iy,iz])
    else:
        print("Unknow format!!!")
        goto .input
    generate_strain_structure(structs,fnames,strain)
    return True

def generate_strain_structure(structs,fnames,strain):
    for struct,fname in zip(structs,fnames):
        i_count=0
        for i_strain in strain:
            outfile_name='strain_'+str(i_count)+'_'+fname+'.vasp'
            struct_cp=struct.copy()
            struct_cp.apply_strain(i_strain)
            struct_cp.to(filename=outfile_name,fmt='poscar')
            i_count+=1
        fis_name='index_strain_'+fname+'.dat'
        data=np.hstack((np.array([range(len(strain))]).T,np.array(strain))) 
        ret=DataIO(data, col_head=['#index','eps_x','eps_y','eps_z'], fmt_all ="%4d %7.4f %7.4f %7.4f\n")
        ret.write(fis_name)

if __name__=="__main__":
    from pymatgen import Structure
    structure=Structure.from_file('POSCAR')
    strain=[[0.01,0.0,0.0],[0.02,0.0,0.0]]
    generate_strain_structure([structure,structure],["POSCAR",'POSCAR1'],strain)

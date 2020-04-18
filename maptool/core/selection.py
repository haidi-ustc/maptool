#!/usr/bin/env python
import os
import numpy as np
from pymatgen import  Element, Structure
from maptool.util.utils import wait,wait_sep,warn_tip

#TODO
# select atom by sphere

def atom_selection(struct):
    in_str=_atom_selection()
    if "|" in in_str.strip():
       atom_index_list=parse_range(in_str,struct)
    else:
       in_str=in_str.strip()
       if in_str[0].isalpha():
          atom_index_list=parse_label(in_str,struct)
       else:
          atom_index_list=parse_index(in_str)
    return atom_index_list,in_str

def _atom_selection():
    '''
    select atoms by three different schemes:
    1. by atomic index
    2. by element symbol
    3. by fractional coordinates range
    '''
    print("")
    print("input data according to tips")
    tip="""
select atoms by following ways:
1. atomic index in POSCAR
   i.e. :  1 2 4-8 10 12-30
   i.e. :  1 2 4 8 10 
2. atomic label
   i.e. :  Si  O
3. atomic position
   i.e. :  0 0.5 | 0.2 0.4 | 0.3 0.7
   this means atoms with 0<x<0.5, 
   0.2<y<0.4 and 0.3<z<0.7 will be seleted
   or just specific the z coordinates,
   i.e. :  ||0.3 0.7
   """
    print(tip)
    wait_sep()
    in_str=wait()
    return in_str

def parse_index(in_str):
    atom_index=[]
    tmp_str=in_str.split()
    for i in tmp_str:
        if '-' not in i:
           atom_index.append(int(i))
        else:
           atom_index.extend(range(int(i.split('-')[0]),int(i.split('-')[1])+1))
    return [i-1 for i in atom_index]

def parse_label(in_str,struct):
    atom_label= [Element(elem) for elem in in_str.split()]
    #atom_label= [Element(elem) for elem in in_str.split()[1:]]
    atom_index=[]
    for i, site in enumerate(struct.sites):
        if site.specie in atom_label:
           atom_index.append(i)
    return atom_index

def parse_range(in_str,struct):
    def check_frac(coord,lim):
        con=[False]*3
        for i in range(3):
           con[i]=coord[i]>=lim[i][0] and coord[i]<=lim[i][1]
        if np.all(con):
           return True
        else:
           return False
    coord_range={}
    tmp_str=in_str.split()
    tmp1_str=' '.join(tmp_str)
    tmp2_str=tmp1_str.split("|")
    icount=0
    for i in tmp2_str:

        if i=='':
          coord_range[icount]=''
        else:
          coord_range[icount]=[float(x) for x in i.split()]
        icount+=1
    for key in coord_range.keys():
        if coord_range[key]=='':
           coord_range[key]=[0,1]
    atom_index=[]
    for i, site in enumerate(struct.sites):
        if check_frac(site.frac_coords,coord_range):
           atom_index.append(i)
    return atom_index

 
if __name__=="__main__":
   ret=parse_index("1 2 4-8 10 12-30")
   print(ret)
   str_st="""2
1.0
       14.4218997955         0.0000000000         0.0000000000
        0.0000000000         5.7512998581         0.0000000000
        0.0000000000         0.0000000000         5.0816001892
   Ta    O   Mn
    8   24    4
Direct
     0.162300004         0.329299995         0.751900039
     0.837699996         0.670700005         0.248099985
     0.337699996         0.170700015         0.251899968
     0.662300004         0.829299954         0.748099961
     0.837699996         0.329299995         0.748099961
     0.162300004         0.670700005         0.251899968
     0.662300004         0.170700015         0.248099985
     0.337699996         0.829299954         0.751900039
     0.099399995         0.404499994         0.430600003
     0.900600022         0.595500006         0.569399997
     0.400600022         0.095499996         0.930600003
     0.599399978         0.904499994         0.069400015
     0.900600022         0.404499994         0.069400015
     0.099399995         0.595500006         0.930600003
     0.599399978         0.095499996         0.569399997
     0.400600022         0.904499994         0.430600003
     0.082899999         0.108600001         0.894700029
     0.917099985         0.891399968         0.105300006
     0.417100018         0.391400009         0.394700029
     0.582900015         0.608600032         0.605300018
     0.917099985         0.108600001         0.605300018
     0.082899999         0.891399968         0.394700029
     0.582900015         0.391400009         0.105300006
     0.417100018         0.608600032         0.894700029
     0.258800008         0.129199996         0.583800027
     0.741199942         0.870800035         0.416199973
     0.241199992         0.370800035         0.083800080
     0.758800058         0.629199965         0.916200020
     0.741199942         0.129199996         0.916200020
     0.258800008         0.870800035         0.083800080
     0.758800058         0.370800035         0.416199973
     0.241199992         0.629199965         0.583800027
     0.000000000         0.173500002         0.250000000
    -0.000000000         0.826500019         0.750000000
     0.500000000         0.326499998         0.750000000
     0.500000000         0.673499981         0.250000000
     """
   st=Structure.from_str(str_st,fmt='POSCAR')
   ret=parse_label("Ta  Mn",st)
   print(ret)

   in_str="0 0.5 | 0.2 0.4 | 0.3 0.7"
   ret=parse_range(in_str,st)
   print(ret)

   in_str="0 0.5 || 0.3 0.7"
   ret=parse_range(in_str,st)
   print(ret)
   
   ret=parse_index("1 2 4")
   print(ret)

   #ret,in_str=atom_selection(st)
   #print(ret,in_str)

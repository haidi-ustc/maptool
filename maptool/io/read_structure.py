import os
import glob
from ase.io import read 
from goto import with_goto
from pymatgen import Structure,Molecule
from pymatgen.io.ase import AseAtomsAdaptor
from maptool import mlog
from maptool.util.utils import sepline,wait_sep,wait

def ase2pmg(atoms):
    '''
    convert ase Atoms. obj. to pymatgen Structure obj.

    Args:
        atoms: Atoms obj.

    Returns:
        structure: Structure obj.
    '''
    aaa=AseAtomsAdaptor()
    return  aaa.get_structure(atoms) 

def pmg2ase(structure):
    '''
    convert Structure obj. to ase Atoms. obj. 

    Args:
        structure: Structure obj.

    Returns:
        atoms: Atoms obj.
    '''
    aaa=AseAtomsAdaptor()
    return aaa.get_atoms(structure)
 
@with_goto
def read_structures(tips=None):
    '''
    Parsing input string from linux command line. 

    Args:
        None
    
    '''
    default_tips='''\
Input the structure filename
supported structure format: xsf .vasp POSCAR .nc .json .xyz ...
paramter format, i.e. :
a.vasp
paramter format, i.e. :
a.vasp b.vasp
paramter format, i.e. :
*.cif
paramter format, i.e. :
NaCl[1-2].cif'''
    if tips:
        print(tips)
    else:
        print(default_tips)
    structs=[] 
    wait_sep()
    label .input
    in_str=wait()

    if in_str=="0":
        return None,None

    if '*' in in_str or '[' in in_str or '?' in in_str:
        fnames=sorted(glob.glob(in_str))
        if len(fnames)==0:
            print("Cannot match any files, check your input format")
            goto .input
    else:
        fnames=in_str.split()
        if not os.path.exists(fnames[0]):
            print("Cannot match any files, check your input format")
            goto .input
    #print(fnames)
    structures,_fnames=read_structures_from_files(fnames)
    if len(structures)==0:
       print("Cannot parse file format, check the file concent")
       goto .input
    return structures,_fnames

def read_structures_from_file(fname):
    '''
    read structure according to filename 

    Args:
        fname: (str) input filename

    Returns:
        structure: Structure obj. or None
    '''
    try:
      atoms=read(fname)
      return  ase2pmg(atoms)
    except:
       try:
           return Molecule.from_file(fname)
       except:
           try:
               return Structure.from_file(fname)
           except:
               print("Parsing error: %s"%fname)
               return None

def read_structures_from_files(fnames):
    '''
    read structures according to filename list

    Args:
        fnames: (list) input filename list

    Returns:
        structure: Structure obj. list and filename list
    '''
    structures=[]
    final_fnames=[]
    assert isinstance(fnames,list)
    for fname in fnames:
        structure=read_structures_from_file(fname)
        mlog.debug(structure)
        if structure is not None:
            structures.append(structure)
            if '/' in fname:
               final_fnames.append(fname.replace('/','_'))
            else:
               final_fnames.append(fname)

    return structures,final_fnames

if __name__=='__main__':
  sts,fn=read_structures()
  print("len sts: ",len(sts))
  print("len fn:",len(fn))
  for f,st in zip(fn,sts):
      print(f)
      print(st)

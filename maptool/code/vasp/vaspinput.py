#!/usr/bin/env python
import os
import re
import itertools
import six
import numpy as np

from math import ceil
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.vasp.sets import MITRelaxSet,MITNEBSet,MITMDSet,MPRelaxSet,MPHSERelaxSet,MPStaticSet,\
              MPHSEBSSet,MPNonSCFSet,MPSOCSet,MVLElasticSet,MVLGWSet,MVLSlabSet,MVLGBSet, MVLNPTMDSet
from pymatgen import Structure
from pymatgen.io.vasp.inputs import Poscar,Kpoints,Potcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from monty.io import zopen
from tabulate import tabulate
from monty.json import MSONable
from pymatgen.util.string import str_delimited
from pymatgen.util.io_utils import clean_lines
from pymatgen.electronic_structure.core import Magmom
from maptool.util.utils import *
from maptool import mlog
from maptool.code.vasp.templates import *

def generate_incar(struct,dirname='.',encut=1.5):
    """
    Generate INCAR according to user's choice. Now these templates are available:
      a >>> Optimization calculation
      b >>> SCF calculation
      c >>> BAND structure calculation
      d >>> DOS calculation
      e >>> ELF calculation
      f >>> Bader charge calculation
      g >>> AIMD NPT calculation
      h >>> AIMD NVT calculation
      i >>> Potential calculation
      j >>> Partial charge calculation
      k >>> STM image calculation
      l >>> optical properties calculation
      m >>> Mechanical properties calculation
      n >>> Frequency calculation
      o >>> Transition state calculation
      p >>> Phonopy + vasp DFPT calculation
      q >>> Phonopy + vasp finite difference calculation
    """
    try:
        pots=Potcar.from_file(os.path.join(dirname, "POTCAR"))
        pot_elems=[]
        max_encut_elems=[]
        for pot in pots:
            pot_elems.append(pot.element)
            max_encut_elems.append(pot.PSCTR['ENMAX'])
        struct_elems=[x.value for x in struct.types_of_specie]
        mlog.debug('Element order in POSCAR %s' % (struct_elems))
        mlog.debug('Element order in POTCAR %s' % (pot_elems))
        if struct_elems==pot_elems:
           pass
        else:
           print("The element order in POTCAR conflicts with POSCAR ")
           os._exit() 
    except:
        warn_tip(0,'\nPOTCAR file not found\n')
        max_encut_elems=[500]
    prop_encut=max(max_encut_elems)*encut
    elec_relax1['ENCUT']=prop_encut
   
    tip='\n'+\
    'for every letter you can append another letters for extra parameters\n'+\
    'The corresponding list are:                                      \n'+\
    'a:  SPIN                                                         \n'+\
    'b:  SOC                                                          \n'+\
    'c:  HSE                                                          \n'+\
    'd:  DIPOLE correction                                            \n'+\
    'e:  Electric filed                                               \n'+\
    'f:  Add grid                                                     \n'+\
    'g:  Add Pressure                                                 \n'+\
    'h:  DFT-D2                                                       \n'+\
    'i:  DFT-D3                                                       \n'+\
    'j:  VDW-DF                                                       \n'+\
    'k:  opt-B86                                                      \n'+\
    'l:  opt-B88                                                      \n'+\
    'm:  LDA+U                                                        \n'+\
    '\nFor exmaple: aai means one optimization by condsidering SPIN and \n'+\
    'DFT-D3 correction.                                               \n'

    warn_tip(1,tip)
    
    sepline(ch=' generate INCAR file ',sp='-')
    print("your choice?")
    print('{} >>> {}'.format('a','Optimization calculation'))
    print('{} >>> {}'.format('b','SCF calculation'))
    print('{} >>> {}'.format('c','BAND structure calculation'))
    print('{} >>> {}'.format('d','DOS calculation'))
    print('{} >>> {}'.format('e','ELF calculation'))
    print('{} >>> {}'.format('f','Bader charge calculation'))
    print('{} >>> {}'.format('g','AIMD NPT calculation'))
    print('{} >>> {}'.format('h','AIMD NVT calculation'))
    print('{} >>> {}'.format('i','Potential calculation'))
    print('{} >>> {}'.format('j','Partial charge calculation'))
    print('{} >>> {}'.format('k','STM image calculation'))
    print('{} >>> {}'.format('l','optical properties calculation'))
    print('{} >>> {}'.format('m','Mechanical properties calculation'))
    print('{} >>> {}'.format('n','Frequency calculation'))
    print('{} >>> {}'.format('o','Transition state calculation'))
    print('{} >>> {}'.format('p','Phonopy + vasp DFPT calculation'))
    print('{} >>> {}'.format('q','Phonopy + vasp finite difference calculation'))


    wait_sep()
    in_str=""
    while in_str=="":
        in_str=input().strip()

    choice=in_str[0]
    in_str=''.join(in_str.split())
#    assert choice in range(1,19)
    ref_incar=MITRelaxSet(struct).incar
    #print(ref_incar)
    ref_incar_cite={}
    spin_paras['MAGMOM']=ref_incar['MAGMOM']
    try:
       LDAU_paras['LDAUJ']=ref_incar['LDAUJ']
       LDAU_paras['LDAUL']=ref_incar['LDAUL']
    except:
       LDAU_paras['LDAUJ']=None
       LDAU_paras['LDAUL']=None

    elec_relax1['LREAL']=ref_incar['LREAL']

    def parse_extra_incar(in_str,modify_val=None):
        incar_str=''
        for i in range(1,len(in_str)):
            if in_str[i] !=' ':
               incar,comment=Incar.from_dict(eval(extra_params[in_str[i]]))
            incar_str+=incar.get_string(pretty=True,comment=comment)
        return incar_str
    incar_str=''
    if choice=='a': # Opt calculation

       for dict_paras in basic_paras:
           incar,comment=Incar.from_dict(eval(dict_paras))
           incar_str+=incar.get_string(pretty=True,comment=comment)
       if len(in_str)>1:
          incar_str+=parse_extra_incar(in_str)
 
    elif choice=='b': # SCF calculation
       ion_relax['NSW']=0
       elec_relax1['EDIFF']=ediff_oth
       output_paras['LCHGARG']=True
       output_paras['LWAVE']=True
       for dict_paras in basic_paras:
           incar,comment=Incar.from_dict(eval(dict_paras))
           incar_str+=incar.get_string(pretty=True,comment=comment)
       if len(in_str)>1:
           incar_str+=parse_extra_incar(in_str)

    elif choice=='c': # band structure calculation
       ion_relax['NSW']=0
       elec_relax1['EDIFF']=ediff_oth
       start_paras['ICHARG']=11
       for dict_paras in basic_paras:
           incar,comment=Incar.from_dict(eval(dict_paras))
           incar_str+=incar.get_string(pretty=True,comment=comment)
       if len(in_str)>1:
           incar_str+=parse_extra_incar(in_str)

    elif choice=='d': # DOS calculation
       ion_relax['NSW']=0
       elec_relax1['EDIFF']=ediff_oth
       start_paras['ICHARG']=11
       for dict_paras in basic_paras:
           incar,comment=Incar.from_dict(eval(dict_paras))
           incar_str+=incar.get_string(pretty=True,comment=comment)
       if len(in_str)>1:
           incar_str+=parse_extra_incar(in_str)

    elif choice=='e': # ELF calculatio
       ion_relax['NSW']=0
       elec_relax1['EDIFF']=ediff_oth
       output_paras['LCHGARG']=True
       output_paras['LELF']=True
       for dict_paras in basic_paras:
           incar,comment=Incar.from_dict(eval(dict_paras))
           incar_str+=incar.get_string(pretty=True,comment=comment)
       if len(in_str)>1:
           incar_str+=parse_extra_incar(in_str)

    elif choice=='f': # Bader charge calculation
       ion_relax['NSW']=0
       elec_relax1['EDIFF']=ediff_oth
       output_paras['LCHGARG']=True
       output_paras['LAECHG']=True
       for dict_paras in basic_paras:
           incar,comment=Incar.from_dict(eval(dict_paras))
           incar_str+=incar.get_string(pretty=True,comment=comment)
       if len(in_str)>1:
           incar_str+=parse_extra_incar(in_str)

    elif choice=='g': # AIMD NPT calculation
       basic_paras.append('md_NPT_paras')
       ion_relax['NSW']=md_step
       ion_relax['IBRION']=0
       ion_relax['POTIM']=1
       ion_relax['ISYM']=0
       for dict_paras in basic_paras:
           incar,comment=Incar.from_dict(eval(dict_paras))
           incar_str+=incar.get_string(pretty=True,comment=comment)
       if len(in_str)>1:
           incar_str+=parse_extra_incar(in_str)

    elif choice=='h': # AIMD NVT calculation
       basic_paras.append('md_NVT_paras')
       ion_relax['NSW']=md_step
       ion_relax['IBRION']=0
       ion_relax['POTIM']=1
       ion_relax['ISYM']=0
       ion_relax['ISIF']=2
       for dict_paras in basic_paras:
           incar,comment=Incar.from_dict(eval(dict_paras))
           incar_str+=incar.get_string(pretty=True,comment=comment)
       if len(in_str)>1:
           incar_str+=parse_extra_incar(in_str)

    elif choice=='i': # Potential calculation
       ion_relax['NSW']=0
       elec_relax1['EDIFF']=ediff_oth
       output_paras['LCHGARG']=True
       output_paras['LVTOT']=True
       for dict_paras in basic_paras:
           incar,comment=Incar.from_dict(eval(dict_paras))
           incar_str+=incar.get_string(pretty=True,comment=comment)
       if len(in_str)>1:
           incar_str+=parse_extra_incar(in_str)

    elif choice=='j': # Partial charge calculation
       basic_paras.append('partial_paras')
       ion_relax['NSW']=0
       start_paras['ISTART']=1
       elec_relax1['EDIFF']=ediff_oth
       output_paras['LCHGARG']=True
       for dict_paras in basic_paras:
           incar,comment=Incar.from_dict(eval(dict_paras))
           incar_str+=incar.get_string(pretty=True,comment=comment)
       if len(in_str)>1:
           incar_str+=parse_extra_incar(in_str)
      
    elif choice=='k': # STM image calculation
       basic_paras.append('stm_paras')
       ion_relax['NSW']=0
       start_paras['ISTART']=1
       elec_relax1['EDIFF']=ediff_oth
       output_paras['LCHGARG']=True
       for dict_paras in basic_paras:
           incar,comment=Incar.from_dict(eval(dict_paras))
           incar_str+=incar.get_string(pretty=True,comment=comment)
       if len(in_str)>1:
           incar_str+=parse_extra_incar(in_str)

    elif choice=='l': # optical properties calculation
       basic_paras.append('optics_paras')
       ion_relax['NSW']=0
       start_paras['ISTART']=1
       elec_relax1['EDIFF']=ediff_oth
       for dict_paras in basic_paras:
           incar,comment=Incar.from_dict(eval(dict_paras))
           incar_str+=incar.get_string(pretty=True,comment=comment)
       if len(in_str)>1:
           incar_str+=parse_extra_incar(in_str)
       
    elif choice=='m': # Mechanical properties calculation
       basic_paras.append('stm_paras')
       ion_relax['NSW']=1
       ion_relax['NFREE']=4
       ion_relax['IBRION']=6
       ion_relax['POTIM']=0.015
       elec_relax1['EDIFF']=ediff_oth
       for dict_paras in basic_paras:
           incar,comment=Incar.from_dict(eval(dict_paras))
           incar_str+=incar.get_string(pretty=True,comment=comment)
       if len(in_str)>1:
           incar_str+=parse_extra_incar(in_str)

    elif choice=='n': # Frequency calculation
       ion_relax['NSW']=1
       ion_relax['NFREE']=4
       ion_relax['IBRION']=5
       ion_relax['POTIM']=0.015
       elec_relax1['EDIFF']=ediff_oth
       for dict_paras in basic_paras:
           incar,comment=Incar.from_dict(eval(dict_paras))
           incar_str+=incar.get_string(pretty=True,comment=comment)
       if len(in_str)>1:
           incar_str+=parse_extra_incar(in_str)

    elif choice=='o': # Transition state calculation
       basic_paras.append('neb_paras')
       ion_relax['POTIM']=0
       ion_relax['EDIFFG']=ediffg_neb
       elec_relax1['EDIFF']=ediff_opt
       for dict_paras in basic_paras:
           
           incar,comment=Incar.from_dict(eval(dict_paras))
           incar_str+=incar.get_string(pretty=True,comment=comment)
       if len(in_str)>1:
           incar_str+=parse_extra_incar(in_str)
     
    elif choice=='p': # Phonopy + vasp DFPT calculation
       ion_relax['IBRION']=8
       elec_relax1['EDIFF']=ediff_phon
       basic_paras.append(grid_paras)
       for dict_paras in basic_paras:
           incar,comment=Incar.from_dict(eval(dict_paras))
           incar_str+=incar.get_string(pretty=True,comment=comment)
       if len(in_str)>1:
           incar_str+=parse_extra_incar(in_str)

    elif 'q' == choice: # Phonopy + vasp finite difference calculation
       ion_relax['NSW']=0
       ion_relax['IBRION']=-1
       elec_relax1['EDIFF']=ediff_phon
       basic_paras.append(grid_paras)
       for dict_paras in basic_paras:
           incar,comment=Incar.from_dict(eval(dict_paras))
           incar_str+=incar.get_string(pretty=True,comment=comment)
       if len(in_str)>1:
           incar_str+=parse_extra_incar(in_str)

    else: # left empty for extend
      raise Exception(f"choice '{choice}' not valid!")

    write_file(os.path.join(dirname, "INCAR"),incar_str)

def generate_kpoint(struct,dirname='.'):
    '''
    Generate KPOINTS file according to user's choice. Currently supports:
      1 : automatic k-grid
      2 : Band structure k-path
      3 : HSE06 band structure k-path
      4 : 3D plot k-grid
    '''
    sepline(ch=' generate KPOINTS file ',sp='-')
    print("your choice?")
    print('{} >>> {}'.format('1','automatic k-grid '))
    print('{} >>> {}'.format('2','Band structure k-path'))
    print('{} >>> {}'.format('3','HSE06 band structure'))
    print('{} >>> {}'.format('4','3D plot k-grid'))
    wait_sep()
    in_str=""
    while in_str=="":
        in_str=input().strip()
    choice=int(in_str)
    assert choice in [1,2,3,4]

    if choice==1:
       auto_kgrid(struct,dirname)
    elif choice==2:
       band_structure_kpath(struct,dirname)
    elif choice==3:
       pass
    else:
      generate_3Dkpoints(struct,dirname)

def auto_kgrid(struct,dirname):
    '''
    Generate KPOINTS according given mesh density.
      input the dimensionality and mesh grid density
      dimensionality can be 0D 1D 2D 3D
      500 for low grid density
      1000 for medium grid density
      2000 for high grid density
      3000 for accurate density
      input format: 1 1000

      note: in 1D system, mesh grids for x & y direction = 1
            in 2D system, mesh grids for z direction = 1
    '''
    print(" input the dimensionality and mesh grid density ")
    print(" dimensionality can be 0D 1D 2D 3D")
    print(" 500 for low grid density")
    print(" 1000 for medium grid density")
    print(" 2000 for high grid density")
    print(" 3000 for accurate density")
    print(" input format: 1 1000")
    wait_sep()
    in_str=""
    while in_str=="":
       in_str=input().strip().split()
    data=[int(x) for x in in_str]
    dim=data[0]
    grid_density=data[1]
    kps=Kpoints.automatic_density(struct,grid_density)
    if dim==0:
       kps.kpts=[[1,1,1]]
    if dim==1:
      # 1 for x and y direction
      # e.g. :
      # ""
      # comment
      # M or G
      # 1 1 <density>
      # ""
       kps.kpts[0][0]=1
       kps.kpts[0][1]=1
    if dim==2:
      # 1 for z direction
      # e.g. :
      # ""
      # comment
      # M or G
      # <dens1> <dens2> 1
      # ""
       kps.kpts[0][2]=1
    kps.write_file(os.path.join(dirname, "KPOINTS"))

def generate_3Dkpoints(struct,dirname):
    """
    Generate 3D system KPOINTS file according to given accuracy
    User should input a float number to set the accuracy.
    """
    #print(" input the dimensionality and mesh grid density ")
    tip='''
    Accuracy Levels: (1) Low:    0.04~0.03;
                     (2) Medium: 0.03~0.02; 
                     (2) Fine:   0.02~0.01; 
    '''
    warn_tip(1,tip)
    print("Input KP-Resolved Value (unit: 2*PI/Ang):")
    wait_sep()
    in_str=""
    while in_str=="":
       in_str=input().strip()
    delta_k=float(in_str)
    (lka,lkb,lkc)=struct.lattice.reciprocal_lattice.abc
    ka=ceil(lka/(2*np.pi)/delta_k)
    kb=ceil(lkb/(2*np.pi)/delta_k)
    ka_dis=np.linspace(-0.5,0.5,ka)
    kb_dis=np.linspace(-0.5,0.5,kb)
    kxx,kyy=np.meshgrid(ka_dis,kb_dis)
    kpts=[[i[0], i[1], 0.0] for i in zip(kxx.flat,kyy.flat)]

    tmp_K="3D K-meshs by maptool with : "+ str(ka)+"x"+ str(kb)+ "\n 1\nReciprocal\n 0 0 0 1"
    # initialize a Kpionts instance from template string
    reciprocal_kpoints=Kpoints.from_string(tmp_K)
    reciprocal_kpoints.labels=None
    reciprocal_kpoints.kpts=kpts
    reciprocal_kpoints.num_kpts=ka*kb
    reciprocal_kpoints.kpts_weights=[1.0]*(ka*kb)
    reciprocal_kpoints.write_file(os.path.join(dirname, "KPOINTS"))
#    print(spk)


def band_structure_kpath(struct,dirname,nkpts=30):
    """
    Generate KPOINTS file for band structure calculation via pymatgen's symmetry
    analyzing system.
    """
    #struct=Structure.from_file('POSCAR')
    #ana_struct=SpacegroupAnalyzer(struct)
    #pst=ana_struct.find_primitive()
    # First brillouin zone
    ibz = HighSymmKpath(struct)
    linemode_kpoints = Kpoints.automatic_linemode(nkpts,ibz)
    linemode_kpoints.write_file(os.path.join(dirname, "KPOINTS"))


def hse06_bandstructure_kpoints(struct, nkpts=20, kmesh=(11, 11, 11),
                                is_shift=(0, 0, 0), ibzkpt_file=False):
    '''
    Generate HSE06 bandstructure KPOINTS
    Append high-symmetry path points to the IBZKPT file and set weight of
    all the high-symmetry path points to zero and then write to "KPOINTS"

    High-symmetry path kpoints is saved as a backup file named 'KPOINTS_bak'

    Note: We asssert the IBZKPT file is valid
    '''
    def chunks(lst, n):
        for i in range(0, len(lst), n):
            yield lst[i: i+n]

    hsk = HighSymmKpath(struct)
    sym_kpts = Kpoints.automatic_linemode(nkpts, hsk)
    sym_kpts.write_file("KPOINTS_bak")

    # high-symmetry k-points
    kpts = sym_kpts.kpts
    nsegs = sym_kpts.num_kpts

    # Line mode k-points
    line_kpts = []
    for rng in chunks(kpts, 2):
        start, end = rng
        line_kpts.append(np.linspace(start, end, nsegs))
    ln_kpts_vec = np.array(line_kpts).reshape((-1, 3))
    ln_kpts_wht = np.zeros(ln_kpts_vec.shape[0])
    ln_kpts = np.c_[ln_kpts_vec, ln_kpts_wht]

    if ibzkpt_file:
        kpoints = open('IBZKPT', "r").readlines()
        kpoints[1] = "{:8}\n".format(
                int(kpoints[1].strip()) + ln_kpts.shape[0]
                )
        with open("KPOINTS", "w") as K:
            K.write("".join(kpoints))
            np.savetxt(K, ln_kpts, fmt="%20.14f%20.14f%20.14f%14d")
    else:
        # generate irreducible k-points in the BZ
        ana_struct = SpacegroupAnalyzer(struct)
        ir_kpts_with_wht = ana_struct.get_ir_reciprocal_mesh(mesh=kmesh,
                                                             is_shift=is_shift)
        ir_kpts_vec = np.array([x[0] for x in ir_kpts_with_wht])
        ir_kpts_wht = np.array([x[1] for x in ir_kpts_with_wht])
        ir_kpts = np.c_[ir_kpts_vec, ir_kpts_wht]

        with open("KPOINTS", "w") as K:
            K.write("Generated by MAPTOOL\n")
            K.write("{:8d}\n".format(ln_kpts.shape[0] + ir_kpts.shape[0]))
            K.write("Reciprocal Lattice\n")

            np.savetxt(K, ir_kpts, fmt="%20.14f%20.14f%20.14f%14d")
            np.savetxt(K, ln_kpts, fmt="%20.14f%20.14f%20.14f%14d")


def generate_potcar(struct,dirname='.'):
    """
    Generate POTCAR according to input structure via pymatgen's POTCAR module
    """
    avail_pot =" ".join(Potcar.FUNCTIONAL_CHOICES)
    tip="""
    Available Pseudo-potentials are:
    """
    tip+=avail_pot+'\n'
    warn_tip(1,tip) 
    print('your choice ?')
    wait_sep()
    in_str=""
    while in_str=="":
       in_str=input().strip()
    assert in_str.upper() in avail_pot

    potcar=Potcar([el.value for el in struct.types_of_specie],functional=in_str)  
    potcar.write_file(os.path.join(dirname, "POTCAR"))


def generate_all_input(struct):
    # select one below
    
    sepline(ch=' generate input files ',sp='-')
    print("your choice?")
    inputsets=['MITRelaxSet','MITNEBSet','MITMDSet','MPRelaxSet','MPHSERelaxSet','MPStaticSet',
              'MPHSEBSSet','MPNonSCFSet','MPSOCSet','MVLElasticSet',
              'MVLGWSet','MVLSlabSet','MVLGBSet',
              'MVLNPTMDSet']
    print('{} >>> {}'.format('1 ','MIT Relax Set'))
    print('{} >>> {}'.format('2 ','MIT NEB Set'))
    print('{} >>> {}'.format('3 ','MIT MD Set'))
    print('{} >>> {}'.format('4 ','MP Relax Set'))
    print('{} >>> {}'.format('5 ','MP HSE Relax Set'))
    print('{} >>> {}'.format('6 ','MP None SCF Set'))
    print('{} >>> {}'.format('7 ','MP SOC Set'))
    print('{} >>> {}'.format('8 ','MVL Elastic Set'))
    print('{} >>> {}'.format('9 ','MVL GW Set'))
    print('{} >>> {}'.format('10','MVL Slab Set'))
    print('{} >>> {}'.format('11','MVL GB Set'))
    print('{} >>> {}'.format('12','MVL NPT Set'))
    wait_sep()
    in_str=wait()
    choice=int(in_str)
    selected_set=inputsets[choice-1]
    cmd=selected_set+'(struct)'
    print(cmd)
    outset=eval(cmd)
    outset.incar['NSW']=100
    print('writting the input files !')
    outset.write_input('./input',include_cif=True)


class Incar(dict, MSONable):
    """
    INCAR object for reading and writing INCAR files. Essentially consists of
    a dictionary with some helper functions
    """

    def __init__(self, params=None):
        """
        Creates an Incar object.

        Args:
            params (dict): A set of input parameters as a dictionary.
        """
        super(Incar, self).__init__()
        if params:

            # if Incar contains vector-like magmoms given as a list
            # of floats, convert to a list of lists
            if (params.get("MAGMOM") and isinstance(params["MAGMOM"][0], (int, float))) \
                    and (params.get("LSORBIT") or params.get("LNONCOLLINEAR")):
                val = []
                for i in range(len(params["MAGMOM"])//3):
                    val.append(params["MAGMOM"][i*3:(i+1)*3])
                params["MAGMOM"] = val

            self.update(params)

    def __setitem__(self, key, val):
        """
        Add parameter-val pair to Incar.  Warns if parameter is not in list of
        valid INCAR tags. Also cleans the parameter and val by stripping
        leading and trailing white spaces.
        """
        super(Incar, self).__setitem__(
            key.strip(), Incar.proc_val(key.strip(), val.strip())
            if isinstance(val, six.string_types) else val)

    def as_dict(self):
        d = dict(self)
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d):
        if d.get("MAGMOM") and isinstance(d["MAGMOM"][0], dict):
            d["MAGMOM"] = [Magmom.from_dict(m) for m in d["MAGMOM"]]
        return Incar({k: v for k, v in d.items() if k not in ("@module","@class",'comment')}),\
               [v for k, v in d.items() if k in ('comment')][0]

    def get_string(self, sort_keys=False, pretty=False,comment='#'):
        """
        Returns a string representation of the INCAR.  The reason why this
        method is different from the __str__ method is to provide options for
        pretty printing.

        Args:
            sort_keys (bool): Set to True to sort the INCAR parameters
                alphabetically. Defaults to False.
            pretty (bool): Set to True for pretty aligned output. Defaults
                to False.
        """
        keys = self.keys()
        if sort_keys:
            keys = sorted(keys)
        lines = []
        for k in keys:
            if k == "MAGMOM" and isinstance(self[k], list):
                value = []

                if (isinstance(self[k][0], list) or isinstance(self[k][0], Magmom)) and \
                        (self.get("LSORBIT") or self.get("LNONCOLLINEAR")):
                    value.append(" ".join(str(i) for j in self[k] for i in j))
                elif self.get("LSORBIT") or self.get("LNONCOLLINEAR"):
                    for m, g in itertools.groupby(self[k]):
                        value.append("3*{}*{}".format(len(tuple(g)), m))
                else:
                    # float() to ensure backwards compatibility between
                    # float magmoms and Magmom objects
                    for m, g in itertools.groupby(self[k], lambda x: float(x)):
                        value.append("{}*{}".format(len(tuple(g)), m))

                lines.append([k, " ".join(value)])
            elif isinstance(self[k], list):
                lines.append([k, " ".join([str(i) for i in self[k]])])
            else:
                lines.append([k, self[k]])

        if pretty:
            return comment+str(tabulate([[l[0], "=", l[1]] for l in lines],
                                tablefmt="plain"))
        else:
            return comment+str_delimited(lines, None, " = ") + "\n"

    def __str__(self):
        return self.get_string(sort_keys=True, pretty=False)

    def write_file(self, filename):
        """
        Write Incar to a file.

        Args:
            filename (str): filename to write to.
        """
        with zopen(filename, "wt") as f:
            f.write(self.__str__())

    @staticmethod
    def from_file(filename):
        """
        Reads an Incar object from a file.

        Args:
            filename (str): Filename for file

        Returns:
            Incar object
        """
        with zopen(filename, "rt") as f:
            return Incar.from_string(f.read())

    @staticmethod
    def from_string(string):
        """
        Reads an Incar object from a string.

        Args:
            string (str): Incar string

        Returns:
            Incar object
        """
        lines = list(clean_lines(string.splitlines()))
        params = {}
        for line in lines:
            m = re.match(r'(\w+)\s*=\s*(.*)', line)
            if m:
                key = m.group(1).strip()
                val = m.group(2).strip()
                val = Incar.proc_val(key, val)
                params[key] = val
        return Incar(params)

    @staticmethod
    def proc_val(key, val):
        """
        Static helper method to convert INCAR parameters to proper types, e.g.,
        integers, floats, lists, etc.

        Args:
            key: INCAR parameter key
            val: Actual value of INCAR parameter.
        """
        list_keys = ("LDAUU", "LDAUL", "LDAUJ", "MAGMOM", "DIPOL",
                     "LANGEVIN_GAMMA", "QUAD_EFG", "EINT")
        bool_keys = ("LDAU", "LWAVE", "LSCALU", "LCHARG", "LPLANE", "LUSE_VDW",
                     "LHFCALC", "ADDGRID", "LSORBIT", "LNONCOLLINEAR")
        float_keys = ("EDIFF", "SIGMA", "TIME", "ENCUTFOCK", "HFSCREEN",
                      "POTIM", "EDIFFG", "AGGAC", "PARAM1", "PARAM2")
        int_keys = ("NSW", "NBANDS", "NELMIN", "ISIF", "IBRION", "ISPIN",
                    "ICHARG", "NELM", "ISMEAR", "NPAR", "LDAUPRINT", "LMAXMIX",
                    "ENCUT", "NSIM", "NKRED", "NUPDOWN", "ISPIND", "LDAUTYPE",
                    "IVDW")

        def smart_int_or_float(numstr):
            if numstr.find(".") != -1 or numstr.lower().find("e") != -1:
                return float(numstr)
            else:
                return int(numstr)

        try:
            if key in list_keys:
                output = []
                toks = re.findall(
                    r"(-?\d+\.?\d*)\*?(-?\d+\.?\d*)?\*?(-?\d+\.?\d*)?", val)
                for tok in toks:
                    if tok[2] and "3" in tok[0]:
                        output.extend(
                            [smart_int_or_float(tok[2])] * int(tok[0])
                            * int(tok[1]))
                    elif tok[1]:
                        output.extend([smart_int_or_float(tok[1])] *
                                      int(tok[0]))
                    else:
                        output.append(smart_int_or_float(tok[0]))
                return output
            if key in bool_keys:
                m = re.match(r"^\.?([T|F|t|f])[A-Za-z]*\.?", val)
                if m:
                    if m.group(1) == "T" or m.group(1) == "t":
                        return True
                    else:
                        return False
                raise ValueError(key + " should be a boolean type!")

            if key in float_keys:
                return float(re.search(r"^-?\d*\.?\d*[e|E]?-?\d*", val).group(0))

            if key in int_keys:
                return int(re.match(r"^-?[0-9]+", val).group(0))

        except ValueError:
            pass

        # Not in standard keys. We will try a hierarchy of conversions.
        try:
            val = int(val)
            return val
        except ValueError:
            pass

        try:
            val = float(val)
            return val
        except ValueError:
            pass

        if "true" in val.lower():
            return True

        if "false" in val.lower():
            return False

        return val.strip().capitalize()

    def diff(self, other):
        """
        Diff function for Incar.  Compares two Incars and indicates which
        parameters are the same and which are not. Useful for checking whether
        two runs were done using the same parameters.

        Args:
            other (Incar): The other Incar object to compare to.

        Returns:
            Dict of the following format:
            {"Same" : parameters_that_are_the_same,
            "Different": parameters_that_are_different}
            Note that the parameters are return as full dictionaries of values.
            E.g. {"ISIF":3}
        """
        similar_param = {}
        different_param = {}
        for k1, v1 in self.items():
            if k1 not in other:
                different_param[k1] = {"INCAR1": v1, "INCAR2": None}
            elif v1 != other[k1]:
                different_param[k1] = {"INCAR1": v1, "INCAR2": other[k1]}
            else:
                similar_param[k1] = v1
        for k2, v2 in other.items():
            if k2 not in similar_param and k2 not in different_param:
                if k2 not in self:
                    different_param[k2] = {"INCAR1": None, "INCAR2": v2}
        return {"Same": similar_param, "Different": different_param}

    def __add__(self, other):
        """
        Add all the values of another INCAR object to this object.
        Facilitates the use of "standard" INCARs.
        """
        params = {k: v for k, v in self.items()}
        for k, v in other.items():
            if k in self and v != self[k]:
                raise ValueError("Incars have conflicting values!")
            else:
                params[k] = v
        return Incar(params)

def write_file(filename,in_str):
    """
    Write Incar to a file.

    Args:
        filename (str): filename to write to.
    """
    with zopen(filename, "wt") as f:
        f.write(in_str)

if __name__=='__main__':
   pass

#!/usr/bin/env python
from  maptool.core.function import *
from  maptool.core.convert import  covert_operation
from  maptool.core.build import  build_operation
from  maptool.core.cleave import  cleave_operation
from  maptool.core.strain import  strain_operation
from  maptool.core.twodim import  twod_operation
from  maptool.core.mpdb import    get_mp_banddos, get_mp_structure,\
                                  get_mp_phase_graph, get_mp_properties
from  maptool.core.analysis import  structure_symmetry,get_primitive_cell,\
                                    get_conventional_cell
from  maptool.code.vasp import   vaspinput,vaspout
from  maptool.util.utils import wait,sepline,wait_sep,your_choice,warn_tip, \
                                multi_structs
from  maptool.io.read_structure import read_structures
from  goto import with_goto


@with_goto
def select_function(choice):
    r"""
    submenu for selecting function
    """
# structure operation
    if choice=="a1":
       print('{} >>> {}'.format('1','random structure generating'))
       print('{} >>> {}'.format('2','random perturbation for atom index'))
       print('{} >>> {}'.format('3','random disturbing for lattice matrix'))
       print('{} >>> {}'.format('4','random disturbing for atom position'))
       your_choice()
       in_str=wait()
       if in_str=="0":
          return None
       return random_operation(in_str)

    elif choice=="a2":
       return covert_operation()

    elif choice=="a3":
         print('{} >>> {}'.format('1','build supercell'))
         print('{} >>> {}'.format('2','build nanotube'))
         print('{} >>> {}'.format('3','build absorption configuration'))
         your_choice()
         in_str=wait()
         if in_str=="0":
             return None
         return build_operation(in_str)

    elif choice=="a4":
         print('{} >>> {}'.format('1','cleave surface'))
         print('{} >>> {}'.format('2','cleave sphere cluster'))
         print('{} >>> {}'.format('3','cleave shell structure'))
         your_choice()
         in_str=wait()
         if in_str=="0":
             return None
         return cleave_operation(in_str)

    elif choice=="a5":
         return strain_operation()

    elif choice=='a6':
         print('{} >>> {}'.format('1','build rippled structure'))
         print('{} >>> {}'.format('2','build multi-layered structure'))
         print('{} >>> {}'.format('3','split multi-layered structure'))
         print('{} >>> {}'.format('4','resize vacuum layer'))
         print('{} >>> {}'.format('5','center atomic-layer along z direction'))
         print('{} >>> {}'.format('6','apply strain along different direction'))
         print('{} >>> {}'.format('7','constrain atom in specific range'))
         print('{} >>> {}'.format('8','get a substrate for 2D material (online!!!)'))
         your_choice()
         in_str=wait()
         if in_str=="0":
             return None
         return twod_operation(in_str)

# structure analysis
    elif choice=="b1":
       return structure_symmetry()
    elif choice=="b2":
       return structure_finger_print()
    elif choice=="b3":
       return structures_difference()   
    elif choice=="b4":
       return get_primitive_cell()
    elif choice=="b5":
       return get_conventional_cell()
    elif choice=="b6":
       return get_xrd()

# vasp in/out tools
    elif choice=="c1":
       structs,fnames=read_structures()
       if structs is None:
           return None
       multi_structs(structs,fnames)
       sepline(ch=' prepare intput files ',sp='-')
       your_choice()
       print('{} >>> {}'.format('1','prepare all files automatically'))
       print('{} >>> {}'.format('2','prepare INCAR file'))
       print('{} >>> {}'.format('3','prepare KPOINTS file'))
       print('{} >>> {}'.format('4','prepare POTCAR file'))
       label .input1
       wait_sep()
       choice=wait()

       if choice=="0":
          return None
       elif choice=="1":
          return generate_all_input(structs,fnames)
       elif choice=="2":
          return generate_incar(structs,fnames)
       elif choice=="3":
          return generate_kpoint(structs,fnames)
       elif choice=="4":
          return generate_potcar(structs,fnames)
       else:
          print("unknown choice, check the input")
          goto .input1
    
    elif choice=="cxx":
       sepline(ch=' summary output files ',sp='=')
       print('{} >>> {}'.format('1','describe OUCAR file'))
       print('{} >>> {}'.format('2','describe OSICAR file'))
       print('{} >>> {}'.format('3','describe vasprun.xml file'))
       label .input2
       wait_sep()
       choice=wait()
       if choice=="0":
          return None
       if choice=="1":
          return describe_OUTCAR()
       elif choice=="2":
          return describe_OSICAR()
       elif choice=="3":
          return describe_vasprun()
       else:
          print("unknown choice, check the input")
          goto .input2

    elif choice=="c2":
       sepline(ch=' vasp output analysis ',sp='-')
       print('{} >>> {}'.format('1 ','total density of states'))
       print('{} >>> {}'.format('2 ','projected density of states'))
       print('{} >>> {}'.format('3 ','band structure'))
       print('{} >>> {}'.format('4 ','projected band structure'))
       print('{} >>> {}'.format('5 ','select one band structure'))
       print('{} >>> {}'.format('6 ','charge density'))
       print('{} >>> {}'.format('7 ','spin density'))
       print('{} >>> {}'.format('8 ','charge density difference'))
       print('{} >>> {}'.format('9 ','spin density component: up/down'))
       print('{} >>> {}'.format('10','average charge density/potential'))
       print('{} >>> {}'.format('11','optics analysis'))
       print('{} >>> {}'.format('12','mechanical analysis'))
       print('{} >>> {}'.format('13','ab initio molecular dynamics analysis'))
       label .input3
       wait_sep()
       choice=wait()
       if choice=="0":
          return None
       if choice=="1":
          return total_dos()
       elif choice=="2":
          return projected_dos()
       elif choice=="3":
          return band_structure()
       elif choice=="4":
          return projected_band_structure()
       elif choice=="5":
          return select_one_band_structure()
       elif choice=="6":
          return charge_density()
       elif choice=="7":
          return spin_density()
       elif choice=="8":
          return charge_density_diff()
       elif choice=="9":
          return spin_density_component()
       elif choice=="10":
          return chg_locp_average()
       elif choice=="11":
          return optics_analysis()
       elif choice=="12":
          return elastic_analysis()
       elif choice=="13":
          return aimd_analysis()
       else:
          print("unknown choice, check the input")
          goto .input3

# online exctraction
    elif choice=="e1":
       return get_mp_banddos()
    elif choice=="e2":
       return get_mp_structure()
    elif choice=="e3":
       return get_mp_properties()
    elif choice=="e4":
       return get_mp_phase_graph()
    elif choice=="88":
         os._exit(0)
    else:
       print("unknown choice, return now")
       return None


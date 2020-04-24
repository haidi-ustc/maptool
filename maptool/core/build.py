#!/usr/bin/env python
import os
from maptool import mlog
from maptool.util.utils import wait,wait_sep,multi_structs
from maptool.io.read_structure import read_structures,ase2pmg,read_structures_from_file
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.core.surface import generate_all_slabs, SlabGenerator
from pymatgen import Structure,Molecule
from ase.build.tube import nanotube
from ase.io import read,write

def build_operation(choice):
    assert choice in ["1","2","3"]
    if choice=="1":
        structs,fnames=read_structures()
        multi_structs(structs,fnames)
        wait_sep()
        tip=        """
Several options are available:

a. A full 3x3 scaling matrix defining the linear combination
   the old lattice vectors. E.g., 2 1 0  0 1 0  0 0 3
   generates a new structure with lattice vectors a' =
   2a + b, b' = 3b, c' = c where a, b, and c are the lattice
   vectors of the original structure.
b. An sequence of three scaling factors. E.g., 2 1 1
   specifies that the supercell should have dimensions 2a x b x
   c.
c. A number, which simply scales all lattice vectors by the
   same factor.
        """
        print(tip)
        wait_sep()
        in_str=wait()
        scaling_list=[int(x) for x in in_str.split()]
        print("scaling list:")
        print(scaling_list)
        for struct,fname in zip(structs,fnames):
            if len(scaling_list)==1:
                scales=scaling_list[0]
                sufix=[scales]
            elif len(scaling_list)==3:
                scales=scaling_list
            elif len(scaling_list)==9:
                scales=[scaling_list[0:3],scaling_list[3:6],scaling_list[6:9]]
            struct_cp=struct.copy()
            struct_cp.make_supercell(scales)
            fname='maptool_SC_'+fname+'.vasp'
            struct_cp.to(filename=fname,fmt='poscar')
        return True
    elif choice=="2":
        print('Only support for CNT now !')
        print('Input the n and m for tube')
        print('Paramter format, i.e. :')
        print('3 3')
        wait_sep()
        in_str=wait()
        m,n=[int(i) for i in in_str.split() ]
        atoms=nanotube(m,n,vacuum=15)
        struct=ase2pmg(atoms)
        struct.to('POSCAR','CNT_'+str(m)+'-'+str(n)+'.vasp')
        return True
    else:
        data={'max_index': 2, 'min_vacum': 20, 'min_slab': 8, 'repeat': [3, 3, 1]}
        def read_adsorb_config(filename):
            with open(filename,'r') as f:
                 datas=f.readlines()
            list_data=[]
            for i in range(len(datas)):
                list_data.append(datas[i][0:datas[i].find('#')].strip().split('='))

            defined_keys=['method','crystal','molecule','max_index','min_vacum','min_slab','repeat']
            data_dict={}
            for key in defined_keys:
                for li in list_data:
                    if key in li[0]:
                       data_dict[key]=li[1]

            data_dict['method']=int(data_dict.get('method').strip())
            data_dict['crystal']=data_dict.get('crystal').strip()
            data_dict['molecule']=data_dict.get('molecule').strip()
            data_dict['max_index']=int(data_dict.get('max_index','1').strip())
            data_dict['min_vacum']=int(data_dict.get('min_vacum','15').strip())
            data_dict['min_slab']=int(data_dict.get('min_slab','5').strip())
            data_dict['repeat']=[int(x) for x in data_dict.get('repeat','1 1 1').strip().split()]
            return data_dict

        def proc_adsorb(cryst,mol,data):
           if data['method'] ==1:
              asf_slab=AdsorbateSiteFinder(cryst)
              ads_sites=asf_slab.find_adsorption_sites()
              ads_structs=asf_slab.generate_adsorption_structures(mol,repeat=data['repeat'])
              for i in range(len(ads_structs)):
                  ads_struct=ads_structs[i]
                  try:
                     miller_str=[str(j) for j in cryst.miller_index]
                  except:
                     miller_str=['adsorb']
                  filename='_'.join(miller_str)+'-'+str(i)+'.vasp'
                  ads_struct.to(filename=filename,fmt='POSCAR')
           else:
              slabs=generate_all_slabs(cryst,max_index=data['max_index'],min_slab_size=data['min_slab'],
                                   min_vacuum_size=data['min_vacum'],lll_reduce=True)
              for slab in slabs:
                  asf_slab=AdsorbateSiteFinder(slab)
                  ads_sites=asf_slab.find_adsorption_sites()
                  ads_structs=asf_slab.generate_adsorption_structures(mol,repeat=data['repeat'])
                  for i in range(len(ads_structs)):
                      ads_struct=ads_structs[i]
                      miller_str=[str(j) for j in slab.miller_index]
                      filename='adsorb'+'_'.join(miller_str)+'-'+str(i)+'.vasp'
                      ads_struct.to(filename=filename,fmt='POSCAR')
 
        filename='adsorb.cfg' 
        if os.path.exists(filename):
           data=read_adsorb_config(filename)
           assert data['method'] in [1,2]
           cryst=read_structures_from_file(data['crystal'])
           mol=read_structures_from_file(data['molecule'])
           proc_adsorb(cryst,mol,data)
        else:
           print('your choice ?')
           print('{} >>> {}'.format('1','read slab from file'))
           print('{} >>> {}'.format('2','build slab by bulk'))
           wait_sep()
           in_str=wait()
           choice=int(in_str)
           assert choice in [1,2]
           data['method']=choice
           tips="""\
Input the structure filename of molecule and substrate
The first file should be molecule and 2nd for crystal
supported structure format: xsf .vasp POSCAR .nc .json .xyz ...
paramter format, i.e. :
mol.xyz POSCAR"""
           structs,fnames=read_structures(tips)
           
           mol=structs[0] 
           mlog.info("read mol from %s"%(fnames[0]))
           mlog.info(mol)
           assert isinstance(mol,Molecule),"the first file should be molecule"
           cryst=structs[1] 
           mlog.info("read crystal from %s"%(fnames[1]))
           mlog.info(cryst)
           assert isinstance(cryst,Structure),"the second file should be crystal"
           proc_adsorb(cryst,mol,data)

        return True

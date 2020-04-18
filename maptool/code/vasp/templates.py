#!/usr/bin/env python3

from maptool import NAME

ediff_opt=1e-5
ediff_oth=1e-6
ediff_phon=1e-7
ediffg   =-0.01
ediffg_neb   =-0.03
md_step=5000

start_paras={
 'comment': '''
#Start parameters for this run
''',
 'SYSTEM' : NAME,
 'NWRITE' : 1,
 'PREC'   :'Accurate',
 'ISART'  : 0,
 'ICHARG' : 2,
}


elec_relax1={
 'comment': '''
#Electronic Relaxation 1
''',
 'ENCUT'  : None,
 'NELM'   : 200,
 'NELMIN' : 6,
 'NELMDL' : -5,
 'EDIFF'  : ediff_opt,
 'LREAL'  : None
}
elec_relax2={
 'comment': '''
#Electronic Relaxation 2
#AMIN     = 0.1
#AMIX     = 0.4
#AMIX_MAG = 1.6
#BMIX     = 1.0
#BMIX_MAG = 1.0
''',
 'ALGO'   : 'Normal',
}

ion_relax={
 'comment': '''
#Ionic relaxation
''',
 'EDIFFG' : -0.01,
 'ISIF'   : 3,
 'IBRION' : 2,
 'POTIM'  : 0.3,
 'ISYM'   : 2,
 'NSW'    : 200
}

pressure_paras={
 'comment': '''
# pressure , unit : Kbar    1Kbar= 0.1 GPa
''',
 'PSTRESS': 10.0
}

dos_paras={
 'comment': '''
# DOS related values
#EMIN     = -20.00
#EMAX     =  20.00
''',
 'ISMEAR': 0,
 'SIGMA' : 0.05,
}

output_paras={
 'comment': '''
# Write flags
''',
 'LWAVE' :False,
 'LCHARG':False,
 'LVTOT' :False,
 'LVHAR' :False,
 'LELF'  :False,
 'LAECHG':False
}

dft_D2={
 'comment': '''
# DFT-D2 correction
''',
 'LVDW'  : True
}

vdw_DF={
 'comment': '''
# DFT-D3 correction
''',
 'GGA'     : 'RE',
 'LUSE_VDW': True,
 'AGGAC'   : 0.0
}

opt_B86={
 'comment': '''
# optB86b-vdw functional  correction
''',
 'GGA'     : 'MK',
 'LUSE_VDW': True,
 'AGGAC'   : 0.0,
 'PARAM1'  : 0.1234,
 'PARAM2'  : 1.0000
}

opt_B88={
 'comment': '''
# optB88-vdw functional correction
''',
 'GGA'     : 'BO',
 'LUSE_VDW': True,
 'AGGAC'   : 0.0,
 'PARAM1'  : 0.18333,
 'PARAM2'  : 0.22
}

paral_paras={
 'comment': '''
# Paralle related parameters
# you have to set a proper value for NPAR ~ sqrt(core)
''',

 'NPAR' : 2
}

hse_paras={
 'comment': '''
# HSE06 related parameters
''',
 'LHFCALC':True,
 'HFSCREEN':0.2,
 'PRECFOCK':'Fast',
 'AEXX'    :0.25,
 'ALGO'    :'All'
}

stm_paras={
 'comment': '''
# STM related parameters
# you have to set a proper value for EINT
''',
 'LPARD'  :True,
 'NBMOD'  :-3,
 'EINT'   :-2
}

partial_paras={
 'comment': '''
# Partial charge related parameters
# you have to set a proper value for IBAND and EINT
''',
 'LPARD'  :True,
 'IBAND'  :32,
 'KPUSE'  :65,
 'LSEPB'  :True,
 'LSEPK'  :True
}

spin_paras={
 'comment': '''
# Spin related parameters
# Guess values are obtained from pymatgen MPRelaxSet
''',
 'ISPIN'  : 2,
 'MAGMOM' : None

}
optics_paras={
 'comment': '''
# Optics related parameters
# you have to set a proper value for NBANDS : ~ 4*NBAND
''',
  'LOPTICS': True,
  'NEDOS'  : 2000,
  'CSHIFT' : 0.1,
  'NBANDS' : 100
}

neb_paras={
 'comment': '''
# NEB related parameters
# you have to set a proper value for IMAGES
''',
  'IMAGES' : 9,
  'ICHAIN' :  0,
  'LCLIMB' : True,
  'SPRING' : -5,
  'IOPT'   : 3
}

dipole_paras={
 'comment': '''
# Dipole correction related parameters
#EPSILON   = 1.0000000  #bulk dielectric constant
''',
 'LDIPOL'  : True,
}

efield_paras={
 'comment': '''
# Electric field related parameters
''',
 'EFIELD'  : 0.5,
 'DIPOL'   : [0.5,0.5,0.5],
 'IDIPOL'  : 1
}

soc_paras={
 'comment': '''
# SOC related parameters
# you have to set a proper value for IMAGES
''',
  'LSORBIT':True
}

LDAU_paras={
 'comment': '''
# LDAU related parameters
# Guess values are obtained from pymatgen MPRelaxSet
''',
 'LDAU'    : True,
 'LDAUJ'   : None,
 'LDAUL'   : None

}

md_NPT_paras={
 'comment': '''
# AIMD related parameters
# you have to set a proper value for TEBEG, TEEND and PMASS
#PMASS     = 10
''',
 'TEBEG'   : 1000,
 'TEEND'   : 1000,
 'NBLOCK'  : 1,
 'KBLOCK'  : 50,
 'SMASS'   : 0,
 'APACO'   : 10,
 'NPACO'   : 500,
 'LANGEVIN_GAMMA_L': 1,
 'LANGEVIN_GAMMA': [10, 10],
 'MDALGO': 3
}

md_NVT_paras={
 'comment': '''
# AIMD related parameters
# you have to set a proper value for TEBEG, TEEND and PMASS
#PMASS     = 10
''',
 'TEBEG'   : 1000,
 'TEEND'   : 1000,
 'NBLOCK'  : 1,
 'KBLOCK'  : 50,
 'SMASS'   : 1,
 'APACO'   : 10,
 'NPACO'   : 500,
}

grid_paras={
 'comment': '''
# parameters for add meshgrid
''',
 'ADDGRID': True
}
basic_paras=['start_paras','elec_relax1','elec_relax2','ion_relax','dos_paras','output_paras']

#
extra_params={'a':'spin_paras',
              'b':'soc_paras',
              'c':'hse_paras',
              'd':'dipole_paras',
              'e':'efield_paras',
              'f':'grid_paras',
              'g':'pressure_paras',
              'h':'dft_D2',
              'i':'dft_D3',
              'j':'vdw_DF',
              'k':'opt_B86',
              'l':'opt_B88',
              'm':'LDAU_paras'
}

#!/usr/bin/env python
# coding: utf-8
import logging
import os

try:
   from ._version import version as __version__
except ImportError:
   from .__about__ import __version__

try:
    from ._date import date as __date__
except ImportError:
    __date__ = 'unkown'

__author__    = "Haidi Wang"
__copyright__ = "Copyright 2018"
__maintainer__= ""
__email__     = "haidi@mail.ustc.edu.cn"
__status__    = "Development"

DEBUG = True
NAME  = 'maptool'

mlog = logging.getLogger(__name__)
mlog.setLevel(logging.DEBUG)
mlogf = logging.FileHandler(os.getcwd()+os.sep+'mpt.log')
mlogf_formatter=logging.Formatter('%(asctime)s - %(levelname)s : %(message)s')
#mlogf_formatter=logging.Formatter('%(asctime)s - %(name)s - [%(filename)s:%(funcName)s - %(lineno)d ] - %(levelname)s \n %(message)s')
mlogf.setFormatter(mlogf_formatter)
mlog.addHandler(mlogf)


modules=['numpy', 'scipy', 'mayavi', 'matplotlib', 'tqdm', 'dpdata',
         'nose', 'coverage', 'spglib', 'pyhull', 'pymatgen', 'qmpy',
         'ase', 'pandas']
try: 
   import pymongo
   PYMONGO=True
except ImportError:
   PYMONGO=False

try:
   from tqdm import tqdm
   TQDM=True
except ImportError:
   TQDM=False

def info():
    """
        Show basic information about """+ NAME + """, its location and version.
        Also information about other libraries used by maptool
        both mandatory and optional
    """

    print(NAME+'\n------------\n')
    print('Version: ' + __version__)
    print('Path:    ' + __path__[0])
    print('Date:    ' + __date__)
    print()
    print('Dependency')
    print('------------')

    import sys
    print('Python version=' + sys.version + '\n')

    if PYMONGO:
        mm = __import__('pymongo')
        print('%10s %10s   %s' % ('pymongo', mm.version, mm.__path__[0]))

    for modui in modules:
        try:
            mm = __import__(modui)
            print('%10s %10s   %s' % (modui, mm.__version__, mm.__path__[0]))
        except ImportError:
            print('%10s %10s Not Found' % (modui, ''))


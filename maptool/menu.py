#!/usr/bin/env python
# coding: utf-8
from goto import with_goto
from maptool.util.utils import sepline,box_center,wait_sep,wait,your_choice
from maptool.submenu import select_function
from maptool.logo import logo
from maptool import NAME

try:
   from ._version import version as __version__
except ImportError:
   from .__about__ import __version__

def head():
    """show the logo and verion info 

    function description:

    Args:
        None

    Returns:
        None

    """
    box_center(ch='-',fill='-',sp="+")
    logo()
    box_center(ch='')
    box_center(__version__)
    box_center(ch='')
    box_center(ch='Written by Wang haidi')
    box_center(ch='URL  https://github.com/haidi-ustc')
    box_center(ch='Bug reports:(haidi@mail.ustc.edu.cn)')
    box_center(ch='-',fill='-',sp="+")

def tail():
    """show the citation info

    function description:

    Args:
        None

    Returns:
        None
    """

    box_center(ch='-',fill='-',sp="+")
    box_center(ch='*BYEBYE*')
    box_center(ch="Thanks for using "+NAME)
    box_center(ch="Have a nice day!!!")
    #box_center(ch="Please cite: Nanoscale 9 (2), 850-855")
    #box_center(ch="https://scholar.google.com/citations?hl=zh-CN&user=9PPScBEAAAAJ")
    box_center(ch='-',fill='-',sp="+")

def structural_operation():
    sepline(ch=" structural operation ",sp='=')
    print('''\
a1 >>> random operation
a2 >>> convert operation
a3 >>> build operation
a4 >>> cleave opeartion
a5 >>> strain operation
a6 >>> 2D structure operation''')

def structural_analysis():
    sepline(ch=" structural analysis ",sp='=')
    print('''\
b1 >>> structure symmetry
b2 >>> structure finger print
b3 >>> structure difference
b4 >>> get primitive cell
b5 >>> get conventional cell
b6 >>> get XRD pattern''')

def vasp_inout():
    sepline(ch=" vasp in/out tools ",sp='=')
    print('''\
c1 >>> prepare input files
c2 >>> analysis output files
c3 >>> summary output files''')

def vasp_workflow():
    sepline(ch=" vasp calclation workflow",sp='=')
    print('''\
d1 >>> optimize structure
d2 >>> calculate band structure
d3 >>> calculate band structure HSE06
d4 >>> calculate dos
d5 >>> calculate dos by HSE06
d6 >>> calculate elastic properties
d7 >>> calculate phonon
d8 >>> execute MD simulation''')

def MP_db():
    sepline(ch=" On-line database ",sp='=')
    print('''\
e1 >>> (MP) get band/dos by mp-ID
e2 >>> (MP) get structure by mp-ID/elements/formula
e3 >>> (MP) get properties by mp-ID
e4 >>> (MP) get phase graph
e5 >>> (OQMD) get structure by ID/elements/formula''',
)

def local_db():
    sepline(ch=" local database ",sp='=')
    print('''\
f1 >>> check local database
f2 >>> get entry by l-ID
f3 >>> get entry by formula
f4 >>> get entry by element
f5 >>> insert entry into database''')

@with_goto
def menu():
    """show the first class menu

    function description:

    Args:
        None

    Returns:
        None
    """

    label .input
    structural_operation()
    structural_analysis()
    vasp_inout()
    vasp_workflow()
    MP_db()
    local_db()
    sepline(ch="=",sp='=')
    print("")
    print("")
    your_choice()
    in_str=wait()
    ret=select_function(in_str) 
    if ret is None:
        goto .input

#TODO
#    sepline(ch=" siesta tools ",sp='=')
#    print('{} >>> {}'.format('si1','prepare input files'))
#    sepline(ch=" gaussian tools ",sp='=')
#    print('{} >>> {}'.format('gi1','prepare input files'))
#    sepline(ch=" nwchem tools ",sp='=')
#    print('{} >>> {}'.format('ni1','prepare input files'))
#    sepline(ch=" quantum espresso tools ",sp='=')
#    print('{} >>> {}'.format('qi1','prepare input files'))
#    sepline(ch=" lammps tools ",sp='=')
#    print('{} >>> {}'.format('li1','prepare input files'))

if __name__=='__main__':
    head()
    menu()
    tail()


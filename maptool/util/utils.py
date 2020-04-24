#!/usr/bin/env python
import os
import numpy as np
from pymatgen import  Element, Structure
from .constants import Len

def sepline(ch='-',sp='-'):
    r'''
    seperate the output by '-'
    '''
    print(ch.center(Len,sp))

def box_center(ch='',fill=' ',sp="|"):
    r'''
    put the string at the center of |  |
    '''
    strs=ch.center(Len,fill)
    print(sp+strs[1:len(strs)-1:]+sp)

def warn_tip(idex=0,instr=''):
    r"""
    wait infomation for waring and tips
    """
    if idex==0:
       stc=" Warning "
    else:
       stc=" Tips "
    box_center(ch=stc,fill='-',sp="+")
    for i_str in instr.split("\n"):
        box_center(ch=' '+i_str+' ',fill=' ',sp=" ")
    box_center(ch='-',fill='-',sp="+")

def wait_sep():
    r'''
    waiting for user input
    '''
    print('--------------->>')

def procs(ch,ind,sp='-->>'):
    r'''
    under processing 
    '''
    if ind !=0:
       nstr=sp+ " ("+str(ind)+') '+ch
    else:
       nstr=sp+"     "+ch
    print(nstr)

def check_file(filename):
    r'''
    checking whether exist the specific file
    '''
    if os.path.exists(filename):
       pass
    else:
       print(filename+" file is not found")
       os._exit(0)

def check_matplotlib():
    r'''
    loading the plot tool: matplot lib
    '''
    try:
       import matplotlib
       matplotlib.use('Agg')
    except ImportError:
       print("you have to install matplotlib")
       os._exit(0)

def wait():
    r"""
    waiting for input
    """
    in_str=""
    while in_str=="":
          in_str=input().strip()
          if in_str=="88":
              os._exit(0)
    return in_str

def your_choice():
    print('0  >>> back')
    print('88 >>> exit')
    print('your choce ?')
    wait_sep()

def multi_structs(structs,fnames):
    if len(structs)>1:
       _fnames='\n'.join(fnames)
       warn_tip(1,"Mutli structures are deteced:\n"+_fnames)

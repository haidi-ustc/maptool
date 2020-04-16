#!/usr/bin/env python
# coding: utf-8
import os
import sys
import string
import copy
from   datetime import datetime
from   optparse import OptionParser, OptionGroup
from   glob import glob
from   maptool.menu import head,menu,tail
from   maptool.util.utils import box_center
from   maptool import info
from   time import time


__author__    = "Haidi Wang"
__copyright__ = "Copyright 2018"
__maintainer__= ""
__email__     = "haidi@mail.ustc.edu.cn"
__status__    = "Development"
__date__      = "May 16, 2018"

def main():

    parser = OptionParser()
    parser.add_option("-i", "--info", dest="info", action="store_true")
    (options, args) = parser.parse_args()

    # shows info
    if options.info:
       info()
       os._exit(0)

    T1=time()
    head()
    menu()
    tail()
    T2=time()
    print("Total Time: %.3f (s) "%(T2-T1))

if __name__=='__main__':
   main()


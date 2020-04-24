import sys,os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from maptool.io.read_structure import   ( ase2pmg, pmg2ase, 
                                          read_structures,
                                          read_structures_from_file, 
                                          read_structures_from_files )
from maptool.core.selection import (parse_range,
                                    parse_index,
                                    parse_label,
                                    parse_sphere)
def setUpModule():
    os.chdir(os.path.abspath(os.path.dirname(__file__)))

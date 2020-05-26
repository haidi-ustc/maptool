import sys,os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from maptool.core.analysis import rdf
from maptool.core.analysis import xrd
from maptool.core.analysis import structure_dedup
from maptool.core.analysis import volume_predict
from maptool.core.oqmd import QMPYRester
from maptool.io.read_structure import read_structures_from_file
from maptool.io.read_structure import read_structures_from_files
def setUpModule():
    os.chdir(os.path.abspath(os.path.dirname(__file__)))

import sys
import os
import unittest
import numpy as np
from pymatgen.io.vasp import Xdatcar

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'core'
from .context import setUpModule
from .context import rdf
from .context import xrd
from .context import read_structures_from_file

class TestRDF(unittest.TestCase):
    def setUp(self):
        self.traj = Xdatcar('test_rdf_XDATCAR.txt').structures
        self.ref_data = np.loadtxt('test_rdf_ref.txt')
        self.threshold = 0.07

    def test_correctness(self):
        traj = self.traj
        ref_data = self.ref_data
        (rdf_vv, radii) = rdf(traj, nbins=100, elem_pair=('V', 'V'))
        (rdf_sese, radii) = rdf(traj, nbins=100, elem_pair=('Se', 'Se'))
        (rdf_sev, radii) = rdf(traj, nbins=100, elem_pair=('Se', 'V'))
        test_data = np.asarray([radii, rdf_sese, rdf_sev, rdf_vv]).T
        diff = np.sum(np.abs((test_data - ref_data).flatten()))
        self.assertGreater(self.threshold, diff)

    def test_invalid_params(self):
        with self.assertRaises(AssertionError):
            rdf(self.traj, r=0)
        with self.assertRaises(AssertionError):
            rdf(self.traj[0:0])
        with self.assertRaises(AssertionError):
            rdf(self.traj, nbins=0)
        with self.assertRaises(AssertionError):
            rdf(self.traj, range=(0, 0))
        with self.assertRaises(Exception):
            rdf(self.traj, elem_pair=('', 'V'))
        with self.assertRaises(Exception):
            rdf(self.traj, elem_pair=('Se', 'H'))
        pass

class TestXRD(unittest.TestCase):
    def setUp(self):
        self.structure = read_structures_from_file('test_xrd.vasp')
        self.threshold = 1e-7
    def test_correctness(self):
        x, y = xrd(self.structure)
        ref_x = np.array([23.87772225, 33.20090362, 41.30842371, 47.66134518, 48.87943686,
                          53.99026029, 60.36971884, 69.69468145, 70.66098376, 74.84482499,
                          76.72042976, 79.41008189, 84.36505928, 86.18940575, 89.73252148])
        ref_y = np.array([ 26.37956616,  23.38287778, 100.        ,  33.07137738,
                           15.34673083,  12.66836049,   9.22679448,  10.48823235,
                           20.14637563,   5.0998731 ,   1.19833501,   4.41542523,
                           25.28893637,  12.07839603,  11.17785816])
        dx = np.sum(np.abs(x - ref_x))
        dy = np.sum(np.abs(y - ref_y))
        self.assertGreater(self.threshold, dx)
        self.assertGreater(self.threshold, dy)
        pass

if '__main__' == __name__:
    unittest.main()

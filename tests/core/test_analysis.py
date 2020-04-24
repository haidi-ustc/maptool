import sys
import os
import unittest
import numpy as np
from pymatgen.io.vasp import Xdatcar

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'core'
from .context import setUpModule
from .context import rdf

THRESHOLD = 0.07

class TestRDF(unittest.TestCase):
    def setUp(self):
        self.traj = Xdatcar('test_rdf_XDATCAR.txt').structures
        self.ref_data = np.loadtxt('test_rdf_ref.txt')

    def test_correctness(self):
        traj = self.traj
        ref_data = self.ref_data
        (rdf_vv, radii) = rdf(traj, nbins=100, elem_pair=('V', 'V'))
        (rdf_sese, radii) = rdf(traj, nbins=100, elem_pair=('Se', 'Se'))
        (rdf_sev, radii) = rdf(traj, nbins=100, elem_pair=('Se', 'V'))
        test_data = np.asarray([radii, rdf_sese, rdf_sev, rdf_vv]).T
        diff = np.sum(np.abs((test_data - ref_data).flatten()))
        self.assertGreater(THRESHOLD, diff)

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


if '__main__' == __name__:
    unittest.main()

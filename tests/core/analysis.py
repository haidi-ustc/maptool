import unittest
import numpy as np
import sys
import os
sys.path.append(os.path.abspath(os.getcwd() + '/../..'))
import maptool
from maptool.core.analysis import rdf
from pymatgen.io.vasp import Xdatcar

THRESHOLD = 0.07


class TestRDF(unittest.TestCase):
    def test_correctness(self):
        traj = Xdatcar('test_rdf_XDATCAR.txt').structures
        ref_data = np.loadtxt('test_rdf_ref.txt')
        (rdf_vv, radii) = rdf(traj, nbins=100, elem_pair=('V', 'V'))
        (rdf_sese, radii) = rdf(traj, nbins=100, elem_pair=('Se', 'Se'))
        (rdf_sev, radii) = rdf(traj, nbins=100, elem_pair=('Se', 'V'))
        test_data = np.asarray([radii, rdf_sese, rdf_sev, rdf_vv]).T
        diff = np.sum(np.abs((test_data - ref_data).flatten()))
        self.assertTrue(diff < THRESHOLD)


if '__main__' == __name__:
    unittest.main()

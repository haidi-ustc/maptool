import sys
import os
import unittest
import hashlib
from glob import glob
import numpy as np
import zipfile
import shutil
from pymatgen.io.vasp import Xdatcar

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'core'
from .context import setUpModule
from .context import rdf
from .context import xrd
from .context import read_structures_from_file
from .context import read_structures_from_files
from .context import structure_dedup


def hash_file(fname: str):
    BLOCK_SIZE = 65536
    file_hash = hashlib.sha256()
    with open(fname, 'rb') as f:
        fb = f.read(BLOCK_SIZE)
        while len(fb) > 0:
            file_hash.update(fb)
            fb = f.read(BLOCK_SIZE)
    return file_hash.hexdigest()


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
        x, y, plot_x, plot_y = xrd(self.structure,
                                   two_theta_range=(0, 90),
                                   fig_name="")
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

    def test_write_file(self):
        fig_name = 'XRD.png'
        peak_raw_fname = 'XRD_peak.txt'
        plot_dat_fname = 'XRD_plot.txt'
        _, _, _, _ = xrd(self.structure,
                         two_theta_range=(0, 120),
                         fig_name=fig_name,
                         peak_raw_fname=peak_raw_fname,
                         plot_dat_fname=plot_dat_fname)
        self.assertTrue(os.path.isfile(fig_name))
        self.assertTrue(os.path.isfile(peak_raw_fname))
        self.assertTrue(os.path.isfile(plot_dat_fname))

        # hash_XRD_png = "ffa0f73b470f8c9531f096d31e4a104b59f52bc7d5d912c42f3cffb71e9ebc16"
        hash_XRD_peak = "6036c73a4f3266b35630c08d22ad19d8c839bfb4aca337c34422b60c20940232"
        hash_XRD_plot = "10fded6ef88ca3d37b4d4c461479d425675ad88a445989b17ed0a67a2f39e041"


        # self.assertEqual(hash_XRD_png, hash_file(fig_name))
        self.assertEqual(hash_XRD_peak, hash_file(peak_raw_fname))
        self.assertEqual(hash_XRD_plot, hash_file(plot_dat_fname))

        os.remove(fig_name)
        os.remove(peak_raw_fname)
        os.remove(plot_dat_fname)


class TestStructureDeduplicate(unittest.TestCase):
    def setUp(self):
        with zipfile.ZipFile('poscars.zip', 'r') as zip_ref:
            zip_ref.extractall(".")
        fnames = glob('poscars/POSCAR*')
        fnames.sort()
        self.structures, self.fnames = read_structures_from_files(fnames)
        self.maxDiff = None  # see the full diff when self.assertEqual fails

    def tearDown(self):
        shutil.rmtree('poscars')

    def test_empty_slist(self):
        slist, flist = structure_dedup([], [])
        self.assertEqual(slist, [])
        self.assertEqual(flist, [])

    def test_empty_flist(self):
        _, flist = structure_dedup(self.structures, [])
        self.assertEqual(flist, [])

    def test_unequal_size_of_slist_and_flist(self):
        with self.assertRaises(AssertionError):
            structure_dedup(self.structures, self.fnames[:4])

    def test_correctness(self):
        slist, flist = structure_dedup(self.structures, self.fnames)
        self.assertEqual(slist[0], self.structures[0])
        flist_ref = ['poscars_POSCAR_1214540',
                     'poscars_POSCAR_1214629',
                     'poscars_POSCAR_1214718',
                     'poscars_POSCAR_1214807',
                     'poscars_POSCAR_1214896',
                     'poscars_POSCAR_1214985',
                     'poscars_POSCAR_1215074',
                     'poscars_POSCAR_1215163',
                     'poscars_POSCAR_1215252',
                     'poscars_POSCAR_1215431',
                     'poscars_POSCAR_1215520',
                     'poscars_POSCAR_1215609',
                     'poscars_POSCAR_1215787',
                     'poscars_POSCAR_1215965',
                     'poscars_POSCAR_1216054',
                     'poscars_POSCAR_1277944',
                     'poscars_POSCAR_18968',
                     'poscars_POSCAR_20426',
                     'poscars_POSCAR_22487']
        self.assertEqual(flist_ref, flist)


if '__main__' == __name__:
    unittest.main()

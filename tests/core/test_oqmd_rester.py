import os
import sys
import unittest
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'core'
from .context import QMPYRester

class TestRester(unittest.TestCase):
    def test_rester_oqmdapi_output_dict(self):
        with QMPYRester() as q:
            kwargs = {'limit': '1'}
            data = q.get_oqmd_phases(verbose=False, **kwargs)
        self.assertTrue(isinstance(data, dict))

    def test_rester_oqmdapi_phase_space_output_dict(self):
        with QMPYRester() as q:
            data = q.get_oqmd_phase_space('Pd-O')
        self.assertTrue(isinstance(data, dict))

    def test_rester_optimade_output_dict(self):
        with QMPYRester() as q:
            kwargs = {'limit': '1'}
            data = q.get_optimade_structures(verbose=False, **kwargs)
        self.assertTrue(isinstance(data, dict))

    def test_rester_oqmdapi_by_id_output(self):
        with QMPYRester() as q:
            data = q.get_oqmd_phase_by_id(fe_id=4061139,fields='name')
        self.assertEqual(data, {'name':'CsHoSiS4'})

    def test_rester_optimade_by_id_output(self):
        with QMPYRester() as q:
            data = q.get_optimade_structure_by_id(id=4061139,fields='id,chemical_formula')
        self.assertEqual(data, {'id': 4061139, 'chemical_formula':'Cs1Ho1S4Si1'})

if '__main__' == __name__:
    unittest.main()

import os
import io
import sys
import unittest
import numpy as np
from ase.io import read
from pymatgen import Structure,Molecule

__package__ = 'mpt_io'
from .context import setUpModule
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from .context import  ( ase2pmg, pmg2ase,
                    read_structures,
                    read_structures_from_file,   
                    read_structures_from_files  )

class TestAse2Pmg(unittest.TestCase):
    print(os.getcwd())
    def test_length(self):
        self.assertEqual(len(self.system_1),len(self.system_2))
    def test_equal(self):
        self.assertEqual(self.system_1,self.system_2)
    def setUp(self):
        self.system_1 = Structure.from_file('POSCAR')
        self.system_2 = ase2pmg(read('POSCAR'))
    def tearDown(self):
        pass


class TestPmg2Ase(unittest.TestCase):
    def setUp(self):
        self.system_1 = read('POSCAR')
        st=Structure.from_file('POSCAR')
        self.system_2 = pmg2ase(st)
        self.places=1e-5
    def test_symbol(self):
        for ii,jj in zip(self.system_1.symbols,self.system_2.symbols):
            self.assertEqual(ii,jj)
    def test_atomic_numbers(self):
        atom_numb1=self.system_1.get_atomic_numbers()
        atom_numb2=self.system_2.get_atomic_numbers()
        for ii,jj in zip(atom_numb1,atom_numb2):
            self.assertEqual(ii,jj)
    def test_positions(self):
        pos1=self.system_1.get_positions()
        pos2=self.system_2.get_positions()
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(pos1[i][j],pos2[i][j],places = self.places)
    def tearDown(self):
        pass

class Std():
    def stub_stdin(self, inputs):
    	stdin = sys.stdin
    	def cleanup():
    		sys.stdin = stdin
    
    	self.addCleanup(cleanup)
    	sys.stdin = io.StringIO(inputs)
    
    def stub_stdout(self):
    	stderr = sys.stderr
    	stdout = sys.stdout
    	def cleanup():
    		sys.stderr = stderr
    		sys.stdout = stdout
    	self.addCleanup(cleanup)
    	sys.stderr = io.StringIO()
    	sys.stdout = io.StringIO()


class ReadStructureTest(unittest.TestCase,Std):

    def setUp(self):
        self.system_ref=Structure.from_file('POSCAR')

    def test_readstructure2(self):
        self.stub_stdin('POSCAR')
        structs,fnames=read_structures(tips='')
        self.assertEqual(len(structs),1)
        self.assertEqual(len(fnames),1)
        self.assertEqual(structs[0],self.system_ref)
       # self.assertEqual(str(sys.stdout.getvalue()), '???')

    def test_readstructure2(self):
        self.stub_stdin('POSCAR*')
        structs,fnames=read_structures(tips='')
        self.assertEqual(len(structs),3)
        self.assertEqual(len(fnames),3)
        for struct,fname in zip(structs,fnames):
            if fname=='POSCAR':
                self.assertEqual(structs[0],self.system_ref)

    def test_readstructure3(self):
        self.stub_stdin('POSCAR[1-2]')
        structs,fnames=read_structures(tips='')
        self.assertEqual(len(structs),2)
        self.assertEqual(len(fnames),2)

    def test_readstructure4(self):
        self.stub_stdin('mol.xyz POSCAR')
        structs,fnames=read_structures(tips='')
        self.assertEqual(len(structs),2)
        self.assertEqual(len(fnames),2)
        self.assertIsInstance(structs[0],Molecule)
        self.assertIsInstance(structs[1],Structure)

if __name__ == '__main__':
       unittest.main()

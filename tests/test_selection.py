import os
import unittest
import numpy as np
from pymatgen import Structure
from context import maptool

from maptool.core.selection import parse_range,parse_index,parse_label
class Compare:
    def test_length(self):
        self.assertEqual(len(self.system_1),len(self.system_2))
    def test_equal(self):
        self.assertEqual(self.system_1,self.system_2)

class TestIndex(unittest.TestCase,Compare):
    def setUp(self):
        self.in_str = "1 4-6 15"
        self.system_1 = parse_index(self.in_str)
        self.system_2 =[0,3,4,5,14]
    def tearDown(self):
        pass

class TestRange1(unittest.TestCase,Compare):
    def setUp(self):
        self.in_str = "0 0.5 || "
        st=Structure.from_file("POSCAR")
        self.system_1 = parse_range(self.in_str,st)
        self.system_2 =[0,2,3,4]
    def tearDown(self):
        pass

class TestRange2(unittest.TestCase,Compare):
    def setUp(self):
        self.in_str = "0 0.5 |0.2 0.3| "
        st=Structure.from_file("POSCAR")
        self.system_1 = parse_range(self.in_str,st)
        self.system_2 =[4]
    def tearDown(self):
        pass

class TestLabel(unittest.TestCase,Compare):
    def setUp(self):
        self.in_str = "Mn"
        st=Structure.from_file("POSCAR")
        self.system_1 = parse_label(self.in_str,st)
        self.system_2 =[3,4,5]
    def tearDown(self):
        pass

if __name__ == '__main__':
   unittest.main()

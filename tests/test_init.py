import unittest

class Test_Math(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass


    @classmethod
    def tearDownClass(cls):
        pass

    def test_fun_01(self):
        print("test1")

    def test_fun_02(self):
        print("test2")

if __name__ == '__main__':  
    unittest.main()


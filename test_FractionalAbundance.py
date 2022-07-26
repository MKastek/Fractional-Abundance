import unittest
unittest.TestLoader.sortTestMethodsUsing = None
import time
import numpy as np
from FractionalAbundance import FractionalAbundance


class TestFA(unittest.TestCase):

    atom_lst = ['Mo']

    def test_concurrent(self):
        for i in range(len(TestFA.atom_lst)):
            with self.subTest(i=i):
                t1 = time.time()
                FA_con = FractionalAbundance(atom=TestFA.atom_lst[i], concurrent=True)
                FA_con.calculate()
                t2 = time.time()
                FA = FractionalAbundance(atom=TestFA.atom_lst[i], concurrent=False)
                t3 = time.time()
                self.assertGreater(t3-t2, t2-t1, f"concurrent is not faster at {TestFA.atom_lst[i]}")

    def test_FA_arr(self):
        for i in range(len(TestFA.atom_lst)):
            with self.subTest(i=i):
                FA = FractionalAbundance(atom=TestFA.atom_lst[i], concurrent=False)
                self.assertIsNotNone(FA.FA_arr, "Test value is not none")

    def test_numba(self):
        for i in range(len(TestFA.atom_lst)):
            with self.subTest(i=i):
                t1 = time.time()
                FA_con = FractionalAbundance(atom=TestFA.atom_lst[i], concurrent=False,  numba=True)
                FA_con.calculate()
                t2 = time.time()
                FA = FractionalAbundance(atom=TestFA.atom_lst[i], concurrent=False)
                t3 = time.time()
                self.assertGreater(np.round(t3-t2,3), np.round(t2-t1,3), f"numba is not faster at {TestFA.atom_lst[i]}")


if __name__ == '__main__':
    unittest.main()
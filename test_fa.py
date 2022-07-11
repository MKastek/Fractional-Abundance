import unittest
import time
from FractionalAbundance import FractionalAbundance


class TestFA(unittest.TestCase):

    atom_lst = ['He', 'Ne', 'Ar', 'Kr', 'Xe']

    def test_concurrent(self):
        for i in range(len(TestFA.atom_lst)):
            with self.subTest(i=i):
                t1 = time.time()
                FA_con = FractionalAbundance(atom=TestFA.atom_lst[i], concurrent=True)
                FA_con.calculate()
                t2 = time.time()
                FA = FractionalAbundance(atom=TestFA.atom_lst[i], concurrent=False)
                t3 = time.time()
                self.assertGreater(t3-t2, t2-t1, "concurrent is not faster")


if __name__ == '__main__':
    unittest.main()
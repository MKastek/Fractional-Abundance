import unittest
import pandas as pd
import os
import time
from FractionalAbundance import FractionalAbundance
import matplotlib.pyplot as plt
unittest.TestLoader.sortTestMethodsUsing = None


class TestFA(unittest.TestCase):

    atom_lst = ['He', 'Ne', 'Ar', 'Kr', 'Xe']

    @classmethod
    def setUpClass(cls):
        cls.time_df = pd.DataFrame()

    @classmethod
    def tearDownClass(cls):
        cls.time_df['element'] = cls.atom_lst
        cls.time_df = cls.time_df.set_index('element')
        cls.time_df.to_csv(os.path.join('tests', 'time-no-numba.csv'))
        cls.time_df.plot.bar(rot=0, title='Computation time [s]').grid()
        plt.savefig(os.path.join('tests','test.png'))

    def test_concurrent(self):
        calc_time_con = []
        calc_time_np = []
        for i in range(len(TestFA.atom_lst)):
            with self.subTest(i=i):
                t1 = time.time()
                FA_con = FractionalAbundance(atom=TestFA.atom_lst[i], concurrent=True)
                t2 = time.time()
                FA = FractionalAbundance(atom=TestFA.atom_lst[i], concurrent=False)
                t3 = time.time()
                # self.assertGreater(t3-t2, t2-t1, f"concurrent is not faster at {TestFA.atom_lst[i]}")
                calc_time_con.append(t2 - t1)
                calc_time_np.append(t3 - t2)

        self.time_df['multi-threading-numba'] = calc_time_con
        self.time_df['numba'] = calc_time_np


if __name__ == '__main__':
    unittest.main()

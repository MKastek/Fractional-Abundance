import unittest
import pandas as pd
import time
from Fractional_Abundace.FractionalAbundance import FractionalAbundance
from Fractional_Abundace.FractionalAbundanceNP import FractionalAbundanceNP
from Fractional_Abundace.FractionalAbundanceCUDA import FractionalAbundanceCUDA
import matplotlib.pyplot as plt
from plasmapy.particles import atomic_number
from pathlib import Path
from timeit import timeit

unittest.TestLoader.sortTestMethodsUsing = None


class TestFA(unittest.TestCase):
    @staticmethod
    def plot_result(time_df):
        time_df.drop("charge-number", axis=1).plot.bar(
            rot=0, title="Computation time [s]"
        ).grid()
        plt.grid(which="minor", alpha=0.3)
        plt.grid(which="major", alpha=0.7)
        plt.legend()
        plt.savefig(Path() / "tests" / "plots" / "test_bar_plot.png")
        plt.clf()

        plt.scatter(time_df["charge-number"], time_df["numba"], label="numba")
        plt.plot(time_df["charge-number"], time_df["numba"], linestyle="--")
        plt.scatter(
            time_df["charge-number"],
            time_df["multi-threading-numba"],
            label="numba with threading",
        )
        plt.plot(
            time_df["charge-number"], time_df["multi-threading-numba"], linestyle="--"
        )
        plt.minorticks_on()
        plt.title("Time comparison: numba vs numba with threading")
        plt.xlabel(r"Charge number")
        plt.ylabel(r"Time (s)")
        plt.grid(which="minor", alpha=0.3)
        plt.grid(which="major", alpha=0.7)
        plt.legend()
        plt.savefig(Path() / "tests" / "plots" / "test_scatter_plot.png")

    @classmethod
    def setUpClass(cls):
        cls.element_list = [
            "He",
            "Li",
            "Ne",
            "Ar",
            "Kr",
            "Mo",
            "Xe",
            "B",
            "Al",
            "Cr",
            "O",
        ]
        cls.element_list = sorted(cls.element_list, key=lambda x: atomic_number(x))
        print(cls.element_list)
        cls.path_to_data = Path() / "data" / "unresolved"
        cls.time_df = pd.DataFrame()

    @classmethod
    def tearDownClass(cls):
        cls.time_df["element"] = cls.element_list
        cls.time_df = cls.time_df.set_index("element")
        cls.time_df.to_csv(Path() / "tests" / "data" / "time.csv")

        cls.plot_result(cls.time_df)

    def test_concurrent(self):
        calc_time_numba_threading = []
        calc_time_numba = []
        calc_time_numpy = []
        calc_time_cuda = []
        charge_number = []
        FA_con = FractionalAbundance(
            element="He", concurrent=True, path_to_data=self.path_to_data
        )
        for i, element in enumerate(TestFA.element_list):
            with self.subTest(i=i):
                trials = 10

                time_numba_threading = (
                    timeit(
                        stmt="FractionalAbundance(element, True)",
                        globals={
                            "FractionalAbundance": FractionalAbundance,
                            "element": element,
                        },
                        number=trials,
                    )
                    / trials
                )

                time_numba = (
                    timeit(
                        stmt="FractionalAbundance(element, False)",
                        globals={
                            "FractionalAbundance": FractionalAbundance,
                            "element": element,
                        },
                        number=trials,
                    )
                    / trials
                )

                time_numpy = (
                    timeit(
                        stmt="FractionalAbundanceNP(element)",
                        globals={
                            "FractionalAbundanceNP": FractionalAbundanceNP,
                            "element": element,
                        },
                        number=trials,
                    )
                    / trials
                )

                time_cuda = (
                    timeit(
                        stmt="FractionalAbundanceCUDA(element)",
                        globals={
                            "FractionalAbundanceCUDA": FractionalAbundanceNP,
                            "element": element,
                        },
                        number=trials,
                    )
                    / trials
                )

                calc_time_numba_threading.append(time_numba_threading)
                calc_time_numba.append(time_numba)
                calc_time_numpy.append(time_numpy)
                calc_time_cuda.append(time_cuda)
                charge_number.append(atomic_number(element))
                print(element)

        self.time_df["multi-threading-numba"] = calc_time_numba_threading
        self.time_df["numba"] = calc_time_numba
        self.time_df["numpy"] = calc_time_numpy
        self.time_df["cuda"] = calc_time_cuda
        self.time_df["charge-number"] = charge_number


if __name__ == "__main__":
    unittest.main()

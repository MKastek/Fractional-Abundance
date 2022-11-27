import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate
from concurrent.futures import ThreadPoolExecutor
import time
from BaseFractionalAbundance import BaseFractionalAbundance
import re
import numba
from numba.core import types
from numba.typed import Dict
from pathlib import Path

float_array = types.float64[:, :]


plt.rcParams["figure.figsize"] = [10.5, 0.65 * 10.5]


class FractionalAbundance(BaseFractionalAbundance):
    def __init__(
        self,
        element,
        concurrent=False,
        path_to_data=Path(__file__).parents[1] / "data" / "unresolved",
    ):

        """

        Parameters
        ----------
        element
        concurrent
        path_to_data
        """

        super().__init__(element, path_to_data)
        self.product_all, self.sum_all = self.calculate_cum_sum_prod(
            self.SCD_matrix, self.ACD_matrix, self.Z
        )
        if not concurrent:
            self.FA_arr = [
                self.get_Fractional_Abundance(
                    ion=ion,
                    product_all=np.array(self.product_all),
                    sum_all=np.array(self.sum_all),
                )
                for ion in range(self.Z + 1)
            ]
        else:
            self.calculate()

    @staticmethod
    @numba.njit(parallel=True)
    def calculate_cum_sum_prod(SCD_matrix, ACD_matrix, Z):
        """

        Parameters
        ----------
        SCD_matrix
        ACD_matrix
        Z

        Returns
        -------

        """
        K = [
            np.divide(
                10 ** SCD_matrix[str(i) + str(i + 1)],
                10 ** ACD_matrix[str(i + 1) + str(i)],
            )
            for i in range(Z)
        ]

        K.insert(0, np.ones_like(K[0]))

        product_all = [K[0]]
        current_product = K[0]
        sum_all = np.zeros_like(K[0])
        for i in range(1, len(K)):
            current_product = np.multiply(K[i], current_product)
            sum_all += current_product
            product_all.append(current_product)

        return product_all, sum_all

    @staticmethod
    @numba.njit(parallel=True)
    def get_Fractional_Abundance(ion, product_all, sum_all):
        """

        Parameters
        ----------
        ion
        product_all
        sum_all

        Returns
        -------

        """
        FA = np.divide(product_all[ion], sum_all)
        return FA

    def worker(self, ion, product_all, sum_all):
        """

        Parameters
        ----------
        ion
        product_all
        sum_all

        Returns
        -------

        """
        fun = self.get_Fractional_Abundance(
            ion=ion, product_all=np.array(product_all), sum_all=np.array(sum_all)
        )
        return fun

    def calculate(self):
        """

        Returns
        -------

        """
        ion_list = list(range(self.Z + 1))
        pool = ThreadPoolExecutor(self.Z + 1)
        pp = [
            pool.submit(
                self.worker, ion=ion, product_all=self.product_all, sum_all=self.sum_all
            )
            for ion in range(len(ion_list))
        ]
        for p in pp:
            self.FA_arr.append(p.result())

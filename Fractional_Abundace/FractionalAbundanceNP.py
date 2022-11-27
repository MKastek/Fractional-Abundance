import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate
from concurrent.futures import ThreadPoolExecutor
from BaseFractionalAbundance import BaseFractionalAbundance
import time
import re
import numba
from numba.core import types
from numba.typed import Dict
from pathlib import Path

float_array = types.float64[:, :]


plt.rcParams["figure.figsize"] = [10.5, 0.65 * 10.5]


class FractionalAbundanceNP(BaseFractionalAbundance):
    def __init__(
        self,
        element,
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
        self.FA_arr = [self.get_Fractional_Abundance(ion) for ion in range(self.Z + 1)]

    def get_Fractional_Abundance(self, ion=None):
        K = [
            np.divide(
                10 ** self.SCD_matrix[str(i) + str(i + 1)],
                10 ** self.ACD_matrix[str(i + 1) + str(i)],
            )
            for i in range(self.Z)
        ]
        K.insert(0, np.ones_like(K[0]))
        K_cumprod = np.cumprod(K, axis=0)
        K_cumprod_sum = np.sum(K_cumprod, axis=0)
        if ion is None:
            return np.divide(K_cumprod, K_cumprod_sum)
        else:
            return np.divide(K_cumprod, K_cumprod_sum)[ion]

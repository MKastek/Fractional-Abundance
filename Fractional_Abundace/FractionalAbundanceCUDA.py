import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate
import time
from numba.core import types
from numba.typed import Dict
from pathlib import Path
from numba import cuda
from BaseFractionalAbundance import BaseFractionalAbundance

float_array = types.float64[:, :]
import math

plt.rcParams["figure.figsize"] = [10.5, 0.65 * 10.5]


class FractionalAbundanceCUDA(BaseFractionalAbundance):
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
        self.K = self.calculate_K()
        self.FA_arr = self.get_Fractional_Abundance()

    def calculate_K(self):
        K = [
            np.divide(
                10 ** self.SCD_matrix[str(i) + str(i + 1)],
                10 ** self.ACD_matrix[str(i + 1) + str(i)],
            )
            for i in range(self.Z)
        ]
        K.insert(0, np.ones_like(K[0]))
        return np.array(K)

    @staticmethod
    @cuda.jit("void(float64[:,:,:], float64[:,:,:])")
    def calculate_FA_gpu(X, Y):
        k, j, i = cuda.grid(3)
        if i < X.shape[2] and j < X.shape[1] and k < X.shape[0]:
            sum_all = 0
            current_product = 1
            for k in range(X.shape[0]):
                current_product *= X[k, j, i]
                sum_all += current_product
                Y[k, j, i] = current_product
                if k == X.shape[0] - 1:
                    for k in range(X.shape[0]):
                        Y[k, j, i] = Y[k, j, i] / sum_all
                        cuda.syncthreads()

    def get_Fractional_Abundance(self):
        threadsperblock = (self.Z + 1, 5, 2)
        blockspergrid_x = math.ceil(self.K.shape[0] / threadsperblock[0])
        blockspergrid_y = math.ceil(self.K.shape[1] / threadsperblock[1])
        blockspergrid_z = math.ceil(self.K.shape[2] / threadsperblock[2])
        blockspergrid = (blockspergrid_x, blockspergrid_y, blockspergrid_z)
        # device array
        device_array = cuda.to_device((self.K).astype(np.float64))
        # result array
        result_array = cuda.device_array_like(np.ones_like(self.K).astype(np.float64))
        # calculate
        self.calculate_FA_gpu[blockspergrid, threadsperblock](
            device_array, result_array
        )
        # copy result
        result_gpu = result_array.copy_to_host()
        return result_gpu

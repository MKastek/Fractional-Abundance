import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate

from numba.core import types
from numba.typed import Dict
from pathlib import Path
from numba import cuda

float_array = types.float64[:, :]
import math

plt.rcParams["figure.figsize"] = [10.5, 0.65 * 10.5]


class FractionalAbundanceCUDA:
    def __init__(
        self,
        element,
        concurrent=False,
        path_to_data=r"C:\Users\marci\Desktop\Projekt NCN\Zadania\1.Styczeń\Fractional_Abundance\data\unresolved",
    ):

        """

        Parameters
        ----------
        element
        concurrent
        path_to_data
        """

        self.element = element
        self.path_to_data = path_to_data
        self.ACD_file, self.SCD_file = self.select_files()

        (
            self.Z,
            self.num_of_Ne_axes,
            self.num_of_Te_axes,
            self.num_of_lines_to_read_with_axes,
            self.sum_of_axes,
        ) = self.read_first_line_of_file(self.SCD_file)
        self.Te, self.Ne = self.read_axes(self.SCD_file)

        self.empty_matrix = np.empty([self.num_of_Te_axes, self.num_of_Ne_axes])
        self.num_of_lines_to_read_with_const_Te = int(np.ceil(self.num_of_Ne_axes / 8))
        self.start_line = int(np.ceil(self.sum_of_axes / 8)) + 3
        self.stop_line = (
            self.start_line
            + self.num_of_lines_to_read_with_const_Te * self.num_of_Te_axes
        )
        self.move = self.stop_line - self.start_line + 1

        self.FA_arr = []
        self.Ne_new, self.Te_new = np.logspace(10, 15, num=100), np.logspace(
            np.log10(5), np.log10(20000), num=800
        )

        self.ACD_matrix = self.read_coefficients_matrices(self.ACD_file, type="ACD")
        self.SCD_matrix = self.read_coefficients_matrices(self.SCD_file, type="SCD")

        self.K = self.calculate_K()
        self.FA_arr = self.get_Fractional_Abundance()

    def select_files(self):
        """

        Returns
        -------

        """
        path_to_data = Path(self.path_to_data)
        ACD_file_path = list(
            path_to_data.glob("acd??_{}*.dat".format(self.element.lower()))
        )[0]
        SCD_file_path = list(
            path_to_data.glob("scd??_{}*.dat".format(self.element.lower()))
        )[0]
        return ACD_file_path, SCD_file_path

    def read_first_line_of_file(self, filepath):
        """

        Parameters
        ----------
        filepath

        Returns
        -------

        """
        with open(filepath) as file:
            first_line = file.readline().strip().split()
            Z, num_of_Ne_axes, num_of_Te_axes = (
                int(first_line[0]),
                int(first_line[1]),
                int(first_line[2]),
            )
            sum_of_axes = num_of_Te_axes + num_of_Ne_axes
            num_of_lines_to_read_with_axes = int(np.ceil(sum_of_axes / 8))
        return (
            Z,
            num_of_Ne_axes,
            num_of_Te_axes,
            num_of_lines_to_read_with_axes,
            sum_of_axes,
        )

    def read_axes(self, filepath):
        """

        Parameters
        ----------
        filepath

        Returns
        -------

        """
        (
            Z,
            num_of_Ne_axes,
            num_of_Te_axes,
            num_of_lines_to_read_with_axes,
            sum_of_axes,
        ) = self.read_first_line_of_file(filepath)

        with open(filepath) as file:
            file.readline()
            file.readline()
            data = []
            while num_of_lines_to_read_with_axes > 0:
                for item in file.readline().strip().split():
                    data.append(item)
                num_of_lines_to_read_with_axes -= 1

        Ne = [10 ** float(item) for item in data[:num_of_Ne_axes]]
        Te = [10 ** float(item) for item in data[num_of_Ne_axes:]]
        return Te, Ne

    def read_data_from(
        self,
        filepath,
        start_line,
        stop_line,
        num_of_lines_to_read_with_const_Te,
        empty_matrix,
    ):
        """

        Parameters
        ----------
        filepath
        start_line
        stop_line
        num_of_lines_to_read_with_const_Te
        empty_matrix

        Returns
        -------

        """

        with open(filepath) as file:
            data = np.array(
                [
                    item.split()
                    for item in file.read().strip().splitlines()[start_line:stop_line]
                ],
                dtype=object,
            )
            data = np.split(data, 1)
            iter = 0
            for i in range(0, int(len(data[0])), num_of_lines_to_read_with_const_Te):
                line_data = np.array([])
                for j in range(num_of_lines_to_read_with_const_Te):
                    line_data = np.concatenate((line_data, data[0][i + j]))
                empty_matrix[iter, :] = line_data
                iter += 1

        f = interpolate.interp2d(
            self.Ne, self.Te, empty_matrix.astype(np.float64), kind="cubic"
        )

        return f(self.Ne_new, self.Te_new)

    def read_coefficients_matrices(self, filepath, type):
        """

        Parameters
        ----------
        filepath
        type

        Returns
        -------

        """

        CD = Dict.empty(key_type=types.unicode_type, value_type=float_array)

        if type == "SCD":
            for i in range(self.Z):
                CD[str(i) + str(i + 1)] = self.read_data_from(
                    filepath,
                    self.start_line + i * self.move,
                    self.stop_line + i * self.move,
                    self.num_of_lines_to_read_with_const_Te,
                    self.empty_matrix.copy(),
                )
        elif type == "ACD":
            for i in range(self.Z):
                CD[str(i + 1) + str(i)] = self.read_data_from(
                    filepath,
                    self.start_line + i * self.move,
                    self.stop_line + i * self.move,
                    self.num_of_lines_to_read_with_const_Te,
                    self.empty_matrix.copy(),
                )
        return CD

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

    def get_Fractional_Abundance(self, ion=None):
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

    def plot_FA_all(self, index_Ne=50):
        """

        Parameters
        ----------
        index_Ne

        Returns
        -------

        """
        for i in range(self.Z + 1):
            x = self.Te_new
            y = self.FA_arr[i][:, index_Ne]
            plt.plot(x, y, label="$" + self.element + "^{" + str(i) + "+}$")
            plt.xscale("log")
            plt.yscale("log")
            plt.ylim((10**-3, 10**0))
            plt.xlim((5, 20000))
            plt.grid()
            plt.title(
                "Fractional Abundance of "
                + self.element
                + " in $N_{e}$  = "
                + "{:.2e}".format(self.Ne_new[index_Ne])
                + " $cm^{-3}$",
                fontsize=16,
            )
            plt.xlabel("$T_{e}$ [eV]", fontsize=16)
            plt.ylabel("FA", fontsize=16)

        plt.show()

    def create_dataset(self, output_filepath=".", filename="fractional_abundance.dat"):
        columns = ["T"]
        for Z in range(self.Z):
            columns.append("Z{}".format(Z + 1))
        FA_output_df = pd.DataFrame(columns=columns)
        for i in range(self.Z):
            FA_output_df[columns[i + 1]] = self.FA_arr[i + 1][:, 40] * 100
        FA_output_df["T"] = self.Te_new
        FA_output_df.to_csv(
            os.path.join(output_filepath, filename), sep=" ", index=False
        )


if __name__ == "__main__":
    path_to_data = r"C:\Users\marci\Desktop\Projekt NCN\Zadania\1.Styczeń\Fractional_Abundance\data\unresolved"
    FA = FractionalAbundanceCUDA(
        element="Xe", concurrent=True, path_to_data=path_to_data
    )
    FA.plot_FA_all()

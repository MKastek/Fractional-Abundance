import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate
from concurrent.futures import ThreadPoolExecutor
import time
import re
import numba
from numba.core import types
from numba.typed import Dict

float_array = types.float64[:,:]


plt.rcParams['figure.figsize'] = [12, 7]


class FractionalAbundance:

    def __init__(self, atom, path_to_data = r"C:\Users\marci\Desktop\Projekt NCN\Zadania\1.Styczeń\Fractional_Abundance\data\unresolved", concurrent = False):

        self.atom = atom
        self.path_to_data = path_to_data
        self.ACD_file, self.SCD_file= self.select_files()

        self.Z, self.num_of_Ne_axes, self.num_of_Te_axes, self.num_of_lines_to_read_with_axes, self.sum_of_axes = self.read_first_line_of_file(self.SCD_file)
        self.Te, self.Ne = self.read_axes(self.SCD_file)

        self.empty_matrix = np.empty([self.num_of_Te_axes, self.num_of_Ne_axes])
        self.num_of_lines_to_read_with_const_Te = int(np.ceil(self.num_of_Ne_axes / 8))
        self.start_line = int(np.ceil(self.sum_of_axes / 8)) + 3
        self.stop_line = self.start_line + self.num_of_lines_to_read_with_const_Te * self.num_of_Te_axes
        self.move = self.stop_line - self.start_line + 1

        self.FA_arr = []
        self.xnew, self.ynew = np.logspace(10, 15, num=100), np.logspace(np.log10(5), np.log10(20000), num=800)

        self.ACD_matrix, self.SCD_matrix = self.read_coefficients_matrices(self.ACD_file,type='ACD'), self.read_coefficients_matrices(self.SCD_file, type='SCD')
        self.product_all, self.sum_all= self.calculate_cum_sum_prod(self.SCD_matrix, self.ACD_matrix, self.Z)

        if not concurrent:
            self.FA_arr = [self.get_Fractional_Abundance(ion = ion, product_all = np.array(self.product_all), sum_all = np.array(self.sum_all)) for ion in range(self.Z +1)]
        else:
            self.calculate()

    def select_files(self):
        filenames = os.listdir(os.path.join(self.path_to_data))
        r_acd = re.compile("acd.*\_{}\.dat".format(self.atom.lower()))
        r_scd = re.compile("scd.*\_{}\.dat".format(self.atom.lower()))
        ACD_file = list(set(list(filter(r_acd.match, filenames))))[0]
        SCD_file = list(set(list(filter(r_scd.match, filenames))))[0]
        return os.path.join(self.path_to_data,ACD_file), os.path.join(self.path_to_data,SCD_file)

    def read_first_line_of_file(self, filepath):
        with open(filepath) as file:
            first_line = file.readline().strip().split()
            Z, num_of_Ne_axes, num_of_Te_axes = int(first_line[0]), int(first_line[1]), int(first_line[2])
            sum_of_axes = num_of_Te_axes + num_of_Ne_axes
            num_of_lines_to_read_with_axes = int(np.ceil(sum_of_axes / 8))
        return Z, num_of_Ne_axes, num_of_Te_axes, num_of_lines_to_read_with_axes, sum_of_axes

    def read_axes(self, filepath):
        Z, num_of_Ne_axes, num_of_Te_axes, num_of_lines_to_read_with_axes, sum_of_axes = self.read_first_line_of_file(filepath)

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

    def read_data_from(self, filepath, start_line, stop_line, num_of_lines_to_read_with_const_Te, empty_matrix):

        with open(filepath) as file:
            data = np.array([item.split() for item in file.read().strip().splitlines()[start_line:stop_line]], dtype=object)
            data = np.split(data, 1)
            iter = 0
            for i in range(0, int(len(data[0])), num_of_lines_to_read_with_const_Te):
                line_data = np.array([])
                for j in range(num_of_lines_to_read_with_const_Te):
                    line_data = np.concatenate((line_data, data[0][i + j]))
                empty_matrix[iter, :] = line_data
                iter += 1

        f = interpolate.interp2d(self.Ne, self.Te, empty_matrix.astype(np.float64), kind='cubic')

        return f(self.xnew, self.ynew)

    def read_coefficients_matrices(self, filepath, type):

        CD =  Dict.empty(
        key_type=types.unicode_type,
        value_type=float_array)

        if type == 'SCD':
            for i in range(self.Z):
                CD[str(i) + str(i + 1)] = self.read_data_from(filepath, self.start_line + i * self.move,
                                                           self.stop_line + i * self.move,
                                                           self.num_of_lines_to_read_with_const_Te,
                                                           self.empty_matrix.copy())
        elif type == 'ACD':
            for i in range(self.Z):
                CD[str(i + 1) + str(i)] = self.read_data_from(filepath, self.start_line + i * self.move,
                                                           self.stop_line + i * self.move,
                                                           self.num_of_lines_to_read_with_const_Te,
                                                           self.empty_matrix.copy())
        return CD

    @staticmethod
    @numba.njit(parallel=True)
    def calculate_cum_sum_prod(SCD_matrix, ACD_matrix, Z):
        K = [np.divide(10 ** SCD_matrix[str(i) + str(i + 1)], 10 ** ACD_matrix[str(i + 1) + str(i)]) for i in range(Z)]

        K.insert(0, np.ones_like(K[0]))

        product_all = [K[0]]
        cur = K[0]
        sum_all = np.zeros_like(K[0])
        for i in range(1, len(K)):
            cur = np.multiply(K[i], cur)
            sum_all += cur
            product_all.append(cur)

        return product_all, sum_all

    @staticmethod
    @numba.njit(parallel=True)
    def get_Fractional_Abundance(ion, product_all, sum_all):
        FA = np.divide(product_all[ion], sum_all)
        return FA

    def plot_FA_all(self, i_Ne=50):
        for i in range(self.Z):
            x = self.ynew
            y = self.FA_arr[i+1][:, i_Ne]
            plt.plot(x, y, label="$" + self.atom + "^{" + str(i) + "+}$")
            plt.xscale("log")
            plt.yscale("log")
            plt.ylim((10 ** -3, 10 ** 0))
            plt.xlim((5, 20000))
            plt.grid()
            plt.title("Fractional Abundance of " + self.atom + " in $N_{e}$  = " + "{:.2e}".format(self.xnew[i_Ne]) + " $cm^{-3}$", fontsize=16)
            plt.xlabel("$T_{e}$ [eV]", fontsize=16)
            plt.ylabel("FA", fontsize=16)

        plt.show()

    def worker(self,ion , product_all, sum_all):
        fun = self.get_Fractional_Abundance(ion = ion, product_all = np.array(product_all), sum_all = np.array(sum_all))
        return fun

    def calculate(self):
        ion_list = list(range(self.Z + 1))
        pool = ThreadPoolExecutor(self.Z + 1)
        pp = [pool.submit(self.worker,ion = ion, product_all = self.product_all, sum_all = self.sum_all) for ion in range(len(ion_list))]
        for p in pp:
            self.FA_arr.append(p.result())

    def create_dataset(self,output_filepath='.',filename='fractional_abundance.dat'):
        columns = ['T']
        for Z in range(self.Z):
            columns.append('Z{}'.format(Z+1))
        FA_output_df = pd.DataFrame(columns=columns)
        for i in range(self.Z):
            FA_output_df[columns[i+1]] = self.FA_arr[i+1][:,40]*100
        FA_output_df['T'] = self.ynew
        FA_output_df.to_csv(os.path.join(output_filepath,filename),sep=" ",index=False)


if __name__ == '__main__':
    path_to_data = r"C:\Users\marci\Desktop\Projekt NCN\Zadania\1.Styczeń\Fractional_Abundance\data\unresolved"
    t1 = time.time()
    FA = FractionalAbundance(atom='Xe',concurrent=True, path_to_data = path_to_data)
    t2 = time.time()
    print(t2-t1)
    t1 = time.time()
    FA = FractionalAbundance(atom='Xe', concurrent=True, path_to_data = path_to_data)
    t2 = time.time()
    print(t2 - t1)



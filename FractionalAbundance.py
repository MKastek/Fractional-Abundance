import numpy as np
import os
import matplotlib.pyplot as plt
from scipy import interpolate
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import time

plt.rcParams['figure.figsize'] = [12, 7]


class FractionalAbundance:

    def __init__(self, atom, SCD_file, ACD_file, connurent = False):
        self.atom = atom
        self.Z, self.num_of_Ne_axes, self.num_of_Te_axes, self.num_of_lines_to_read_with_axes, self.sum_of_axes = self.read_first_line_of_file(SCD_file)
        self.Te, self.Ne = self.read_axes(SCD_file)
        self.empty_matrix = np.empty([self.num_of_Te_axes, self.num_of_Ne_axes])
        self.num_of_lines_to_read_with_const_Te = int(np.ceil(self.num_of_Ne_axes / 8))
        self.start_line = int(np.ceil(self.sum_of_axes / 8)) + 3
        self.stop_line = self.start_line + self.num_of_lines_to_read_with_const_Te * self.num_of_Te_axes
        self.move = self.stop_line - self.start_line + 1
        self.FA_arr = []
        self.xnew, self.ynew = np.logspace(10, 15, num=100), np.logspace(1, 4, num=1000)

        if not connurent:
            self.FA_arr = [self.get_Fractional_Abundance(i, ACD_file, SCD_file) for i in range(self.Z + 1)]

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

        if type == 'SCD':
            CD = {str(i) + str(i + 1): self.read_data_from(filepath, self.start_line + i * self.move,
                                                           self.stop_line + i * self.move,
                                                           self.num_of_lines_to_read_with_const_Te,
                                                           self.empty_matrix.copy()) for i in range(self.Z)}
        elif type == 'ACD':
            CD = {str(i + 1) + str(i): self.read_data_from(filepath, self.start_line + i * self.move,
                                                           self.stop_line + i * self.move,
                                                           self.num_of_lines_to_read_with_const_Te,
                                                           self.empty_matrix.copy()) for i in range(self.Z)}
        return CD

    def get_Fractional_Abundance(self, ion, ACD_file, SCD_file):

        ACD, SCD = self.read_coefficients_matrices(ACD_file, type='ACD'), self.read_coefficients_matrices(SCD_file, type='SCD')

        K = [np.divide(10 ** SCD[str(i) + str(i + 1)], 10 ** ACD[str(i + 1) + str(i)]) for i in range(self.Z)]
        K.insert(0, np.ones_like(K[0]))

        product_all = np.cumprod(K, axis=0)
        FA = np.divide(product_all[ion], np.sum(product_all, axis=0))

        return FA

    def plot_FA_all(self, i_Ne=50):
        for i in range(self.Z + 1):
            x = self.ynew
            y = self.FA_arr[i][:, i_Ne]
            plt.plot(x, y, label="$" + self.atom + "^{" + str(i) + "+}$")
            plt.xscale("log")
            plt.yscale("log")
            plt.ylim((10 ** -3, 10 ** 0))
            plt.xlim((10 ** 1, 10 ** 4))
            plt.grid()
            plt.title("Fractional Abundance of " + self.atom + " in $N_{e}$  = " + "{:.2e}".format(self.xnew[i_Ne]) + " $cm^{-3}$", fontsize=16)
            plt.xlabel("$T_{e}$ [eV]", fontsize=16)
            plt.ylabel("FA", fontsize=16)

        plt.show()

    def worker(self, ion, ACD_file, SCD_file):
        fun = self.get_Fractional_Abundance(ion, ACD_file, SCD_file)
        return fun

    def calculate(self, ACD_file, SCD_file):
        ion_list = list(range(self.Z + 1))
        pool = ProcessPoolExecutor(self.Z + 1)
        pp = [pool.submit(self.worker, ion, ACD_file, SCD_file) for ion in range(len(ion_list))]
        for p in pp:
            self.FA_arr.append(p.result())


t1 = time.time()
if __name__ == '__main__':
    FA = FractionalAbundance(atom='Xe',SCD_file='scd89_xe.dat',ACD_file='acd89_xe.dat',connurent=True)
    FA.calculate(ACD_file='acd89_xe.dat', SCD_file='scd89_xe.dat')
    t2 = time.time()
    print(f"Execution time: {t2 - t1} s")
    FA.plot_FA_all(i_Ne=40)
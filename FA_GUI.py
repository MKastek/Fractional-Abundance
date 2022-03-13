from tkinter import *
import numpy as np
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
# Get the newest widget themes from Tk 8.5
from tkinter import ttk
from tkinter import filedialog as fd
# Create the main window that holds all the widgets
from matplotlib.figure import Figure
import pandas as pd
from tkinter.messagebox import showerror
import tkinter as tk
from FractionalAbundance import FractionalAbundance
import os
import threading


class FA_GUI(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title("Spectral data")
        self.geometry("1000x750")

        self.top_frame = ttk.Frame(self, padding="10 10 10 10")
        self.top_frame.pack(side=TOP)

        self.ACD_button = Button(self.top_frame, text="ACD", command=self.open_ACD_file, width=10, height=2)
        self.SCD_button = Button(self.top_frame, text="SCD", command=self.open_SCD_file, width=10, height=2)
        self.plot_button = Button(self.top_frame, text="plot", command=threading.Thread(target=self.plot).start, width=10, height=2)

        self.ACD_button.pack(side=LEFT)
        self.SCD_button.pack(side=LEFT)
        self.plot_button.pack(side=LEFT)

        self.fig = Figure(figsize=(5, 4), dpi=100)
        self.subplot = self.fig.add_subplot(111)
        self.subplot.grid()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        self.toolbar = NavigationToolbar2Tk(self.canvas, self)
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    def open_ACD_file(self):
        self.ACD_file = os.path.basename(fd.askopenfilename(initialdir='.'))


    def open_SCD_file(self):
        self.SCD_file = os.path.basename(fd.askopenfilename(initialdir='.'))

    def plot(self):
        FA = FractionalAbundance(atom='Xe', SCD_file=self.SCD_file, ACD_file=self.ACD_file)
        self.fig.clf()
        self.subplot = self.fig.add_subplot(111)
        for i in range(FA.Z + 1):
            x = FA.ynew
            y = FA.FA_arr[i][:, 50]
            self.subplot.plot(x, y, label="$" + FA.atom + "^{" + str(i) + "+}$")
            self.subplot.grid()
            self.subplot.set_xscale("log")
            self.subplot.set_yscale("log")
            self.subplot.set_ylim((10 ** -3, 10 ** 0))
            self.subplot.set_xlim((10 ** 1, 10 ** 4))
        self.fig.canvas.draw_idle()

if __name__ == "__main__":
    Fa_GUI = FA_GUI()
    Fa_GUI.mainloop()
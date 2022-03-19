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
import FractionalAbundance
import os
import threading


class FA_GUI(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title("Spectral data")
        self.geometry("1000x750")

        self.top_frame = ttk.Frame(self, padding="10 10 10 10")
        self.top_frame.pack(side=TOP)

        self.status_variable = StringVar()
        self.status_variable.set("Status: Load the data")
        self.status_label = Label(self.top_frame,textvariable=self.status_variable)

        self.ACD_button = Button(self.top_frame, text="ACD", command=self.open_ACD_file, width=10, height=2)
        self.SCD_button = Button(self.top_frame, text="SCD", command=self.open_SCD_file, width=10, height=2)
        self.plot_button = Button(self.top_frame, text="plot", command=lambda: threading.Thread(target=self.plot).start(), width=10, height=2)

        self.status_label.pack(side=LEFT)
        self.ACD_button.pack(side=LEFT)
        self.SCD_button.pack(side=LEFT)
        self.plot_button.pack(side=LEFT)

        self.fig = Figure(figsize=(5, 4), dpi=100)
        self.plot_data()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        self.toolbar = NavigationToolbar2Tk(self.canvas, self)
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    def open_ACD_file(self):
        self.ACD_file = fd.askopenfilename(initialdir='.')
        print(self.ACD_file)


    def open_SCD_file(self):
        self.SCD_file = fd.askopenfilename(initialdir='.')
        print(self.SCD_file)

    def plot_data(self):
        self.subplot = self.fig.add_subplot(111)
        self.subplot.grid()
        self.subplot.set_title("Fractional Abundance")
        self.subplot.set_xlabel("$T_{e} [eV]$")
        self.subplot.set_ylabel("FA")
        self.subplot.set_xscale("log")
        self.subplot.set_yscale("log")
        self.subplot.set_ylim((10 ** -3, 10 ** 0))
        self.subplot.set_xlim((10 ** 1, 10 ** 4))

    def plot(self):
        self.status_variable.set("Status: Calculating...")
        FA = FractionalAbundance.FractionalAbundance(atom='Xe', SCD_file=self.SCD_file, ACD_file=self.ACD_file)
        self.fig.clf()
        self.plot_data()
        for i in range(FA.Z + 1):
            x = FA.ynew
            y = FA.FA_arr[i][:, 50]
            self.subplot.plot(x, y, label="$" + FA.atom + "^{" + str(i) + "+}$")
        self.fig.canvas.draw_idle()
        self.status_variable.set("Status: Load the data")


if __name__ == "__main__":
    FA_GUI = FA_GUI()
    FA_GUI.mainloop()
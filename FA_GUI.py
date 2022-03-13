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
from FractionalAbundance import FractionalAbundance
import os
import threading


def open_ACD_file():
    global ACD_file
    ACD_file = os.path.basename(fd.askopenfilename(initialdir='.'))


def open_SCD_file():
    global SCD_file
    SCD_file = os.path.basename(fd.askopenfilename(initialdir='.'))

def plot():
    FA = FractionalAbundance(atom='Xe', SCD_file=SCD_file, ACD_file=ACD_file)
    fig.clf()
    subplot = fig.add_subplot(111)
    for i in range(FA.Z + 1):
        x = FA.ynew
        y = FA.FA_arr[i][:, 50]
        subplot.plot(x, y, label="$" + FA.atom + "^{" + str(i) + "+}$")
        subplot.grid()
        subplot.set_xscale("log")
        subplot.set_yscale("log")
        subplot.set_ylim((10 ** -3, 10 ** 0))
        subplot.set_xlim((10 ** 1, 10 ** 4))

        fig.canvas.draw_idle()

root = Tk()
root.geometry('1000x750')
# Define the title for the window
root.title("Spectral data")


top_frame = ttk.Frame(root, padding="10 10 10 10")
top_frame.pack(side=TOP)

ACD_button = Button(top_frame, text="ACD",  command=open_ACD_file, width=10, height=2)
SCD_button = Button(top_frame, text="SCD", command=open_SCD_file,width=10, height=2)
plot_button = Button(top_frame, text="plot", command=threading.Thread(target=plot).start, width=10, height=2)

ACD_button.pack(side=LEFT)
SCD_button.pack(side=LEFT)
plot_button.pack(side=LEFT)




fig = Figure(figsize=(5, 4), dpi=100)
subplot = fig.add_subplot(111)
subplot.grid()
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.draw()
canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

toolbar = NavigationToolbar2Tk(canvas, root)
toolbar.update()
canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

root.mainloop()
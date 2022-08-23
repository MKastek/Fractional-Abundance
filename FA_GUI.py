from tkinter import *
import numpy as np
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
# Get the newest widget themes from Tk 8.5
from tkinter import ttk
from tkinter import filedialog as fd
# Create the main window that holds all the widgets
from matplotlib.figure import Figure
import pandas as pd
import tkinter as tk
from Fractional_Abundace import FractionalAbundance
import os
import threading


class FA_GUI(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title("Fractional Abundance")
        self.geometry("1000x750")

        self.top_frame = ttk.Frame(self, padding="10 10 10 10")
        self.top_frame.pack(side=TOP)

        self.status_variable = StringVar()
        self.status_variable.set("Status: Load the ACD, SCD files")
        self.status_label = Label(self.top_frame,textvariable=self.status_variable)

        self.ACD_button = Button(self.top_frame, text="ACD", command=self.open_ACD_file, width=10, height=2)
        self.SCD_button = Button(self.top_frame, text="SCD", command=self.open_SCD_file, width=10, height=2)
        self.plot_button = Button(self.top_frame, text="plot", command=lambda: threading.Thread(target=self.plot).start(), width=10, height=2)
        self.save_button = Button(self.top_frame, text="save", command=self.save,width=10, height=2)

        self.status_label.pack(side=LEFT)
        self.ACD_button.pack(side=LEFT)
        self.SCD_button.pack(side=LEFT)
        self.plot_button.pack(side=LEFT)
        self.save_button.pack(side=LEFT)

        self.slider_frame_Te = ttk.Frame(self.top_frame, padding="10 10 10 10")
        self.slider_frame_Te.pack(side=RIGHT)

        self.Te_label = Label(self.slider_frame_Te, text="Te [eV]")
        self.Te_label.pack(side=TOP)

        self.Te_scale = DoubleVar()
        self.Te_scale.set(4.0)
        self.slider_Te = Scale(self.slider_frame_Te, orient='horizontal', from_=1.1, to=4.0, resolution=0.1, variable=self.Te_scale,command=lambda x: self.replot())
        self.slider_Te.pack(side=TOP)

        self.fig = Figure(figsize=(5, 4), dpi=100)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_plot_hover)
        self.plot_data()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        self.toolbar = NavigationToolbar2Tk(self.canvas, self)
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    def open_ACD_file(self):
        self.ACD_file = fd.askopenfilename(initialdir='.')
        self.element = os.path.basename(self.ACD_file)[-6:-4]

    def open_SCD_file(self):
        self.SCD_file = fd.askopenfilename(initialdir='.')

    def plot_data(self, scale_Te=4.0):
        self.subplot = self.fig.add_subplot(111)
        self.annot = self.subplot.annotate("", xy=(300, 0.01), xytext=(0, 0), textcoords="offset points")
        self.annot.set_visible(False)
        self.subplot.grid()
        self.subplot.set_title("Fractional Abundance")
        self.subplot.set_xlabel("$T_{e}$ [eV]")
        self.subplot.set_ylabel("FA")
        self.subplot.set_xscale("log")
        self.subplot.set_yscale("log")
        self.subplot.set_ylim((10 ** -3, 10 ** 0))
        self.subplot.set_xlim((10 ** 1, 10 ** scale_Te))

    def plot(self, scale_Te=4.0):
        self.status_label.config(font=("Segoe UI",12, "bold"))
        self.status_variable.set("Status: Calculating...")
        self.FA = FractionalAbundance.FractionalAbundance(element=self.element)
        self.fig.clf()
        self.plot_data(scale_Te=scale_Te)
        for i in range(self.FA.Z + 1):
            x = self.FA.ynew
            y = self.FA.FA_arr[i][:, 50]
            self.subplot.plot(x, y, label="$" + self.FA.element + "^{" + str(i) + "+}$")
        self.fig.canvas.draw_idle()
        self.status_label.config(font=("Segoe UI", 10))
        self.status_variable.set("Status: Load the data")

    def replot(self):
        self.subplot.set_xlim((10 ** 1, 10 ** self.Te_scale.get()))
        self.fig.canvas.draw_idle()

    def save(self):
        df_FA = pd.DataFrame(columns=[i for i in range(self.FA.Z)])
        for i in range(self.FA.Z + 1):
            df_FA[i] = self.FA.FA_arr[i][:, 50]
        df_FA['T[eV]'] = np.logspace(1, 4, num=1000)
        file = fd.asksaveasfilename(initialdir='.',defaultextension=".csv")
        df_FA.to_csv(os.path.basename(file))

    def on_plot_hover(self,event):
        vis = self.annot.get_visible()
        for curve,num in zip(self.subplot.get_lines(), range(len(self.subplot.get_lines()))):
            if curve.contains(event)[0]:
                #print(num)
                #print(event.xdata,event.ydata)
                self.annot.xy = (event.xdata, event.ydata)
                self.annot.set_text(self.element.capitalize()+" "+str(num)+"+")
                self.annot.set_visible(True)
                self.fig.canvas.draw_idle()
            else:
                if vis:
                    self.annot.set_visible(False)
                    self.fig.canvas.draw_idle()


if __name__ == "__main__":
    FA_GUI = FA_GUI()
    FA_GUI.mainloop()
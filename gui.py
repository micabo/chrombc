# MB
"""
Evaluation of chromatographic data.
Write a mock-up of the ChemStation Integrator (as can be inferred from the manual).
All times in seconds.
"""

from cutil import *
from integrator import CSIntegrator

import sys

import tkinter as tk
from tkinter import ttk
import tkinter.filedialog as tkfile

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases import key_press_handler


class ChromGUI(tk.Frame):


    def __init__(self, master = None, plotting_fn = lambda x: None):
        tk.Frame.__init__(self, master)
        self.master.wm_title("ChroMBC")
        self.frm_plot = tk.Frame(self)
        self.frm_btn = tk.Frame(self)

        self.fig = Figure(figsize=(6,4), dpi=200)
        #self.fig = Figure()
        self.ax = self.fig.add_subplot(111)

        plotting_fn(self.ax)

        self.canvas = FigureCanvasTkAgg(self.fig, master = self.frm_plot)
        self.canvas.draw()

        self.toolbar = NavigationToolbar2Tk(self.canvas, self.frm_plot)
        self.toolbar.update()

        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.canvas.mpl_connect("key_press_event", key_press_handler)

        self.integral = [[],[]]

        self.canvas.mpl_connect("button_press_event", self.draw_integral)
        self.canvas.mpl_connect("button_release_event", self.draw_integral)

# =============================================================================
#         self.canvas.mpl_connect("button_press_event", self.draw_blobb('go'))
#         self.canvas.mpl_connect("button_release_event", self.draw_blobb('ro'))
# =============================================================================

        self.btn_int_tgl = tk.Button(master=self.frm_btn, text="Int Tgl", command=self.toggle_int)
        self.btn_int_tgl.pack(side=tk.LEFT)

        self.btn_quit = tk.Button(master = self.frm_btn, text="Quit", command = self._quit)
        self.btn_quit.pack(side=tk.LEFT)

        self.frm_btn.pack(side=tk.TOP)
        self.frm_plot.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=1)

        self.integrate = False



    def _quit(self):
        self.master.quit()     # stops mainloop
        self.master.destroy()  # this is necessary on Windows to prevent
                        # Fatal Python Error: PyEval_RestoreThread: NULL tstate


    def draw_integral(self, e):
        if not self.integrate:
            return
        if len(self.integral[0]) == 0:
            self.integral[0].append(e.xdata)
            self.integral[1].append(e.ydata)
        elif len(self.integral[0]) == 1:
            self.integral[0].append(e.xdata)
            self.integral[1].append(e.ydata)
            self.ax.plot(self.integral[0], self.integral[1], 'k', linewidth=.5)
            self.canvas.draw()
        else:
            self.integral = [[],[]]
            self.draw_integral(e)


    def draw_blobb(self, style):
        def blobb(e):
            if self.integrate:
                self.ax.plot(e.xdata, e.ydata, style)
                self.canvas.draw()
        return blobb


    def toggle_int(self):
        self.integrate = not self.integrate



class ChromGUI2(tk.Frame):
    # do not use matplotlib

    def __init__(self, master = None, plotting_fn = lambda x: None):
        tk.Frame.__init__(self, master)
        self.master.wm_title("ChroMBC")

        self.frm_plot = tk.Frame(self)
        self.frm_btn = tk.Frame(self)

        self.menubar = tk.Menu(master = self)
        self.filemenu = tk.Menu(self.menubar)
        self.filemenu.add_command(label="Load", command = self.load)
        self.filemenu.add_command(label="Save As", command = self.saveas)
        self.menubar.add_cascade(label="Data", menu = self.filemenu)
        self.master.config(menu = self.menubar)

        self.canvas = tk.Canvas(self.frm_plot, width=400, height=400, bg="white")
        self.canvas.create_line(0, 0, 50, 50)
        self.canvas.bind("<Button-1>", self.click)
        self.canvas.bind("<Configure>", self.resize)
        self.canvas.pack(fill=tk.BOTH, expand=1)

        self.btn_int_tgl = ttk.Button(master=self.frm_btn, text="Int Tgl", command=self.toggle_int)
        self.btn_int_tgl.pack(side=tk.LEFT)

        self.btn_quit = ttk.Button(master=self.frm_btn, text="Quit", command=self._quit)
        self.btn_quit.pack(side=tk.LEFT)

        self.frm_btn.pack(side=tk.TOP)
        self.frm_plot.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=1)

        self.integrate = False

    def load(self):
        self.file = tkfile.askopenfile()
        if self.file:
            print("opening: ", self.file)

    def saveas(self):
        # problem: does not save immediately
        self.file = tkfile.asksaveasfile()
        if self.file:
            self.file.write("saved")


    def click(self, event):
        c = event.widget
        x = c.canvasx(event.x)
        y = c.canvasy(event.y)
        self.canvas.create_line(0, 0, x, y)


    def resize(self, event):
        w = event.widget
        print(w, event.width, event.height)


    def _quit(self):
        self.master.destroy()
        self.master.quit()


    def toggle_int(self):
        self.integrate = not self.integrate


if __name__ == "__main__":
    c = ChromData("./SST.txt")
    csint = CSIntegrator(c)
    csint.find_peaks()

    # start GUI
    root = tk.Tk()
    ChromGUI(root, csint.plot_on).pack()
    root.mainloop()

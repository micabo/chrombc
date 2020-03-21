# -*- coding: utf-8 -*-
"""
Evaluation of chromatographic data.
GUI
"""

from cutil import ChromData
from integrator import CSIntegrator, CSIntegrator

import tkinter as tk
from tkinter import ttk
import tkinter.filedialog as tkfile

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases import key_press_handler



class RibbonFrame(ttk.Frame):
    "Frame containing several buttons etc. on top of screen == 'Ribbon'"
    def __init__(self, master=None):
        ttk.Frame.__init__(self, master)

        #self.ent_1_str = tk.StringVar()

        self.btn_1 = ttk.Button(master=self, text="Toggle Integration", command=master._toggle_int)
        self.btn_2 = ttk.Button(master=self, text="Run Integration", command=self.master.run_integration)
        self.btn_3 = ttk.Button(master=self, text="Quit", command=master._quit)
        #self.ent_1 = tk.Entry(master=self, textvariable=self.ent_1_str)

        #self.ent_1.pack(side=tk.LEFT)
        self.btn_2.pack(side=tk.LEFT)
        self.btn_1.pack(side=tk.LEFT)
        self.btn_3.pack(side=tk.LEFT)

        #self.ent_1.focus_set()
        #self.ent_1_str.set("default")
        #self.ent_1_str.get()



class SidePaneFrame(ttk.Frame):
    "SidePane contains a listbox populated with data files"
    def __init__(self, master=None, datafiles=None):
        ttk.Frame.__init__(self, master)

        self.lib_1 = tk.Listbox(self, listvariable=datafiles, height=20)
        self.scroll = ttk.Scrollbar(self, orient=tk.VERTICAL, command=self.lib_1.yview)
        self.lib_1['yscrollcommand'] = self.scroll.set
        self.lib_1.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
        self.scroll.pack(side=tk.LEFT, fill=tk.Y, expand=1)


    def add_entry(self, index, entry):
        self.lib_1.insert(index, entry)



class PlotFrame(ttk.Frame):
    "Handles the display of data and interactive integration"
    def __init__(self, master=None, use_mpl=False):
        ttk.Frame.__init__(self, master)
        if use_mpl:
            self._build_matplotlib()
        else:
            self._build_native()

        self.integral = [[],[]]


    def _build_native(self):
        self.canvas = tk.Canvas(self, width=1000, height=600, bg="white")
        self.canvas.create_line(0, 0, 50, 50)
        self.canvas.bind("<Button-1>", self._click)
        self.canvas.bind("<Configure>", self._resize)
        self.canvas.pack(fill=tk.BOTH, expand=1)


    def _build_matplotlib(self):
        self.fig = Figure(figsize=(5,3), dpi=200)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self)
        self.toolbar.update()

        self.canvas.mpl_connect("key_press_event", key_press_handler)
        self.canvas.mpl_connect("button_press_event", self.draw_integral)
        self.canvas.mpl_connect("button_release_event", self.draw_integral)

        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)


    def _click(self, event):
        c = event.widget
        x = c.canvasx(event.x)
        y = c.canvasy(event.y)
        self.canvas.create_line(0, 0, x, y)


    def _resize(self, event):
        w = event.widget
        print(w, event.width, event.height)


    def plot_curve(self, data):
        self.ax.clear()  # clears everything - also the zoom level...
        self.ax.plot(data.x, data.y)
        self.canvas.draw()


    def plot_peaks(self, peaks):
        for peak in peaks:
            self.ax.plot(
                [peak.start.x, peak.end.x],
                [peak.start.y, peak.end.y],
                'k', linewidth = 0.5)
        self.canvas.draw()


    def draw_integral(self, e):
        if not self.master.integrate:
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



class TestGUI(ttk.Frame):
    "Main Application"
    def __init__(self, master=None):
        ttk.Frame.__init__(self, master)

        #self.datafiles = tk.StringVar(value=["file {}".format(i) for i in range(100)])
        self.datafiles = tk.StringVar(value=[])
        self.data = CSIntegrator()
        self.current_data = None

        self._setup_menubar()
        self._setup_panes()

        self.bind_all("<<ListboxSelect>>", self.set_current_data)
        self.bind_all("<Escape>", lambda x: print("esc", x))

        self.integrate = False


    def set_current_data(self, event):
        # get index of currently selected data
        self.current_data = event.widget.curselection()[0]
        self.plotdisplay.plot_curve(self.data.cdata[self.current_data])
        if self.data.cdata[self.current_data].peaks:
            self.plotdisplay.plot_peaks(
                self.data.cdata[self.current_data].peaks)


    def _setup_menubar(self):
        self.master.wm_title("ChroMBC")
        self.menubar = tk.Menu(master=self)
        self.filemenu = tk.Menu(self.menubar)
        self.filemenu.add_command(label="Load", command=self._load)
        self.filemenu.add_command(label="Save As", command=self._saveas)
        self.menubar.add_cascade(label="File", menu=self.filemenu)
        self.master.config(menu=self.menubar)


    def _setup_panes(self):
        self.ribbon = RibbonFrame(self)
        self.sidepane = SidePaneFrame(self, datafiles=self.datafiles)
        self.plotdisplay = PlotFrame(self, use_mpl=True)

        self.ribbon.pack(side=tk.TOP)
        self.sidepane.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
        self.plotdisplay.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)


    def _load(self):
        file = tkfile.askopenfile()
        if file:
            self.current_data = len(self.data.cdata)
            self.data.add_chromatogramm(ChromData(file.name))
            self.plotdisplay.plot_curve(self.data.cdata[self.current_data])
            self.sidepane.add_entry(self.current_data, file.name.split("/")[-1])


    def _saveas(self):
        self.file = tkfile.asksaveasfile()
        if self.file:
            self.file.write("saved")


    def _quit(self):
        self.master.destroy()
        self.master.quit()


    def _toggle_int(self):
        self.integrate = not self.integrate


    def run_integration(self):
        self.data.find_peaks(self.current_data)
        self.plotdisplay.plot_peaks(self.data.cdata[self.current_data].peaks)


if __name__ == "__main__":
    root = tk.Tk()
    TestGUI(root).pack(fill=tk.BOTH, expand=1)
    root.mainloop()

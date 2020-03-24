# -*- coding: utf-8 -*-
"""
Evaluation of chromatographic data.
GUI
"""
import os

from cutil import ChromData
from integrator import CSIntegrator, MBIntegrator

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases import key_press_handler

import itertools as it


#-----------------------------------------------------------------------------
class RibbonFrame(ttk.Frame):
    "Frame containing several buttons etc. on top of screen == 'Ribbon'"
    def __init__(self, master=None, target=None):
        ttk.Frame.__init__(self, master)
        self.target = target

        self.smooth_val = tk.StringVar()
        self.smooth_val.set("0")

        self.ent_smooth = tk.Entry(master=self, textvariable=self.smooth_val)
        self.btn_smooth = ttk.Button(master=self, text="Smooth",
                                     command=self.target.apply_smoothing)
        self.btn_baseline = ttk.Button(master=self, text="Baseline",
                                       command=self.target.draw_baseline)
        self.btn_int = ttk.Button(master=self, text="Run Integration",
                                  command=self.target.run_integration)
        self.btn_tgl = ttk.Button(master=self, text="Toggle Integration",
                                  command=self.target.toggle_int)

        self.btn_quit = ttk.Button(master=self, text="Quit", command=self.target._quit)

        self.ent_smooth.pack(side=tk.LEFT)
        self.btn_smooth.pack(side=tk.LEFT)
        self.btn_baseline.pack(side=tk.LEFT)
        self.btn_int.pack(side=tk.LEFT)
        self.btn_tgl.pack(side=tk.LEFT)
        self.btn_quit.pack(side=tk.LEFT)

        self.ent_smooth.focus_set()


    def get_sval(self):
        return int(self.smooth_val.get())


#-----------------------------------------------------------------------------
class DataChoose(ttk.Frame):
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


#-----------------------------------------------------------------------------
class PlotFrame(ttk.Frame):
    "Handles the display of data and interactive integration"
    def __init__(self, master=None, target=None):
        ttk.Frame.__init__(self, master)
        self.target = target
        self._build_matplotlib()

        self.integral = [[],[]]
        self.color = it.cycle(["black", "red"])


    def _build_matplotlib(self):
        self.fig = Figure(figsize=(5,3), dpi=200)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.tk_canvas = self.canvas.get_tk_widget()
        self.toolbar = NavigationToolbar2Tk(self.canvas, self)
        self.toolbar.update()

        self.tk_canvas = self.canvas.get_tk_widget()

        self.canvas.mpl_connect("key_press_event", key_press_handler)
        self.canvas.mpl_connect("button_press_event", self.draw_integral)
        self.canvas.mpl_connect("button_release_event", self.draw_integral)
# =============================================================================
#         self.tk_canvas.bind("<Configure>", self._resize)
#         self.tk_canvas.bind("<Button-1>", self._click)
# =============================================================================
        self.canvas.draw()
        self.tk_canvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)


    def _click(self, event):
        c = event.widget
        x = c.canvasx(event.x)
        y = c.canvasy(event.y)
        c.create_line(0, 0, x, y)


    def draw_integral(self, e):
        if not self.target.integrate:
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


    def _resize(self, event):
        w = event.widget
        print(w, event.width, event.height)


    def plot_curve(self, x, y, clear=True, **pparam):
        if clear:
            self.ax.clear()
        self.ax.plot(x, y, **pparam)
        self.canvas.draw()


    def plot_peaks(self, peaks, clear=False):
        if clear:
            self.ax.clear()
        color = next(self.color)
        for peak in peaks:
            self.ax.plot(
                [peak.start.x, peak.end.x],
                [peak.start.y, peak.end.y],
                linewidth = 0.5,
                color=color)
        self.canvas.draw()


#-----------------------------------------------------------------------------
class ResultFrame(ttk.Frame):
    "Display the results"
    def __init__(self, master=None):
        ttk.Frame.__init__(self, master)
        self.result_scroll = tk.Scrollbar(master=self)
        self.result_txt = tk.Text(self, yscrollcommand=self.result_scroll.set, height=10)
        self.result_scroll.config(command=self.result_txt.yview)
        self.result_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.result_txt.pack(side=tk.TOP, fill=tk.X)


    def print_results(self, peak_table):
        self.result_txt.delete("1.0", tk.END)

        l = ["{No:>3}|{RT:>10.2f}|{Area:>15.3f}|{Area%:>10.2f}".format(**peak) for peak in peak_table]
        s = "\n".join(l)
        header = "{:3}|{:10}|{:15}|{:10}\n".format("No", "RT", "Area", "Area%")
        header += "{:3}|{:10}|{:15}|{:10}\n".format("-" * 3, "-" * 10, "-" * 15, "-" * 10)
        self.result_txt.insert("1.0", header + s)


#-----------------------------------------------------------------------------
class IntegratorFrame(ttk.Frame):
    "Interaction with integrator parameters"
    def __init__(self, master=None, target=None):
        ttk.Frame.__init__(self, master)
        self.target = target
        # gui interface to the parameters of the integrator
        self.frames = []
        self.entries = []
        self.buttons = []
        for label in ["threshold", "min_height", "min_width"]:
            self.frames.append(ttk.Frame(self))
            self.entries.append(tk.Entry(self.frames[-1]))
            self.buttons.append(ttk.Button(
                self.frames[-1],
                text=label,
                command=self.change_param(self.entries[-1], label)))

            self.entries[-1].pack(side=tk.LEFT)
            self.buttons[-1].pack(side=tk.LEFT)
            self.frames[-1].pack(side=tk.TOP)


    def change_param(self, entry, name):
        def _change():
            self.target.change_parameter(name, float(entry.get()))
        return _change


#-----------------------------------------------------------------------------
class TestGUI(ttk.Frame):
    "Main Application"
    def __init__(self, master=None):
        ttk.Frame.__init__(self, master)

        #self.datafiles = tk.StringVar(value=["file {}".format(i) for i in range(100)])
        self.datafiles = tk.StringVar(value=[])
        self.data = MBIntegrator()
        self.current_data = 0

        self._setup_menubar()
        self._setup_panes()

        self.bind_all("<<ListboxSelect>>", self.set_current_data)
        self.bind_all("<Escape>", lambda x: print("esc", x))

        self.integrate = False


    def _setup_menubar(self):
        self.master.wm_title("ChroMBC")
        self.menubar = tk.Menu(master=self)
        self.filemenu = tk.Menu(self.menubar)
        self.filemenu.add_command(label="Load", command=self._load)
        self.filemenu.add_command(label="Save As", command=self._saveas)
        self.menubar.add_cascade(label="File", menu=self.filemenu)
        self.master.config(menu=self.menubar)


    def _setup_panes(self):
        self.frame1 = ttk.Frame(self)
        self.ribbon = RibbonFrame(self.frame1, target=self)
        self.ribbon.pack(side=tk.TOP)

        self.frame2 = ttk.Frame(self)
        self.data_choose = DataChoose(self.frame2, datafiles=self.datafiles)
        self.data_choose.pack(fill=tk.BOTH, expand=1)

        self.frame3 = ttk.Frame(self)
        self.plotdisplay = PlotFrame(self.frame3, target=self)
        self.plotdisplay.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        self.tabbed_view = ttk.Notebook(self.frame3)
        self.result_view = ResultFrame(self.tabbed_view)
        self.tabbed_view.add(self.result_view, text="Results")
        self.int_view = IntegratorFrame(self.tabbed_view, self)
        self.tabbed_view.add(self.int_view, text="Integrator")
        self.tabbed_view.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        self.frame1.pack(side=tk.TOP)
        self.frame2.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
        self.frame3.pack(side=tk.LEFT, fill=tk.X, expand=1)


    def _load(self):
        file = filedialog.askopenfile(initialdir=os.getcwd() + "\\data")
        if file:
            self.current_data = len(self.data)
            self.data.add_chromatogramm(ChromData(file.name))
            self.plotdisplay.plot_curve(*self.data[self.current_data])
            self.data_choose.add_entry(self.current_data, file.name.split("/")[-1])


    def _saveas(self):
        self.file = filedialog.asksaveasfile()
        if self.file:
            self.file.write("saved")


    def _quit(self):
        self.master.destroy()
        self.master.quit()


    def toggle_int(self):
        self.integrate = not self.integrate


    def set_current_data(self, event):
        # get index of currently selected data
        try:
            self.current_data = event.widget.curselection()[0]
        except IndexError:
            # no data loaded yet, tuple returned from curselection() is empty
            return

        self.plotdisplay.plot_curve(*self.data[self.current_data])
        if self.data[self.current_data].peaks:
            self.plotdisplay.plot_peaks(
                *self.data[self.current_data].peaks)


    def change_parameter(self, name, value):
        self.data.settings[name] = value


    def apply_smoothing(self):
        width = self.ribbon.get_sval()
        if width == 0:
            return
        self.data[self.current_data].smooth(width)
        self.plotdisplay.plot_curve(*self.data[self.current_data])


    def draw_baseline(self):
        baseline = self.data.find_baseline(self.current_data)
        self.plotdisplay.plot_curve(*baseline, clear=False, linewidth=0.5)


    def run_integration(self):
        self.data.find_peaks(self.current_data)
        self.plotdisplay.plot_peaks(self.data[self.current_data].peaks)
        self.display_results()


    def display_results(self):
        peak_table = self.data[self.current_data].get_peak_table()
        self.result_view.print_results(peak_table)


#-----------------------------------------------------------------------------
if __name__ == "__main__":
    root = tk.Tk()
    root.geometry("1400x800+0+0")
    TestGUI(root).pack(fill=tk.BOTH, expand=1)
    root.mainloop()

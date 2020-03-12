# MB
"""
Evaluation of chromatographic data.
Write a mock-up of the ChemStation Integrator (as can be inferred from the manual).
All times in seconds.
"""

import numpy as np
from functools import partial
import pandas as pd
from bisect import bisect
import scipy.signal as sig
from scipy.integrate import trapz, simps, quad
from scipy.optimize import curve_fit
from scipy.special import erfc
import matplotlib.pyplot as plt
from netCDF4 import Dataset

import tkinter
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

#-----------------------------------------------------------------------------
# constants

SQRT_PI = np.sqrt(np.pi)
SQRT_1_2 = np.sqrt(2)/2

#-----------------------------------------------------------------------------
# global functions

def gaussian(x, h, mu, sigma):
    """Return value of gaussian function at point x
    h: 'height'
    mu: mean
    sigma: std dev
    """
    z = (x - mu)/sigma
    return h * np.exp(-z*z/2)


def emg(x, h, mu, sigma, tau):
    """Return value of exponentially-modified gaussian at point x
    h: 'height'
    mu: mean
    sigma: std dev
    tau: decay
    """
    dx = x - mu
    s_t = sigma/tau
    return (h * s_t * SQRT_PI * SQRT_1_2 * np.exp(s_t**2/2 - dx/tau) *
            erfc((s_t - dx/sigma) * SQRT_1_2))


def fit(x, y, start_time, stop_time, f, initial_guess):
    start_index = bisect(x, start_time)
    stop_index = bisect(x, stop_time)
    params = curve_fit(f, x[start_index:stop_index], y[start_index:stop_index], initial_guess)
    return params[0]

#-----------------------------------------------------------------------------
# classes

class ChromData:
    """Interface to chromatographic data stored as either
    an AIA/NETCDF-File or .dat/.txt file, where the data is given as x-y data.
    An equal spacing in time is assumed
    """

    @staticmethod
    def construct_time(n, dt):
        "Returns a generator which construct a sequence with a given length n and the step dt"
        result = 0.0
        while n:
            yield result
            result += dt
            n -= 1


    @staticmethod
    def bytearray2string(x):
        return "".join([element.decode('ascii') for element in x])

    def __init__(self, path):
        if ".cdf" == path[-4:]:
            self._build_from_cdf(path)
        else:
            self._build_from_xy(path)


    def _build_from_cdf(self, path):
        d = Dataset(path, "r", format="NETCDF3_CLASSIC")
        self.y = np.array(d.variables["ordinate_values"])
        self.dt = float(d.variables["actual_sampling_interval"].getValue())
        self.x = np.fromiter(ChromData.construct_time(len(self.y), self.dt), np.float32)
        self.peak_names = [ChromData.bytearray2string(name.data) for name in d.variables["peak_name"]]
        self.peak_rts = [float(rt) for rt in d.variables["peak_retention_time"]]
        self.baseline = {"start_time": [float(x) for x in d.variables["baseline_start_time"]],
                         "stop_time": [float(x) for x in d.variables["baseline_stop_time"]],
                         "start_value": [float(x) for x in d.variables["baseline_start_value"]],
                         "stop_value": [float(x) for x in d.variables["baseline_stop_value"]]}
        self.y_unit = str(d.detector_unit)
        self.x_unit = str(d.retention_unit)
        self.raw = False
        d.close()


    def _build_from_xy(self, path):
        self.x, self.y = np.loadtxt(path, unpack=True)
        self.dt = self.x[1] - self.x[0]
        self.raw = True  # no additional information (baseline etc.) given


    def _write_xy(self, path):
        np.savetxt(path, np.array((self.x, self.y)).T)



class MBIntegrator:

    win_width = 21 # window_width for smoothing

    # point_types - a dictionary of possible point types
    point_types = {
            'undefined': 0,
            'baseline': 1,
            'peak-start': 2,
            'inflection-1': 3,
            'apex': 4,
            'inflection-2': 5,
            'peak-end': 6,
            'valley': 7,
            'skim-start': 8,
            'skim-end': 9,
            'shoulder-start': 10,
            'shoulder-end': 11
            }


    def __init__(self, cdata, **parameters):
        self.dt = cdata.dt
        self.x = np.array(cdata.x)
        self.y = np.array(cdata.y)
        self.N = len(self.x)

        # get the first and second derivative with a rolling average smoothing
        self.dy = pd.Series(np.gradient(self.y, self.dt)).rolling(
            window=MBIntegrator.win_width, center=True).mean()
        self.ddy = pd.Series(np.gradient(self.dy, self.dt)).rolling(
            window=MBIntegrator.win_width, center=True).mean()

        # point_type: an array containing a number for each point
        # this number defines the type of the point
        self.point_type = np.full_like(self.x, 0, dtype=np.int8)

        default_values = {
                "threshold": 0.1,
                "min_width": 0.5,
                "min_height": 0.05
                }
        self.settings = {**default_values, **parameters}


    def find_peaks(self):
        N_start = int((MBIntegrator.win_width + 1)/2)
        N_end = self.N - N_start
        i = N_start
        inflection_1_found = False
        inflection_2_found = False
        critical_points = []
        while i < N_end:
            #print(i)
            possible_peak = []
            if self.dy[i]*self.dy[i-1] < 0 and self.ddy[i] < 0:
                self.point_type[i-1] = self.point_types["apex"]
                possible_peak.append([self.x[i-1], self.y[i-1]])
                j = i - 2
                # go to the left and find peak start
                while j > N_start:
                    #print(j)
                    if self.ddy[j-1]*self.ddy[j] < 0: # found inflection
                        self.point_type[j-1] = self.point_types["inflection-1"]
                        inflection_1_found = True
                    if self.dy[j-1] < self.settings["threshold"] and inflection_1_found:
                        self.point_type[j] = self.point_types["peak-start"]
                        possible_peak.append([self.x[j], self.y[j]])
                        break
                    j -= 1
                # go to the right and find peak end
                k = i
                while k < N_end:
                    #print(k)
                    if self.ddy[k+1]*self.ddy[k] < 0: # found inflection
                        self.point_type[k+1] = self.point_types["inflection-2"]
                        inflection_2_found = True
                    if self.dy[k+1] > -self.settings["threshold"] and inflection_2_found:
                        self.point_type[k] = self.point_types["peak-end"]
                        possible_peak.append([self.x[k], self.y[k]])
                        break
                    k += 1
                # now go and find the next peaks start, if close enough, merge end and start and use the next peak end for height determination

                # reject peaks which are too small
                if len(possible_peak) == 3:
                    time_list = sorted(possible_peak, key=lambda x: x[0])
                    height_list = sorted(possible_peak, key=lambda x: x[1])
                    width = time_list[2][0] - time_list[0][0]
                    height = height_list[2][1] - height_list[0][1]
                    if width > self.settings["min_width"] and height > self.settings["min_height"]:
                        for point in possible_peak:
                            critical_points.append(point)

                # TODO: when start and end have been found, continue searching for peaks after the end of the current peak
                inflection_1_found = False
                inflection_2_found = False
                i = k + 1
            else:
                i += 1

        # last order of business - merge close lying end/start points of adjacent peaks
        # now start with fitting peaks!
        return list(zip(*critical_points))

    ## TODO:
    ## - Construct baseline from all peak-start to peak-end


class CSIntegrator:
    "Mock-up of the ChemStation Integrator."
    sampling_ratio = 15


    def __init__(self, cdata, **parameters):
        self.peak_start_indices = []
        self.dt = cdata.dt
        self.x = np.array(cdata.x)
        self.y = np.array(cdata.y)

        # get the first and second derivative with a rolling average smoothing
        self.dy = pd.Series(np.gradient(self.y, self.dt)).rolling(window=20).mean()
        self.ddy = pd.Series(np.gradient(self.dy, self.dt)).rolling(window=20).mean()

        # define the default values for the initial events
        self.events = {"slope_sensitivity": 1.0,
                       "peak_width": 3,
                       "area_reject": 0.0,
                       "area%_reject": 0.0,
                       "height_reject": 0.0}

        # overwrite these default values with the values given to the constructor
        for key, value in parameters.items():
            self.events[key] = value

        # timed events are given as tuples (event_name, parameter) if no parameter is required it is set to None
        self.timed_events = [
            [200, ("peak_width", 15)],
            [250, ("peak_width", 3)],
            [270, ("peak_width", 15)],
            [300, ("peak_width", 30)]]
        self.timed_events.sort(key=lambda x: x[0])


    def run(self, return_sampled_data=False):
        "Apply the integration algorithm. Returns the sampled data for inspection, analytical results are saved in the class instance."
        # go through the data from left to right and sample the data according to peak_width
        i = 0
        timed_event_index = 0
        step = 0
        samples = []
        while i + 2*step < len(self.x):
            # check if a timed event needs to be activated
            while timed_event_index < len(self.timed_events) and \
                  self.timed_events[timed_event_index][0] < self.x[int(i+step/2)]:
                self.events[self.timed_events[timed_event_index][1][0]] = self.timed_events[timed_event_index][1][1]
                timed_event_index += 1
            # define the step size according to the current peak_width
            step = int(self.events["peak_width"]/CSIntegrator.sampling_ratio/self.dt) or 1
            this_point = np.mean(self.y[i:i+step])
            next_point = np.mean(self.y[i+step:i+2*step])
            slope = (next_point - this_point)/(step*self.dt)
            if return_sampled_data:
                samples.append([self.x[int(i+step/2)], this_point])
            if slope >= self.events["slope_sensitivity"]:
                # at the moment this yields only possible peak starts, compare to the JS version for tricks etc.
                self.peak_start_indices.append(int(i + step/2))
            i += step
        return list(zip(*samples)) if return_sampled_data else None


if __name__ == "__main__":
    ##c = ChromData("./_007_008-0301_DAD1A.cdf")
    c = ChromData("./SST.txt")
    ##c_int = CSIntegrator(c, slope_sensitivity=0.5, peak_width=100)
    ##samples = c_int.run(True)
    ##indices = c_int.peak_start_indices
    ##params = fit(c.x, c.y, 290, 315, gaussian_bl, [50, 303, 1, 0, 0])
    ##fit_y = [gaussian_bl(xi, *params) for xi in c.x]

    d = MBIntegrator(c)
    points = d.find_peaks()

    root = tkinter.Tk()
    fig = plt.Figure(figsize=(6,4), dpi=300)
    ax = fig.add_subplot(111)
    ax.plot(d.x, d.y)
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
    tkinter.mainloop()
# =============================================================================
#     params_gauss = fit(d.x, d.y, 1040, 1080, gaussian, [50, 1060, 5])
#     params_emg = fit(d.x, d.y, 1040, 1080, emg, [50, 1060, 5, 1])
#     fit_gauss = [gaussian(xi, *params_gauss) for xi in d.x]
#     fit_emg = [emg(xi, *params_emg) for xi in d.x]
#
#     plt.plot(d.x, d.y)
#     plt.plot(d.x, fit_gauss)
#     plt.plot(d.x, fit_emg)
#     ##plt.plot(d.x, d.dy, linewidth=0.5)
#     ##plt.plot(d.x, d.ddy, linewidth=0.5)
#     plt.plot(*points, 'rx')
#     plt.show()
# =============================================================================

    #plt.ylabel(c.y_unit)
    #plt.xlabel(c.x_unit)
    ##plt.plot(c.x.take(indices), c.y.take(indices), 'gx')
    ##plt.plot(*samples)
    #plt.show()

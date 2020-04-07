# -*- coding: utf-8 -*-
"""Helper Constants/Functions/Classes
"""

from bisect import bisect_left

from scipy.integrate import trapz
from scipy.optimize import curve_fit
from scipy.special import erfc
from scipy.signal import find_peaks

from collections import namedtuple
from netCDF4 import Dataset

import numpy as np
import pandas as pd

#-----------------------------------------------------------------------------
# constants

SQRT_PI = np.sqrt(np.pi)
SQRT_1_2 = np.sqrt(2)/2


#-----------------------------------------------------------------------------
# named tuple datastructures

Point = namedtuple('Point', ['x', 'y'])
Peak = namedtuple('Peak', ['start', 'apex', 'end'])


#-----------------------------------------------------------------------------
# functions

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


def fit_peak(x, y, start, stop, f, initial_guess):
    "Fit x-y data with a function f within the specified range (start/stop)"
    start_i = bisect_left(x, start)
    stop_i = bisect_left(x, stop)
    params = curve_fit(f, x[start_i:stop_i], y[start_i:stop_i], initial_guess)
    return list(params[0])


def integrate_peak(peak, x, y):
    "Integrate a peak (given as named tuple -> need start/end)"
    # TODO: 0 is taken as the baseline -> needs to be adjusted / corrected!!!
    start_i = bisect_left(x, peak.start.x)
    stop_i = bisect_left(x, peak.end.x)
    return trapz(y[start_i:stop_i], x[start_i:stop_i])


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
        if path[-4:] == ".cdf":
            self._build_from_cdf(path)
        else:
            self._build_from_xy(path)

        self.peaks = []
        self.peak_table = []


    def __len__(self):
        return len(self.x)


    def __iter__(self):
        return iter((self.x, self.y))


    def _build_from_cdf(self, path):
        d = Dataset(path, "r", format="NETCDF3_CLASSIC")
        self.y = np.array(d.variables["ordinate_values"])
        self.dt = float(d.variables["actual_sampling_interval"].getValue())
        self.x = np.fromiter(ChromData.construct_time(len(self.y), self.dt), np.float64)
        self.peak_names = [ChromData.bytearray2string(name.data) for name in d.variables["peak_name"]]
        self.peak_rts = [float(rt) for rt in d.variables["peak_retention_time"]]
        self.baseline = {"start_time": [float(x) for x in d.variables["baseline_start_time"]],
                         "stop_time": [float(x) for x in d.variables["baseline_stop_time"]],
                         "start_value": [float(x) for x in d.variables["baseline_start_value"]],
                         "stop_value": [float(x) for x in d.variables["baseline_stop_value"]]}
        self.y_unit = str(d.detector_unit)
        self.x_unit = str(d.retention_unit)
        self.raw = False
        #TODO: the read data on peak names etc. so far is unused -> do sth with it
        d.close()


    def _build_from_xy(self, path):
        self.x, self.y = np.loadtxt(path, unpack=True)
        self.dt = self.x[1] - self.x[0]
        self.raw = True  # no additional information (baseline etc.) given


    def _write_xy(self, path):
        np.savetxt(path, np.array((self.x, self.y)).T)


    def clear_peaks(self):
        self.peaks = []
        self.peak_table = []


    def add_peak(self, peak):
        self.peaks.append(peak)


    def get_derivatives(self, width):
        """Returns first and second derivative of y with respect to x
        The derivatives are smoothed with a rolling average of width 'width'
        """
        self.dy = pd.Series(np.gradient(self.y, self.dt)).rolling(
            window=width, center=True).mean()
        self.ddy = pd.Series(np.gradient(self.dy, self.dt)).rolling(
            window=width, center=True).mean()
        return self.dy, self.ddy


    def smooth(self, width):
        "Smooth y data with a rolling average of width 'width'"
        self.y = pd.Series(self.y).rolling(window=width, center=True).mean()
        self.peaks = []
        self.peak_table = []
        return self.y


    def build_peak_table(self):
        cumulative_area = 0
        peak_no = 0
        for peak in self.peaks:
            peak_area = integrate_peak(peak, self.x, self.y)
            cumulative_area += peak_area
            peak_no += 1
            self.peak_table.append(dict(No=peak_no, RT=peak.apex.x, Area=peak_area))
        for peak in self.peak_table:
            peak['Area%'] = 100 * peak['Area'] / cumulative_area


    def get_peak_table(self):
        if len(self.peak_table) == 0:
            self.build_peak_table()
        return self.peak_table


#-----------------------------------------------------------------------------
if __name__ == "__main__":
    c = ChromData("./data/SST.txt")
    x, y = c

    # testing scipy.signal.find_peaks
    indices = find_peaks(y, prominence=(0.1, None))[0]
    y_max = np.take(y, indices, 0)
    x_max = np.take(x, indices, 0)
    import matplotlib.pyplot as plt
    plt.plot(x, y, '-', x_max, y_max, 'x')
    plt.show()

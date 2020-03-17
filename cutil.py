# -*- coding: utf-8 -*-
"""
Helper Constants/Functions/Classes
"""

import numpy as np
from bisect import bisect_left

from scipy.integrate import trapz
from scipy.optimize import curve_fit
from scipy.special import erfc

from collections import namedtuple

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


def fit(x, y, start, stop, f, initial_guess):
    "Fit x-y data with a function f within the specified range (start/stop)"
    start_i = bisect_left(x, start)
    stop_i = bisect_left(x, stop)
    params = curve_fit(f, x[start_i:stop_i], y[start_i:stop_i], initial_guess)
    return list(params[0])


def integrate_peak(peak, x, y):
    "Integrate a peak (given as named tuple -> need start/end)"
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
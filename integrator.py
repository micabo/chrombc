# MB
"""
Evaluation of chromatographic data.
Write a mock-up of the ChemStation Integrator (as can be inferred from the manual).
All times in seconds.
"""

from cutil import Point, Peak, gaussian, fit, integrate_peak
from bisect import bisect_left

import sys
import numpy as np
import pandas as pd


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
        'peak-rise': 8,
        'peak-fall': 9
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
            "min_width": 1,
            "min_height": 0.5
            }
        self.settings = {**default_values, **parameters}

        self.peaks = []
        self.peak_fits = []
        self.y_fit = None
        self.baseline = None
        self.y_bc = None


    def find_peaks(self):
        N_start = int((MBIntegrator.win_width + 1)/2)
        N_end = self.N - N_start
        i = N_start
        inflection_1_found = False
        inflection_2_found = False
        self.peaks = []
        while i < N_end:
            apex = start = end = False
            if self.dy[i]*self.dy[i-1] < 0 and self.ddy[i] < 0:
                self.point_type[i-1] = self.point_types["apex"]
                apex = Point(self.x[i-1], self.y[i-1])
                j = i - 2
                # go to the left and find peak start
                while j > N_start:
                    self.point_type[j] = self.point_types["peak-rise"]
                    if self.ddy[j-1]*self.ddy[j] < 0: # found inflection
                        self.point_type[j-1] = self.point_types["inflection-1"]
                        inflection_1_found = True
                    if self.dy[j-1] < self.settings["threshold"] and inflection_1_found:
                        self.point_type[j] = self.point_types["peak-start"]
                        start = Point(self.x[j], self.y[j])
                        break
                    j -= 1
                # go to the right and find peak end
                k = i
                while k < N_end:
                    self.point_type[k] = self.point_types["peak-fall"]
                    if self.ddy[k+1]*self.ddy[k] < 0: # found inflection
                        self.point_type[k+1] = self.point_types["inflection-2"]
                        inflection_2_found = True
                    if self.dy[k+1] > -self.settings["threshold"] and inflection_2_found:
                        self.point_type[k] = self.point_types["peak-end"]
                        end = Point(self.x[k], self.y[k])
                        break
                    k += 1
                # TODO: now go and find the next peaks start, if close enough, merge end and start and use the next peak end for height determination

                # reject peaks which are too small
                if apex and start and end:
                    width = end.x - start.x
                    height = apex.y - (start.y + end.y)/2
                    if width > self.settings["min_width"] and height > self.settings["min_height"]:
                        self.peaks.append(Peak(start, apex, end))

                # when start and end have been found:
                # continue searching for peaks after the end of the current
                # peak (k + 1)
                inflection_1_found = False
                inflection_2_found = False
                i = k + 1
            else:
                i += 1

        # TODO: last order of business - merge close lying end/start points of adjacent peaks
        return self.peaks


    def create_baseline(self):
        # create a crudely interpolated baseline
        # peaks which were not detected are taken as baseline
        self.baseline = np.copy(self.y)
        for p in self.peaks:
            start = bisect_left(self.x, p.start.x)
            end = bisect_left(self.x, p.end.x)
            slope = (self.y[end] - self.y[start])/(self.x[end] - self.x[start])
            i = start
            while i <= end:
                self.baseline[i] = self.y[start] + slope * (self.x[i] - self.x[start])
                i += 1
        return self.baseline


    def fit_peaks(self):
        "fit peaks with gaussian"
        self.create_baseline()
        self.y_bc = self.y - self.baseline
        for peak in self.peaks:
            height = peak.apex.y - (peak.start.y + peak.end.y)/2
            width = peak.end.x - peak.start.x
            try:
                self.peak_fits.append(
                    fit(self.x, self.y_bc, peak.start.x - width*0.1,
                        peak.end.x + width*0.1, gaussian,
                        [height, peak.apex.x, width]))
            except RuntimeError:
                print("Could not fit peak at:", peak.apex.x, file=sys.stderr)
                self.peak_fits.append(None)


    def generate_y_fit(self):
        self.y_fit = np.full_like(self.y, 0)
        for fit_params in self.peak_fits:
            if fit_params == None:
                break
            self.y_fit += np.array([gaussian(xi, *fit_params) for xi in self.x])
        return self.y_fit


    def plot_on(self, ax):
        ax.plot(self.x, self.y, label="original")
        ax.legend()


class CSIntegrator:
    "Mock-up of the ChemStation Integrator."
    sampling_ratio = 15  # number of points per peak


    def __init__(self, cdata, **parameters):
        self.dt = cdata.dt
        self.x = np.array(cdata.x)
        self.y = np.array(cdata.y)
        self.peaks = []
        self.peak_table = []

        # define the default values for the initial events
        self.events = {"slope_sensitivity": 0.1,
                       "peak_width": 3,
                       "area_reject": 0.0,
                       "area%_reject": 0.0,
                       "height_reject": 0.0}

        # overwrite these default values with the values given to the constructor
        for key, value in parameters.items():
            self.events[key] = value

        # timed events are given as tuples (event_name, parameter) if no parameter is required it is set to None
        self.timed_events = [[300, ("peak_width", 5)]]
        self.timed_events.sort(key=lambda x: x[0])


    def find_peaks(self):
        """Apply the integration algorithm.
        Returns the sampled data for inspection,
        analytical results are saved in the class instance.
        """
        # go through the data from left to right and sample the data according to peak_width
        i = 0
        step = int(self.events["peak_width"]/CSIntegrator.sampling_ratio/self.dt) or 1
        halfstep = int(step/2) or 1

        in_peak = past_apex = past_inflection = False
        start = apex = end = False
        last_slope = -1

        while i + 2*step < len(self.x):
            # check if a timed event needs to be activated
            for e in self.timed_events:
                if e[0] < self.x[int(i+halfstep)]:
                    self.events[e[1][0]] = e[1][1]

            # define the step size according to the current peak_width
            step = int(self.events["peak_width"]/CSIntegrator.sampling_ratio/self.dt) or 1
            this_y = np.mean(self.y[i:i+step])
            next_y = np.mean(self.y[i+step:i+2*step])
            slope = (next_y - this_y)/(step*self.dt)

            if in_peak:
                if past_inflection: # detect end
                    if slope > -self.events["slope_sensitivity"]:
                        end = Point(self.x[i+halfstep], self.y[i+halfstep])
                        self.peaks.append(Peak(start, apex, end))
                        in_peak = past_apex = past_inflection = False
                elif past_apex: # detect peak inflection after apex
                    if last_slope < slope:
                        past_inflection = True
                else: # detect peak apex
                    if slope <= 0:
                        apex = Point(self.x[i+halfstep], self.y[i+halfstep])
                        past_apex = True
            else: # detect peak start
                if slope >= self.events["slope_sensitivity"]:
                    start = Point(self.x[i+halfstep], self.y[i+halfstep])
                    in_peak = True

            last_slope = slope
            i += step
        self._build_peak_table()
        return self.peaks


    def _build_peak_table(self):
        cumulative_area = 0
        for peak in self.peaks:
            peak_area = integrate_peak(peak, self.x, self.y)
            cumulative_area += peak_area
            self.peak_table.append(dict(RT=peak.apex.x, Area=peak_area))
        for peak in self.peak_table:
            peak['Area%'] = 100 * peak['Area'] / cumulative_area


    def print_peak_table(self):
        print("Peak Table")
        for peak in self.peak_table:
            print("RT: {RT:6.2}\tArea: {Area:6.3}\tArea%: {Area%:6.2}".format(**peak))


    def plot_on(self, ax):
        ax.plot(self.x, self.y)
        for i, peak in enumerate(self.peaks):
            ax.plot([peak.start.x, peak.end.x], [peak.start.y, peak.end.y],
                    'k', linewidth = 0.5)
            ax.text(peak.apex.x, peak.apex.y + 0.1,
                    "{:.0f}\n{:4.2f}".format(peak.apex.x, self.peak_table[i]["Area%"]),
                    rotation = 90, horizontalalignment='center')

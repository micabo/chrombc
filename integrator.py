# -*- coding: utf-8 -*-
"""Evaluation of chromatographic data with different integrator.
All times in seconds.
"""

from cutil import ChromData, Point, Peak, gaussian, fit_peak
from collections import namedtuple
from scipy.signal import find_peaks as sfind_peaks
import numpy as np
import pdb


#-----------------------------------------------------------------------------
Event = namedtuple("Event", ["name", "time", "value"])


#-----------------------------------------------------------------------------
class Integrator:
    """The typical interface of an integrator"""
    def __init__(self):
        self.cdata = []
        self.parameters = {}
        self.bl_parameters = {}


    def __getitem__(self, index):
        return self.cdata[index]
    
    
    def __iter__(self):
        return iter(self.cdata)


    def __len__(self):
        return len(self.cdata)


    def add_chromatogramm(self, data):
        if isinstance(data, ChromData):
            self.cdata.append(data)
        elif isinstance(data, list) and len(data) > 0:
            self.add_chromatogramm(data[0])
            self.add_chromatogramm(data[1:])
        else:
            raise TypeError("Function expects argument of type ChromData "
                            "or a list of ChromData")


    def find_peaks(self, index):
        pass


    def find_baseline(self, index):
        pass
        
    
    def subtract_baseline(self, index):
        pass
    
    
    def fit_peaks(self, index):
        pass
    
    
    def build_chromatogram_from_fits(self, index):
        pass
        


#-----------------------------------------------------------------------------
class CSIntegrator(Integrator):
    "Mock-up of the ChemStation Integrator."
    sampling_ratio = 15  # number of points per peak


    def __init__(self, **kwargs):
        Integrator.__init__(self)

        # define the default values for the initial events
        default_events = {"slope_sensitivity": 1.0 * 60,
                          "peak_width": 0.05 * 60,
                          "area_reject": 0.0,
                          "area%_reject": 0.0,
                          "height_reject": 0.0}

        # overwrite these default values with the values given to the constructor
        self.parameters = {**default_events, **kwargs}

        # timed events are given as tuples (event_name, parameter) if no parameter is required it is set to None
        self.timed_events = [[300, ("peak_width", 5)]]
        self.timed_events.sort(key=lambda x: x[0])


    def find_peaks(self, index):
        """Apply the integration algorithm.
        Returns the sampled data for inspection,
        analytical results are saved in the class instance.
        """
        assert index < len(self.cdata)
        x, y = self.cdata[index]
        dt = self.cdata[index].dt

        # go through the data from left to right and sample the data according to peak_width
        i = 0
        step = int(self.parameters["peak_width"]/CSIntegrator.sampling_ratio/dt) or 1
        halfstep = int(step/2) or 1

        in_peak = past_apex = past_inflection = False
        start = apex = end = False
        last_slope = -1

        while i + 2*step < len(x):
            # check if a timed event needs to be activated
            for event in self.timed_events:
                time, event = event
                event, parameter = event
                if time < x[int(i+halfstep)]:
                    self.parameters[event] = parameter

            # define the step size according to the current peak_width
            step = int(self.parameters["peak_width"]/CSIntegrator.sampling_ratio/dt) or 1
            this_y = np.mean(y[i:i+step])
            next_y = np.mean(y[i+step:i+2*step])
            slope = (next_y - this_y)/(step*dt)

            if in_peak:
                if past_inflection: # detect end
                    if slope > -self.parameters["slope_sensitivity"]:
                        end = Point(x[i+halfstep], y[i+halfstep])
                        self.cdata[index].add_peak(Peak(start, apex, end))
                        in_peak = past_apex = past_inflection = False
                elif past_apex: # detect peak inflection after apex
                    if last_slope < slope:
                        past_inflection = True
                else: # detect peak apex
                    if slope <= 0:
                        apex = Point(x[i+halfstep], y[i+halfstep])
                        past_apex = True
            else: # detect peak start
                if slope >= self.parameters["slope_sensitivity"]:
                    start = Point(x[i+halfstep], y[i+halfstep])
                    in_peak = True

            last_slope = slope
            i += step

        self.cdata[index].build_peak_table()
        return self.cdata[index].peaks


#-----------------------------------------------------------------------------
class MBIntegrator(Integrator):

    win_width = 21 # window_width for smoothing of derivatives

    # self.point_types - a dictionary of possible point types
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


    def __init__(self, **kwargs):
        Integrator.__init__(self)

        default_values = {
            "threshold": 0.1,
            "min_width": 1,
            "min_height": 0.5
            }
        self.parameters = {**default_values, **kwargs}

        self.bl_events = []
        self.bl_events.append(Event("threshold", 0, 0.01))
        self.bl_events.append(Event("width", 0, 500))
        self.bl_events.append(Event("threshold", 1200, 1))


    def find_peaks(self, index):
        assert index < len(self.cdata)
        x, y = self.cdata[index]
        self.cdata[index].clear_peaks()
        # get the first and second derivative
        dy, ddy = self.cdata[index].get_derivatives(self.win_width)

        point_type = np.full_like(x, 0, dtype=np.int8)

        N = len(self.cdata[index])
        N_start = int((MBIntegrator.win_width + 1)/2) or 1
        N_end = N - N_start

        inflection_1_found = False
        inflection_2_found = False

        i = N_start
        while i < N_end:
            apex = start = end = False
            if dy[i]*dy[i-1] < 0 and ddy[i] < 0:
                point_type[i-1] = self.point_types["apex"]
                apex = Point(x[i-1], y[i-1])
                j = i - 2
                # go to the left and find peak start
                while j > N_start:
                    point_type[j] = self.point_types["peak-rise"]
                    if ddy[j-1]*ddy[j] < 0: # found inflection
                        point_type[j-1] = self.point_types["inflection-1"]
                        inflection_1_found = True
                    if dy[j-1] < self.parameters["threshold"] and inflection_1_found:
                        point_type[j] = self.point_types["peak-start"]
                        start = Point(x[j], y[j])
                        break
                    j -= 1
                # go to the right and find peak end
                k = i
                while k < N_end:
                    point_type[k] = self.point_types["peak-fall"]
                    if ddy[k+1]*ddy[k] < 0: # found inflection
                        point_type[k+1] = self.point_types["inflection-2"]
                        inflection_2_found = True
                    if dy[k+1] > -self.parameters["threshold"] and inflection_2_found:
                        point_type[k] = self.point_types["peak-end"]
                        end = Point(x[k], y[k])
                        break
                    k += 1
                # TODO: now go and find the next peaks start, if close enough, merge end and start and use the next peak end for height determination

                # reject peaks which are too small
                if apex and start and end:
                    width = end.x - start.x
                    height = apex.y - (start.y + end.y)/2
                    if width > self.parameters["min_width"] and height > self.parameters["min_height"]:
                        self.cdata[index].add_peak(Peak(start, apex, end))

                # when start and end have been found:
                # continue searching for peaks after the end of the current
                # peak (k + 1)
                inflection_1_found = False
                inflection_2_found = False
                i = k + 1

            i += 1

        # TODO: last order of business - merge close lying end/start points of adjacent peaks

        self.cdata[index].build_peak_table()
        return self.cdata[index].peaks


    def find_baseline(self, index):
        # should implement timed events
        width = 500

        x, y = self.cdata[index]
        x_min = [x[0]]
        y_min = [y[0]]
        i_min = 0
        i = 1
        while i < len(y):
            if y[i] < y[i_min]:
                i_min = i
            if i % width == 0:
                x_min.append(x[i_min])
                y_min.append(y[i_min])
                i_min = i
            i += 1

        # clean baseline -> i.e. search for "peaks"
        i = 0
        #pdb.set_trace()

        while i < len(y_min) - 1:
            j = i + 1

            for event in self.bl_events:
                if event.time <= x_min[i]:
                    # assume strictly ordered events
                    self.bl_parameters[event.name] = event.value
                else:
                    break

            while j < len(y_min):
                slope = (y_min[j] - y_min[i]) / (x_min[j] - x_min[i])
                if abs(slope) > self.bl_parameters["threshold"]:
                    j += 1
                    continue

                # slope is less than threshold
                if j - i > 1:
                    k = i + 1
                    while k < j:  # fill in interpolated values
                        y_min[k] = y_min[i] + slope * (x_min[k] - x_min[i])
                        k += 1
                i = j
                break
            else:
                # end of j loop is reached without fulfilling the above criteria
                break

        return x_min, y_min
        
        
    def subtract_baseline(self, index):
        x, y = self.cdata[index]
        xb, yb = self.find_baseline(index)
        ybc = np.full_like(y, 0)
        
        i = 0
        while i < len(ybc):
            # TODO: implement baseline correction
            pass
        
        return x, ybc
        
        
    def fit_peaks(self, index):
         "fit peaks with gaussian"
         xbc, ybc = subtract_baseline(index)
         for peak in self.cdata[index].peaks:
             height = peak.apex.y - (peak.start.y + peak.end.y)/2
             width = peak.end.x - peak.start.x
             try:
                # TODO: adapt to new situation -> where are the fits stored?
                 self.peak_fits.append(
                     fit_peak(xbc, ybc, peak.start.x - width*0.1,
                         peak.end.x + width*0.1, gaussian,
                         [height, peak.apex.x, width]))
             except RuntimeError:
                 print("Could not fit peak at:", peak.apex.x, file=sys.stderr)
                 self.peak_fits.append(None)


    def build_chromatogram_from_fits(self, index):
        x, y = self.cdata[index]
        y_fit = np.full_like(y, 0)
        for fit_params in self.peak_fits:
            if fit_params == None:
                break
            y_fit += np.array([gaussian(xi, *fit_params) for xi in x])
        return y_fit
    


class ScIntegrator(Integrator):
    
    def __init__(self):
        Integrator.__init__(self)
        
    
    def find_peaks(self, index):
        assert index < len(self.cdata)
        x, y = self.cdata[index]
        self.cdata[index].clear_peaks()
        
        # does not do what I want yet...
        peak_indices, details = sfind_peaks(y, prominence=0.1)
        p = details["prominences"]
        l = details["left_bases"]
        r = details["right_bases"]
        
        for start_i, apex_i, end_i in zip(l, peak_indices, r):
            start = Point(x[start_i], y[start_i])
            apex = Point(x[apex_i], y[apex_i])
            end = Point(x[end_i], y[end_i])
            self.cdata[index].add_peak(Peak(start, apex, end))
        
        self.cdata[index].build_peak_table()
        return self.cdata[index].peaks


#-----------------------------------------------------------------------------
if __name__ == "__main__":
    m = ScIntegrator()
    m.add_chromatogramm(ChromData("./data/SST.txt"))
    x, y, peak_indices, details = m.find_peaks(0)
    
    p = details["prominences"]
    l = details["left_bases"]
    r = details["right_bases"]
    
    import matplotlib.pyplot as plt
    plt.plot(x, y)
    for i, index in enumerate(peak_indices):
        x_p = x[index]
        y_p = y[index]
        plt.annotate(i,
                     xy=(x_p, y_p),
                     xytext=(x_p, y_p * 1.1),
                     fontsize='small',
                     horizontalalignment='center',
                     rotation=90)
    for index in l:
        plt.plot(x[index], y[index], 'rx')
    
    for index in r:
        plt.plot(x[index], y[index], 'go')
        
    plt.show()
    
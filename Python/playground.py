# LIBRARIES IMPORT

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

from scipy import interpolate

# LIBRARIES IMPORT

import os
import pandas as pd
import numpy as np
import scipy as sp
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation

from scipy import signal                    # for filtering
from scipy import interpolate               # for interpolation
from mpl_toolkits.mplot3d import Axes3D     # for 3D visualization

# SETTINGS

#   Main Directory and Data type
DATADIR = "C:/Users/bdour/Documents/Data"
data_type =  'Vicon'

#   Participant number
participant = "01"

#   Trial
trial = "RDist_0"

#   Extensions
ext_data = '.csv'       # extension for data files

#   Paths

#       Import
# path for participant folder
path = os.path.join(DATADIR, data_type + '\\' + participant)
# path for dynamic file
dynamic = os.path.join(path, 'Nexus/' + data_type + '- ' + participant + '_' + trial + ext_data)

#       Dynamic
d_data = pd.read_csv(dynamic, skiprows=5)
d_data = np.array(d_data)
d_data = d_data[1:,2:]/1000
idx = int(np.shape(d_data)[1]/3)

pos_d = d_data[:,0:idx]

for i in range(1,int(idx/3)+1):
    exec ('pd_m' + str(i) + '=' + str('np.array(pos_d[:,') + str(3*(i-1)) + ':' + str(3*i) + '])')

# FUNCTION

def nan_find(y):
    '''
    Generates a NaNs logical array where the indices of each NaN observation is True.
    Generates a local function that can extract the indices of each NaN observation as a list.
    Input:
        - y: nx1 array that contains NaNs
    Output:
        - nan_logic: logical array where the indices of each NaN observation is True
        - find_true: function that returns the indices of all True observations in an array.
    Example:
        nan_logic, find_true = nan_find(y)
        find_true(nan_logic) -> returns array with indices of all NaN in y
        find_true(~nan_logic) -> returns array with indices of all non-NaN in y
    '''
    nan_logic = np.isnan(y)
    # lambda k: k.nonzero()[0] defines the function find_true, and k represents the corresponding input
    find_true = lambda k: k.nonzero()[0]

    return nan_logic, find_true

def cubic_spline_fill(x, y):
    '''
    Assesses if time series y has missing observations (i.e. NaN).
    If yes, interpolates missing observations using a cubic spline fitted to the data.
    Input:
        x: nx1 array corresponding to the frames indices
        y: nx1 array: tested time series (e.g. x coordinates of a marker)
    Output:
        y_interp: new nx1 array with interpolated values replacing NaNs (only returned if data has NaNs)
        Notes:
            - Only interpolates time series with NaN (otherwise return original time series)
            - Does not interpolate empty time series (returns same empty time series)
            - Does not interpolate the beginning and/or end of the time series if it has missing observations (i.e. only interpolates between edges)
    Dependencies:
        nan_find
    '''
    # if no missing observation, return original signal
    if np.shape(np.where(np.isnan(y) == True))[1] == 0:

        return y

    # if empty observation, return original signal
    elif np.shape(np.where(np.isnan(y) == True))[1] == np.shape(y)[0]:

        return y

    else:
        # generate a NaNs logical array where the indices of each NaN observation is True
        nan_logic, find_true = nan_find(y)
        # find corresponding indices of missing observations
        obs = find_true(~nan_logic)
        # isolate non-missing portion of the signal
        a = x[obs]
        b = y[obs]
        # find the equation of the cubic spline that best fits the corresponding signal
        cs = interpolate.CubicSpline(a, b)
        # initiate y_interp as an empty array
        y_interp = np.array(np.empty(np.shape(y)[0]))
        # fill y_interp with NaNs
        y_interp[:] = np.nan
        # apply cubic spline equation to interpolate the whole signal (- missing observations at edges)
        y_interp[obs[0]:obs[-1]] = cs(x[obs[0]:obs[-1]])
        # to avoid unstable edges, NaNs are applied to the two data points at the two edges of the signal
        y_interp[obs[0]] = np.nan
        y_interp[obs[-1]] = np.nan

        return y_interp

x = np.array(range(0,int(len(pos_d))))
y = pd_m24[:,0]
y[43:76] = np.nan
y[154:192] = np.nan
y[245:173] = np.nan
y[347:387] = np.nan
y[423:455] = np.nan
y_interp = cubic_spline_fill(x, y)

# Create a figure canvas
fig, ax = plt.subplots()
# Plot the interpolated signal
ax.plot(x, y_interp, label='Interpolated', color='orange')
# Plot the original data with missing observations
ax.plot(x, y, label='Original')
# Add legend to plot
ax.legend()
plt.show()

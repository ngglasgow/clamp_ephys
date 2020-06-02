import os
import numpy as np
import pandas as pd
from neo.io import IgorIO


def igor_to_pandas(path_to_file):
    '''This function opens an igor binary file (.ibw), extracts the time
    series data, and returns a pandas DataFrame'''

    data_raw = IgorIO(filename=path_to_file)
    data_neo = data_raw.read_block()
    data_neo_array = data_neo.segments[0].analogsignals[0]
    data_df = pd.DataFrame(data_neo_array.as_array().squeeze())

    return data_df


def new_mean_baseline(data, fs, baseline_start, baseline_end=6000):
    '''
    Find the mean baseline in a given time series, defined as the last 3s of the sweep
    Parameters
    ----------
    data: pandas.Series or pandas.DataFrame
        The time series data for which you want a baseline.
    fs: int or float
        The sampling frequency in kHz.
    baseline_start: int or float
        The time in ms when baseline starts.
    baseline_end: int or float
        The length of the sweep and the time when baseline ends.

    Returns
    -------
    baseline: float or pandas.Series
        The mean baseline over the defined window
    '''
    start = baseline_start * fs
    stop = baseline_end * fs
    window = data.iloc[start:stop]
    baseline = window.mean()

    return baseline


def new_std_baseline(data, fs, baseline_start, baseline_end=6000):
    '''
    Find the mean baseline in a given time series
    Parameters
    ----------
    data: pandas.Series or pandas.DataFrame
        The time series data for which you want a baseline.
    fs: int or float
        The sampling frequency in kHz.
    baseline_start: int or float
        The time in ms when baseline starts.
    baseline_end: int or float
        The length of the sweep and the time when baseline ends.

    Returns
    -------
    baseline: float or pandas.Series
        The mean baseline over the defined window
    '''
    start = baseline_start * fs
    stop = baseline_end * fs
    window = data.iloc[start:stop]
    std = window.std()

    return std


def mean_baseline(data, fs, stim_time, pre_stim=100):
    '''
    Find the mean baseline in a given time series
    Parameters
    ----------
    data: pandas.Series or pandas.DataFrame
        The time series data for which you want a baseline.
    fs: int or float
        The sampling frequency in kHz.
    stim_time: int or float
        The time in ms when stimulus is triggered.
    pre_stim: int or float
        Time in ms before the stimulus trigger over which baseline is measured.

    Returns
    -------
    baseline: float or pandas.Series
        The mean baseline over the defined window
    '''
    start = (stim_time - pre_stim) * fs
    stop = (stim_time - 1) * fs
    window = data.iloc[start:stop]
    baseline = window.mean()

    return baseline


def std_baseline(data, fs, stim_time):
    '''
    Find the mean baseline in a given time series
    Parameters
    ----------
    data: pandas.Series or pandas.DataFrame
        The time series data for which you want a baseline.
    fs: int or float
        The sampling frequency in kHz.
    stim_time: int or float
        The time in ms when stimulus is triggered.
    pre_stim: int or float
        Time in ms before the stimulus trigger over which baseline is measured.

    Returns
    -------
    baseline: float or pandas.Series
        The mean baseline over the defined window
    '''
    start = 100 * fs
    stop = (stim_time - 1) * fs
    window = data.iloc[start:stop]
    std = window.std()

    return std


def epsc_peak(data, baseline, fs, stim_time, post_stim=100, polarity='-', index=False):
    '''
    Find the peak EPSC value for a pandas.Series or for each sweep (column) of
    a pandas.DataFrame. This finds the absolute peak value of mean baseline
    subtracted data.

    Parameters:
    -----------
    
    data: pandas.Series or pandas.DataFrame
        Time series data with stimulated synaptic response triggered at the
        same time for each sweep.
    baseline: scalar or pandas.Series
        Mean baseline values used to subtract for each sweep.
    fs: int or float
        The sampling frequency in kHz. Default is 10 kHz.
    stim_time: int or float
        Time in ms at which stimulus is triggered each sweep.
    polarity: str
        The expected polarity of the EPSC; negative: '-'; postitive: '+'.
        Default is '-'.
    post_stim: int or float
        Time in ms that marks the end of the sampling window post stimulus.
        Default is 100 ms.
    index: bool
        Determines whether or not to return the peak index in addition to the peak. 


    Returns
    -------
    epsc_peaks: pandas.Series
        The absolute peak of mean baseline subtracted time series data.
    epsc_peak_index: int
        The time at which the peak occurs
    '''

    subtracted_data = data - baseline
    start = stim_time * fs
    end = (stim_time + post_stim) * fs
    peak_window = subtracted_data.iloc[start:end]

    if index is True:
        if polarity == '-':
            epsc_peaks = peak_window.min()
            epsc_peaks_index = peak_window.idxmin()
        elif polarity == '+':
            epsc_peaks = peak_window.max()
            epsc_peaks_index = peak_window.idxmax()
        else:
            raise ValueError(
                "polarity must either be + or -"
            )    
        return epsc_peaks, epsc_peaks_index

    elif index is False:
        if polarity == '-':
            epsc_peaks = peak_window.min()
        elif polarity == '+':
            epsc_peaks = peak_window.max()
        else:
            raise ValueError(
                "polarity must either be + or -"
            )    
        return epsc_peaks
        

def series_resistance(data, fs, tp_start=5, vm_jump=10, pre_tp=3, unit_scaler=-12):
    '''
    Calculate the approximate series resistance (Rs) from a test pulse (tp).

    Parameters
    ----------
    data: pandas.Series or pandas.DataFrame
        Raw time series daata of the v-clamp recording in nA.
    fs: int or float
        Sampling frequency in kHz.
    tp_start: int or float
        Time in ms when test pulse begins. Default is 5.
    vm_jump: int or float
        Amplitude of whatever windows needs here the test pulse voltage command in mV.
        This is 10 mV in MIES and -5 in Nathan's Igor program.
    pre_tp: int or float
        Time in ms before start of test pulse by which to measure the baseline.
    unit_scaler: int or float
        Scaling factor to convert reported current as Amp; e.g. if units are in pA value is -12.

    Returns:
    rs: pandas.Series of float
        The series resistance for each sweep in MOhms.
    '''
    # find the baseline 5 ms pre test pulse and subtract from raw data
    rs_baseline = mean_baseline(fs=fs, data=data, stim_time=tp_start, pre_stim=pre_tp)
    rs_subtracted = data - rs_baseline

    # set up indices for starting and ending peak window
    start = tp_start * fs
    end = (tp_start + 2) * fs
    rs_window = rs_subtracted.iloc[start:end]

    if vm_jump > 0:
        rs_peak = rs_window.max()
    else:
        rs_peak = rs_window.min()

    # calculate Rs via V=IR -> Rs = V/I
    rs = ((vm_jump * 10**-3) / (rs_peak * 10**unit_scaler)) * 10**-6

    return rs


def decay_func(time, current_peak, tau, offset):
    '''
    Exponential decay function for calculating tau
    Parameters
    ----------
    time: array
        x array of time
    current_peak: scalar/float
        value of the starting peak for the exponential decay
    tau: scalar/float
        the time constant of decay of the exponential
    offset: scalar/float
        the asymptote that the exponential will decay to

    Returns
    -------
    y: array
        the y values for the exponential decay
    '''
    y = current_peak * np.exp(-time/tau) + offset
    return y
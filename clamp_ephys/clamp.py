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


def mean_baseline(sf, data, stim_time, pre_stim=100):
    '''
    Find the mean baseline in a given time series
    Parameters
    ----------
    sf: int or float
        The sampling frequency in kHz.
    data: pandas.Series or pandas.DataFrame
        The time series data for which you want a baseline.
    stim_time: int or float
        The time in ms when stimulus is triggered.
    pre_stim: int or float
        Time in ms before the stimulus trigger over which baseline is measured.

    Returns
    -------
    baseline: float or pandas.Series
        The mean baseline over the defined window
    '''
    start = (stim_time - pre_stim) * sf
    stop = (stim_time - 1) * sf
    window = data.iloc[start:stop]
    baseline = window.mean()

    return baseline


def std_baseline(sf, data, stim_time, pre_stim=100):
    '''
    Find the baseline st. dev. in a given time series
    Parameters
    ----------
    data: pandas.Series or pandas.DataFrame
        The time series data for which you want a baseline.
    stim_time: int or float
        The time in ms when stimulus is triggered.
    pre_stim: int or float
        Time in ms before the stimulus trigger over which baseline is measured.
    sf: int or float
        The sampling frequency in kHz.

    Returns
    -------
    std: float or pandas.Series
        The baseline st. dev. over the defined window
    '''
    start = (stim_time - pre_stim) * sf
    stop = (stim_time - 1) * sf
    window = data.iloc[start:stop]
    std = window.std()

    return std


def epsc_peak(sf, data, baseline, stim_time, polarity='-', post_stim=100):
    '''
    Find the peak EPSC value for a pandas.Series or for each sweep (column) of
    a pandas.DataFrame. This finds the absolute peak value of mean baseline
    subtracted data.

    Parameters:
    -----------
    sf: int or float
        The sampling frequency in kHz. Default is 10 kHz.
    data: pandas.Series or pandas.DataFrame
        Time series data with stimulated synaptic response triggered at the
        same time for each sweep.
    baseline: scalar or pandas.Series
        Mean baseline values used to subtract for each sweep.
    stim_time: int or float
        Time in ms at which stimulus is triggered each sweep.
    polarity: str
        The expected polarity of the EPSC; negative: '-'; postitive: '+'.
        Default is '-'.
    post_stim: int or float
        Time in ms that marks the end of the sampling window post stimulus.
        Default is 100 ms.


    Returns
    -------
    epsc_peaks: pandas.Series
        The absolute peak of mean baseline subtracted time series data.
    '''

    subtracted_data = data - baseline
    start = stim_time * sf
    end = (stim_time + post_stim) * sf
    peak_window = subtracted_data.iloc[start:end]

    if polarity == '-':
        epsc_peaks = peak_window.min()
    elif polarity == '+':
        epsc_peaks =peak_window.max()
    else:
        raise ValueError(
            "polarity must either be + or -"
        )    

    return epsc_peaks


def series_resistance(sf, data, tp_start=5, vm_jump=10, pre_tp=3, peak_factor=-12):
    '''
    Calculate the approximate series resistance (Rs) from a test pulse (tp).

    Parameters
    ----------
    sf: int or float
        Sampling frequency in kHz. Default is 10 kHz.
    data: pandas.Series or pandas.DataFrame
        Raw time series daata of the v-clamp recording in nA.
    tp_start: int or float
        Time in ms when test pulse begins. Default is 5.
    vm_jump: int or float
        Amplitude of whatever windows needs here the test pulse voltage command in mV.
        This is 10 mV in MIES and -5 in Nathan's Igor program.
    pre_tp: int or float
        Time in ms before start of test pulse by which to measure the baseline.
    peak_factor: int or float
        Factor to multiply rs_peak by to get pA. Default is -12.

    Returns:
    rs: pandas.Series of float
        The series resistance for each sweep in MOhms.
    '''
    # find the baseline 5 ms pre test pulse and subtract from raw data
    rs_baseline = mean_baseline(sf=sf, data=data, stim_time=tp_start, pre_stim=pre_tp)
    rs_subtracted = data - rs_baseline

    # set up indices for starting and ending peak window
    start = tp_start * sf
    end = (tp_start + 2) * sf
    rs_window = rs_subtracted.iloc[start:end]

    if vm_jump > 0:
        rs_peak = rs_window.max()
    else:
        rs_peak = rs_window.min()

    # calculate Rs via V=IR -> Rs = V/I
    rs = ((vm_jump * 10**-3) / (rs_peak * 10**peak_factor)) * 10**-6

    return rs
# -*- coding: utf-8 -*-
"""
Created 25 March 2019
Nathan Glasgow

This module defines simple voltage clamp analysis functions for use with Igor
binary files.

"""
# from __main__ import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import elephant
from neo.io import IgorIO
import os
from collections import OrderedDict
import math


def get_metadata(file, data_notes):
    '''Takes a filename and parses it for metadata, and returns metadata in an
    orderedDict as a pandas DataFrame for saving later
    Also takes information from the cell spreadsheet in data_notes'''

    # pull out cell id, cell number, date and condition
    file_split = file.split('_')
    cell_id = file_split[0]+'_'+file_split[1]
    cell_num = cell_id[-1:]
    date = '20'+cell_id[2:4]+'-'+cell_id[4:6]+'-'+cell_id[6:8]

    if 'drug' in file:
        condition = 'TTX+4-AP'
    else:
        condition = 'control'

    # grab metadata from data notes spreadsheet
    file_data = data_notes[data_notes['Cell name'] == file]
    cell_path = file_data['File Path'].tolist()[0]
    genotype = file_data['Genotype'].tolist()[0]
    cell_type = file_data['Cell type'].tolist()[0]
    depol_sweep_start = file_data['Depol sweeps start'].tolist()[0]
    depol_sweep_stop = file_data['Depol sweeps stop'].tolist()[0]

    # save metadate into orderedDict pandas DataFrame
    dict = OrderedDict()
    dict['Date'] = date
    dict['Cell ID'] = cell_id
    dict['Cell Number'] = cell_num
    dict['Cell Path'] = cell_path
    # dict['Condition'] = condition
    dict['Genotype'] = genotype
    dict['Cell Type'] = cell_type
    dict['Exclude Sweep Start'] = depol_sweep_start
    dict['Exclude Sweep Stop'] = depol_sweep_stop
    metadata = pd.DataFrame(dict, index=range(1))

    return metadata


def igor_to_pandas(file, data_dir):
    '''This function opens an igor binary file (.ibw), extracts the time
    series data, and returns a pandas DataFrame'''

    file_path = os.path.join(data_dir, file)
    data_raw = IgorIO(filename=file_path)
    data_neo = data_raw.read_block()
    data_neo_array = data_neo.segments[0].analogsignals[0]
    data_df = pd.DataFrame(data_neo_array.as_array())

    return data_df


def mean_baseline(data, stim_time, pre_stim=100, sf=10):
    '''
    Find the mean baseline in a given time series
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
    baseline: float or pandas.Series
        The mean baseline over the defined window
    '''
    start = (stim_time - pre_stim) * sf
    stop = (stim_time - 1) * sf
    window = data.iloc[start:stop]
    baseline = window.mean()

    return baseline


def epsc_peak(data, baseline, stim_time, polarity='-', post_stim=100, sf=10):
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
    stim_time: int or float
        Time in ms at which stimulus is triggered each sweep.
    polarity: str
        The expected polarity of the EPSC; negative: '-'; postitive: '+'.
        Default is '-'.
    post_stim: int or float
        Time in ms that marks the end of the sampling window post stimulus.
        Default is 100 ms.
    sf: int or float
        The sampling frequency in kHz. Default is 10 kHz.

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


def series_resistance(data, tp_start, vm_jump, sf=10):
    '''
    Calculate the approximate series resistance (Rs) from a test pulse (tp).

    Parameters
    ----------
    data: pandas.Series or pandas.DataFrame
        Raw time series daata of the v-clamp recording in nA.
    tp_start: int or float
        Time in ms when test pulse begins.
    vm_jump: int or float
        Amplitude ofwhatever windows needs here the test pulse voltage command in mV..
    sf: int or float
        Sampling frequency in kHz. Default is 10 kHz.

    Returns:
    rs: pandas.Series of float
        The series resistance for each sweep in MOhms.
    '''
    # find the baseline 10 ms pre test pulse and subtract from raw data
    rs_baseline = mean_baseline(data, stim_time=tp_start, pre_stim=11)
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
    rs = ((vm_jump * 10**-3) / (rs_peak * 10**-9)) * 10**-6

    return rs

# -*- coding: utf-8 -*-
"""
Created 5 March 2019
epsc_peak_x.y.z.py

"""
# from __main__ import *

import pandas as pd
import numpy as np
import elephant
from neo.io import IgorIO
import os
import platform


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


def std_baseline(data, stim_time, pre_stim=100, sf=10):
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


def epsc_peak(data, baseline, stim_time, polarity='-', post_stim=100, sf=10):
    '''
    Find the peak EPSC value for a pandas.Series or for each sweep (column) of
    a pandas.DataFrame. This finds the absolute peak value of mean baseline
    subtracted data.

    Creates a 10 ms window from response start time and asks whether the polarity
    of the response changes during that time window. If so, the response is not
    sustained, and sustained returns "False".

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

    sustained: str
        Whether or not the response is sustained over 10 ms, determined by 
        whether a negative value exists if the polarity of the amplitude
        is positive, or a positive value exists if the polarity is negative.
    '''

    subtracted_data = data - baseline
    start = stim_time * sf
    end = (stim_time + post_stim) * sf
    peak_window = subtracted_data.iloc[start:end]

    if polarity == '-':
        epsc_peaks = peak_window.min()
        response_start = peak_window.idxmin()
    elif polarity == '+':
        epsc_peaks = peak_window.max()
        response_start = peak_window.idxmax()
    else:
        raise ValueError(
            "polarity must either be + or -"
        )

    response_end = response_start + (10*sf)
    response_window = subtracted_data.iloc[response_start:response_end]
    response_window.columns = ['Amplitude']

    isthere_negative = response_window < 0
    isthere_positive = response_window > 0

    if polarity == '-':
       boolean_counts = isthere_positive['Amplitude'].value_counts()
       
    elif polarity == '+':
       boolean_counts = isthere_negative['Amplitude'].value_counts()
    else:
        raise ValueError(
            "polarity must either be + or -"
        )
    
    cross_0 = boolean_counts[True]

    if cross_0 >= 1:
        sustained = 'False'
    else:
        sustained = 'True'

    return epsc_peaks, sustained


''' *********************************************************************** '''

''' ################## Define file structure on server #################### '''
# home_dir will depend on the OS, but the rest will not
# query machine identity and set home_dir from there
machine = platform.uname()[0]

if machine == 'Darwin':
    home_dir = '/Volumes/Urban'

elif machine == 'Linux':
    home_dir = '/run/user/1000/gvfs/smb-share:server=130.49.237.41,share=urban'

elif machine == 'Windows':
    home_dir = r"N:\urban\Huang"

else:
    print("OS not recognized. \nPlease see Nate for correction.")

project_dir = os.path.join(home_dir, 'Injected_GC_data', 'New_VC_pairs')
figure_dir = os.path.join(project_dir, 'figures')
table_dir = os.path.join(project_dir, 'tables')
data_dir = os.path.join(project_dir, 'data')

''' ## Open the notes spreadsheet and parse for what we want to analyze ## '''
# open metadata file
allcells_data_notes = pd.read_csv(os.path.join(table_dir, 'MC_TC_summary.csv'))

# pull out cell_id for directory, file name, and make the full path
file_name_list = allcells_data_notes['Cell name'].tolist()
cell_id_list = []

for file in file_name_list:
    file_split = file.split('_')
    cell_id = file_split[0]+'_'+file_split[1]
    cell_id_list.append(cell_id)

file_path_list = []

for cell, file in zip(cell_id_list, file_name_list):
    file_path = os.path.join(file + '.ibw')
    file_path_list.append(file_path)

allcells_data_notes = pd.concat([pd.DataFrame({'File Path': file_path_list}), 
    allcells_data_notes], axis=1)
light_data_notes = allcells_data_notes[allcells_data_notes['Cell name'].str.contains("light")]
spontaneous_data_notes = allcells_data_notes[allcells_data_notes['Cell name'].str.contains("spontaneous")]

allcells_data_notes.to_csv(os.path.join(table_dir, 'allcells_data_notes.csv'))
light_data_notes.to_csv(os.path.join(table_dir, 'light_data_notes.csv'))
spontaneous_data_notes.to_csv(os.path.join(table_dir, 'spontaneous_data_notes.csv'))

''' ##########################################################################
This is all the analysis, figures, saving
Read in file metadata, open file from igor, convert to pandas
##############################################################################
'''
# make empty lists for response answers for each file
responses_df = pd.DataFrame()

# loop through all the files in file_name_list for plots and saving
for file_name, file_path in zip(file_name_list, file_path_list):
    # set file name from list
    file = file_name

    # open igor file and convert to pandas
    data = igor_to_pandas(file_path, data_dir)

    '''
        Pull out EPSC peak from unfiltered signals
        Baseline 100 ms preceding blue light
        Peak within 150 ms of blue light
    '''

    baseline = mean_baseline(data, 500)
    #peaks = epsc_peak(data, baseline, stim_time=500, post_stim=150)
    peaks, sustained_unfilt = epsc_peak(data, baseline, stim_time=500, post_stim=150)

    '''
        Pull out EPSC peaks from filtered signals
        Baseline 100 ms preceding blue light
        Peak within 150 ms of blue light
    '''

    # filter signal with butterworth filter at 500 Hz for data
    filt_data = elephant.signal_processing.butter(data.T,
                                                  lowpass_freq=500.0,
                                                  fs=10000.0)
    filt_data = pd.DataFrame(filt_data).T

    filt_baseline = mean_baseline(filt_data, 500)
    filt_baseline_std = std_baseline(filt_data, 500)
    #filt_peaks = epsc_peak(filt_data, filt_baseline, stim_time=500, post_stim=150)
    filtpeaks, sustained_filt = epsc_peak(filt_data, filt_baseline, stim_time=500, 
                                                post_stim=150)

    # make binary choice of whether the response is more than 3x std
    mean_std = filt_baseline_std.mean()
    mean_peaks = filt_peaks.mean()

    response_2x = abs(mean_peaks) > mean_std * 2
    response_3x = abs(mean_peaks) > mean_std * 3
    response_4x = abs(mean_peaks) > mean_std * 4

    if response_3x == 'True' and sustained_filt == 'True':
        true_EPSC = 'True'
    else:
        true_EPSC = 'False'

    responses = pd.DataFrame({'Cell name': file_name,
                              'Mean Peaks (nA)': mean_peaks,
                              'Mean STD (nA)': mean_std,
                              'Response 2x STD': response_2x,
                              'Response 3x STD': response_3x,
                              'Response 4x STD': response_4x,
                              'Sustained?': sustained_filt,
                              'True EPSC?': true_EPSC
                              },
                              index=range(1))

    responses_df = pd.concat([responses_df, responses], ignore_index=True)

# save responses as csv
responses_df.to_csv(os.path.join(table_dir, 'responses.csv'))

# merge with data_notes and save
responses_data_notes = pd.merge(allcells_data_notes, responses_df)
responses_data_notes.to_csv(os.path.join(table_dir, 'responses_data_notes.csv'),
                            float_format='%8.4f', index=False)

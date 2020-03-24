# -*- coding: utf-8 -*-
"""
Created 5 March 2019
epsc_peak_x.y.z.py

"""
# from __main__ import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style(rc = {'figure.facecolor':'white'})
import elephant
from neo.io import IgorIO
import os
import platform
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

    if 'light' in file:
        condition = 'light'
    else:
        condition = 'spontaneous'

    # grab metadata from data notes spreadsheet
    file_data = data_notes[data_notes['Cell name'] == file]
    cell_path = file_data['File Path'].tolist()[0]
    cell_type = file_data['Cell Type'].tolist()[0]

    # save metadate into orderedDict pandas DataFrame
    dict = OrderedDict()
    dict['Date'] = date
    dict['Cell ID'] = cell_id
    dict['Cell Number'] = cell_num
    dict['Cell Path'] = cell_path
    dict['Condition'] = condition
    dict['Cell Type'] = cell_type
    metadata = pd.DataFrame(dict, index=range(1))

    return metadata

def igor_to_pandas(file, data_dir):
    '''This function opens an igor binary file (.ibw), extracts the time
    series data, and returns a pandas DataFrame'''

    file_path = os.path.join(data_dir, file)
    data_raw = IgorIO(filename=file_path)
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

def create_data_notes(timepoint, summary_file, ibw_file_paths):
    '''
    Create data_notes summary spreadsheets
    
    Parameters
    ----------
    timepoint: str
        Name of the injection timepoint used in the analysis
    summary_file: .csv file
        Manually-generated summary file of all cells in the dataset of a timepoint
    ibw_file_paths:  list
        List of ibw files found for a timepoint

    Returns:
    file_name_list: list
        List of all the file names in the timepoint data set
    data_notes: pandas.DataFrame
        DataFrame of parsed notes for each cell from manually-inputted summary_file
    '''
    # Open the notes spreadsheet and parse for what we want to analyze ## '''
    # open metadata file
    data_notes = pd.read_csv(os.path.join(paths.tables, summary_file))

    # pull out cell_id for directory, file name, and make the full path
    file_name_list = data_notes['Cell name'].tolist()

    data_notes = pd.concat([pd.DataFrame({'File Path': ibw_file_paths}), 
        data_notes], axis=1)

    light_data_notes = data_notes[data_notes['Cell name'].str.contains("light")]
    spontaneous_data_notes = data_notes[data_notes['Cell name'].str.contains("spontaneous")]

    # data_notes.to_csv(os.path.join(paths.tables, '{}_data_notes.csv'.format(timepoint)))
    # light_data_notes.to_csv(os.path.join(paths.tables, '{}_light_data_notes.csv'.format(timepoint)))
    # spontaneous_data_notes.to_csv(os.path.join(paths.tables, '{}_spontaneous_data_notes.csv'.format(timepoint)))

    return file_name_list, data_notes

def cell_response_checker(timepoint, file_name_list, data_dir, sf=25):
    responses_df = pd.DataFrame()

    for file in file_name_list:
        metadata = get_metadata(file, data_notes)
        file_path = metadata['Cell Path'][0]
        cell_id = metadata['Cell ID'][0]
        cell_type = metadata['Cell Type'][0]
        condition = metadata['Condition'][0]

        # open igor file and convert to pandas
        data = igor_to_pandas(file_path, data_dir)

        '''
            Pull out EPSC peak from unfiltered signals
            Baseline 100 ms preceding blue light
            Peak within 250 ms of blue light
        '''

        baseline = mean_baseline(sf=25, data=data, stim_time=500)
        peaks = epsc_peak(sf=25, data=data, baseline=baseline, stim_time=500, post_stim=150)

        '''
            Pull out EPSC peaks from filtered signals
            Baseline 100 ms preceding blue light
            Peak within 250 ms of blue light
        '''

        # filter signal with butterworth filter at 500 Hz for data
        filt_data = elephant.signal_processing.butter(data.T,
                                                    lowpass_freq=500.0,
                                                    fs=10000.0)
        filt_data = pd.DataFrame(filt_data).T

        filt_baseline = mean_baseline(sf=sf, data=filt_data, stim_time=500)
        filt_baseline_std = std_baseline(sf=sf, data=filt_data, stim_time=500)
        filt_peaks = epsc_peak(sf=sf, data=filt_data, baseline=filt_baseline, stim_time=500, post_stim=150)


        # make binary choice of whether the response is more than 3x std
        mean_std = filt_baseline_std.mean()
        mean_peaks = filt_peaks.mean()

        response_2x = abs(mean_peaks) > mean_std * 2
        response_3x = abs(mean_peaks) > mean_std * 3

        responses = pd.DataFrame({'Cell name': file,
                                'Mean Peaks (pA)': mean_peaks,
                                'Mean STD (pA)': mean_std,
                                'Response 2x STD': response_2x,
                                'Response 3x STD': response_3x,
                                },
                                index=range(1))

        responses_df = pd.concat([responses_df, responses], ignore_index=True)
    return responses_df


''' *********************************************************************** '''

''' ################## Define file structure on server #################### '''
class file_structure:
    def __init__(self, location, project_path):
        '''
        Creates an object with paths as attributes:
        location:   str value 'local' or 'server' only, refers to where you are
                    doing the actual work, 'local' by default.
        project_path:   str of the root project path from wherever your home dir is
        '''
        machine = platform.uname()[0]

        if location == 'local':
            if machine == 'Darwin':
                home_dir = '/Volumes/Urban'

            elif machine == 'Linux':
                home_dir = os.path.join(os.path.expanduser('~'), 'urban/neurobio/Huang')

            elif machine == 'Windows':
                home_dir = r"C:\Users\jhuang\Documents\phd_projects"

            else:
                print("OS not recognized. \nPlease see Nate for correction.")

        elif location == 'server':
            if machine == 'Darwin':
                home_dir = '/Volumes/Urban'

            elif machine == 'Linux':
                home_dir = os.path.join(os.path.expanduser('~'), 'urban/neurobio/Huang')

            elif machine == 'Windows':
                home_dir = r"N:\Huang"

            else:
                print("OS not recognized. \nPlease see Nate for correction.")

        self.project = os.path.join(home_dir, project_path)
        self.figures = os.path.join(self.project, 'figures')
        self.tables = os.path.join(self.project, 'tables')
        self.p2 = os.path.join(self.project, 'data/p2')
        self.p2_files = os.listdir(self.p2)
        self.p14 = os.path.join(self.project, 'data/p14')
        self.p14_files = os.listdir(self.p14)

    def __repr__(self):
        return 'Project file structure and file lists for {}'.format(self.project)

paths = file_structure('local', 'Injected_GC_data/VC_pairs')

# analysis for p2 data set
# Parameters (sampling frequency = 25 kHz, data collected as pA, test pulse is +10 mV jump
#                  and starts at t = 5 ms)
 
file_name_list, data_notes = create_data_notes('p2', 'p2_MC_TC_summary.csv', paths.p2_files)
responses_df = cell_response_checker('p2', file_name_list, paths.p2, sf=25)


# merge with data_notes and save
responses_data_notes = pd.merge(data_notes, responses_df)
responses_data_notes.to_csv(os.path.join(paths.tables, 'p2_responses_data_notes.csv'),
                            float_format='%8.4f', index=False)


# analysis for p14 data set
# Parameters (sampling frequency = 10 kHz, data collected as nA, test pulse is -5 mV jump
#                  and starts at t = 50 ms)
 
file_name_list, data_notes = create_data_notes('p14', 'p14_MC_TC_summary.csv', paths.p14_files)
responses_df = cell_response_checker('14', file_name_list, paths.p14, sf=10)


# merge with data_notes and save
responses_data_notes = pd.merge(data_notes, responses_df)
responses_data_notes.to_csv(os.path.join(paths.tables, 'p14_responses_data_notes.csv'),
                            float_format='%8.4f', index=False)

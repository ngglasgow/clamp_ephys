# -*- coding: utf-8 -*-
"""
Created 14 November 2019
ipsc_peak_x.y.z.py

Modified from EPSC_0.1.6_meantrace.py

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

    data_notes.to_csv(os.path.join(paths.tables, '{}_data_notes.csv'.format(timepoint)))
    light_data_notes.to_csv(os.path.join(paths.tables, '{}_light_data_notes.csv'.format(timepoint)))
    spontaneous_data_notes.to_csv(os.path.join(paths.tables, '{}_spontaneous_data_notes.csv'.format(timepoint)))

    return file_name_list, data_notes

def indiv_cell_analysis(timepoint, file, data_dir, sf=25, amp_factor=1, peak_factor=-12,
                         tp_start=5, vm_jump=10, pre_tp=3):
    '''
    Make and save plots and tables for an individual file

    Parameters
    ----------
    timepoint: str
        Name of the injection timepoint used in the analysis
    file: str
        The file being analyzed
    data_dir: str
        Path to the data folder holding ibw files to be analyzed for that time point
    sf: int or float
        The sampling frequency in kHz. Default is 10 kHz.
    amp_factor: int
        Factor by which to multiply current values by to get pA. Default is 1.
    peak_factor: int or float
        Factor to multiply rs_peak by to get pA. Default is -12.
    tp_start: int or float
        Time in ms when test pulse begins.
    vm_jump: int or float
        Amplitude of whatever windows needs here the test pulse voltage command in mV.
        This is 10 mV in MIES and -5 in Nathan's Igor program.
    pre_tp: int or float
        Time in ms before start of test pulse by which to measure the baseline.

    '''
       
    # gather metadata and set some key parameters for use later on in loop
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
        Blue light comes on at 500 ms
        Peak within 150 ms of blue light
    
    '''
    baseline = mean_baseline(sf, data, 500)
    peaks = epsc_peak(sf, data, baseline, stim_time=500, post_stim=250)
    

    '''
        Pull out EPSC peaks from FILTERED signals
        Baseline 100 ms preceding blue light
        Blue light comes on at 500 ms
        Peak within 150 ms of blue light
    '''

    # filter signal with butterworth filter at 500 Hz for data
    filt_data = elephant.signal_processing.butter(data.T,
                                                lowpass_freq=500.0,
                                                fs=10000.0)
    filt_data = pd.DataFrame(filt_data).T

    filt_baseline = mean_baseline(sf=sf, data=filt_data, stim_time=500)
    filt_peaks = epsc_peak(sf, filt_data, filt_baseline, stim_time=500, post_stim=250)

    ''' Calculating Series Resistance (rs) from test pulse (tp) '''
    rs = series_resistance(sf, data, tp_start, vm_jump, pre_tp, peak_factor)

    ''' Plot EPSC peaks and Rs over time of experiemnt '''
    # set up index markers for data | drug line and drug stimuli
    # pull out number of sweeps for both conditions and all
    n_control_sweeps = len(peaks)

    # set up auto y max for peak plots (min since negative)
    y_min = peaks.min()
    y_min_lim = y_min * 1.15 * amp_factor

    # set up logic for Rs y scaling: if < 20 MOhms, don't scale, if > scale
    if rs.max() <= 20:
        rs_y_min = 0
        rs_y_max = 20
    else:
        rs_y_min = rs.min() * 0.5
        rs_y_max = rs.max() * 1.2

    # make a figure with 2 plots
    fig, axs = plt.subplots(2, 2, figsize=(6, 6), constrained_layout=True)
    fig.suptitle('Summary for ' + timepoint + ' ' + cell_type + ' ' + cell_id + ' ' + condition)

    # optional for plotting unfiltered on same graph for comparison
    axs[0, 0].plot(peaks*amp_factor, marker='.', color='darkgray', linestyle='', label='raw')

    # plot the filterd peak currents NOTE: convert peak values to pA
    axs[0, 0].plot(filt_peaks*amp_factor, color='k', marker='.', linestyle='', label='filtered')
    axs[0, 0].set_xlabel('Stimulus Number')
    axs[0, 0].set_ylabel('EPSC Peak (pA)')
    axs[0, 0].set_ylim(0, y_min_lim)
    axs[0, 0].legend()

    # plot the series resistance values
    axs[0, 1].plot(rs, marker='.', color='k', linestyle='')
    axs[0, 1].set_xlabel('Stimulus Number')
    axs[0, 1].set_ylabel('Rs (MOhm)')
    axs[0, 1].set_ylim(rs_y_min, rs_y_max)

    ''' Plot averaged EPSC trace overlaying all the individual traces '''
    # calculate the mean and the SEM of the entire time series
    filt_subtracted = filt_data - filt_baseline
    filt_data_mean = filt_subtracted.mean(axis=1)
    filt_data_sem = filt_subtracted.sem(axis=1)
    filt_data_std = filt_subtracted.std(axis=1)

    # calculate auto y min limit for mean + std
    mean_std = (filt_data_mean - filt_data_std)
    y_min_mean_std = mean_std[5000:].min()
    y_min_mean_lim = y_min_mean_std * 1.1 * amp_factor

    # set up time value for length of traces and window of what to plot
    sweep_length = len(data)                  # allow for different sweep length
    sweep_time = np.arange(0, sweep_length/sf, 1/sf)     # time of sweeps in ms

    # set up length of line for light stimulation
    blue_light = np.arange(500, 550, 1/sf)

    # plot mean data trace with all traces in gray behind
    axs[1, 0].plot(sweep_time, filt_subtracted*amp_factor, color='darkgray', linewidth=0.5)
    axs[1, 0].plot(sweep_time, filt_data_mean*amp_factor, color='k')
    axs[1, 0].hlines(75, 500, 550, color='deepskyblue')
    axs[1, 0].set_xlabel('Time (ms)')
    axs[1, 0].set_ylabel('Current (pA)')
    axs[1, 0].set_xlim(450, 1000)
    axs[1, 0].set_ylim(y_min_lim, 100)

    # plot mean data trace with shaded SEM gray behind
    axs[1, 1].plot(sweep_time, filt_data_mean*amp_factor, color='k', label='mean')
    axs[1, 1].fill_between(sweep_time,
                        (filt_data_mean - filt_data_std) * amp_factor,
                        (filt_data_mean + filt_data_std) * amp_factor,
                        color='darkgray',
                        label='st. dev.')
    axs[1, 1].hlines(75, 500, 550, color='deepskyblue')
    axs[1, 1].set_xlabel('Time (ms)')
    axs[1, 1].set_ylabel('Current (pA)')
    axs[1, 1].set_xlim(450, 1000)
    axs[1, 1].set_ylim(y_min_mean_lim, 100)
    axs[1, 1].legend(loc=1)

    # fig

    # save figure to file

    fig_save_path = os.path.join(paths.figures, timepoint, cell_type, condition, cell_id)
    fig.savefig(fig_save_path +  '_' + timepoint + '_' + cell_type + '_' + condition + '_summary.png', 
        dpi=300, format='png')
    
    plt.close()

    ''' Save all sweeps metadata for raw, filtered and rs to a csv file '''
    # save each sweep raw peak, filtered peak, and Rs with metadata to summary file
    data_dict = OrderedDict()
    data_dict['Raw Peaks (nA)'] = peaks
    data_dict['Filtered Peaks (nA)'] = filt_peaks
    data_dict['Rs (MOhms)'] = rs
    sweep_data = pd.DataFrame(data_dict)

    # fill in n=sweeps of metadata_df so that can join with peaks for clean df
    metadata_df = pd.DataFrame()
    for i in range(len(peaks)):
        metadata_df = pd.concat([metadata_df, metadata], ignore_index=True)

    # join summary data with metadata
    sweep_meta_data = metadata_df.join(sweep_data)

    # if we decide later on to combine all sweeps and mean, stdev, sem into one
    # file then the following code would be used to make an initial all sweeps
    # df that we can concat the summary measures at the bottom, would need to
    # change it a bit too. I don't think it's necessary now, OK to have 2 files
    # per cell and just open each individually
    # sweep_number = pd.DataFrame(range(len(sweep_data)), columns=['Sweep #'])
    # sweep_number

    # save summary_data to file
    sweep_meta_path = os.path.join(paths.tables, timepoint, cell_type, condition, cell_id)
    sweep_meta_data.to_csv(sweep_meta_path + '_' + timepoint + '_' + cell_type + '_' + condition + 
        '_all_sweeps_data.csv', float_format='%8.4f', index=False, header=True)


    ''' Find the mean peak epscs for raw, filtered, and rs, and save '''
    # make metadata into a series for easy appending
    metaseries = metadata.loc[0]

    # find mean, st dev, and sem of all sweeps for raw, filt, and rs
    mean = metaseries.append(sweep_data.mean())
    std = metaseries.append(sweep_data.std())
    sem = metaseries.append(sweep_data.sem())

    # combine into dataframe, add measure type string and # of sweeps
    summary_data = pd.DataFrame([mean, std, sem])
    measures = pd.DataFrame([['mean'],['st. dev.'],['sem']], columns=['Measure'])
    n_sweeps = pd.DataFrame(len(sweep_data), index=range(3), columns=['# Sweeps'])

    summary_data = pd.concat([summary_data, measures, n_sweeps], axis=1)

    # define path for saving file and save it
    summary_path = os.path.join(paths.tables, timepoint, cell_type, condition, cell_id)
    summary_data.to_csv(summary_path + '_' + timepoint + '_' + cell_type + '_' + condition + 
        '_summary_data.csv', float_format='%8.4f', index=False)


    '''###### save filtered subtracted mean time series data to file #######'''
    mean_timeseries_path = os.path.join(paths.tables, timepoint, cell_type, condition, cell_id)
    filt_data_mean.to_csv(mean_timeseries_path + '_' + timepoint + '_' + cell_type + '_' + condition + 
        '_mean_timeseries.csv', float_format='%8.4f', index=False, header=False)


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
for file in file_name_list:
    indiv_cell_analysis('p2', file, paths.p2, sf=25, amp_factor=1, peak_factor=-12,
                         tp_start=5, vm_jump=10, pre_tp=3)

# analysis for p14 data set
# p14 parameters (sampling frequency = 10 kHz, data collected as nA, test pulse is -5 mV jump
#                  and starts at t = 50 ms)

file_name_list, data_notes = create_data_notes('p14', 'p14_MC_TC_summary.csv', paths.p14_files)
for file in file_name_list:
    indiv_cell_analysis('p14', file, paths.p14, sf=10, amp_factor=1000, peak_factor=-9,
                         tp_start=50, vm_jump=-5, pre_tp=11)
